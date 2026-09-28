/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#include "vio.h"
#include <Eigen/Eigenvalues>
#include <omp.h>
#include <exception>
#include <algorithm>
#include <cmath>

namespace {
// cake_slam pattern: independent slots, exception capture inside OpenMP,
// then serial ordered state/map writes. Cameras themselves remain sequential.
template<class Function> void parallelPatches(int count,int threads,Function&& function)
{
  if(count<=0) return;
  threads=std::min(threads,count);
  const int chunk=count<4*threads ? 1 : 4;
  std::vector<std::exception_ptr> failures(count);
  #pragma omp parallel for num_threads(threads) schedule(dynamic,chunk) if(threads>1)
  for(int i=0;i<count;++i) {
    try { function(i); }
    catch(...) { failures[i]=std::current_exception(); }
  }
  for(const auto& failure:failures) if(failure) std::rethrow_exception(failure);
}
// Only bandwidth changes between pyramid levels at a fixed pose. Reuse all
// ray-plane intersections, camera lookups and angular warp differentials.
void scaleWarpGeometry(const SphericalWarp& source,int level_delta,SphericalWarp& target)
{
  const double scale=std::ldexp(1.0,level_delta);
  target.world=source.world;target.bearings=source.bearings;
  target.reference_radii=source.reference_radii;target.covariance=source.covariance;
  for(double& radius:target.reference_radii) radius*=scale;
  for(auto& covariance:target.covariance) covariance*=scale*scale;
  target.pose_R=source.pose_R;target.pose_t=source.pose_t;
  // Sensor-pyramid footprints are configured separately at every scale.
  target.cover=0;target.geometry_ready=false;target.ready=false;
}
int depthCell(double coordinate) {
  return std::max(0,std::min(20,int(std::floor((coordinate+1.0)*10.0))));
}
int depthKey(int x,int y,int z) { return (z*21+y)*21+x; }
}

VIOManager::VIOManager()
{
  visual_threads=std::max(1,std::min(8,omp_get_max_threads()));
}

VIOManager::~VIOManager()
{
  delete visual_submap;
  for (auto& pair : feat_map) delete pair.second;
  feat_map.clear();
}

void VIOManager::setImuToLidarExtrinsic(const V3D &transl, const M3D &rot)
{
  Pli = -rot.transpose() * transl;
  Rli = rot.transpose();
}

void VIOManager::setLidarToCameraExtrinsic(vector<double> &R, vector<double> &P)
{
  Rcl << MAT_FROM_ARRAY(R);
  Pcl << VEC_FROM_ARRAY(P);
}

void VIOManager::initializeVIO()
{
  if (!visual_submap) visual_submap = new SubSparseMap;

  fx = cam->fx();
  fy = cam->fy();
  cx = cam->cx();
  cy = cam->cy();
  image_resize_factor = cam->scale();

  width = cam->width();
  height = cam->height();

  Rci = Rcl * Rli;
  Pci = Rcl * Pli + Pcl;

  V3D Pic;
  M3D tmp;
  Jdphi_dR = Rci;
  Pic = -Rci.transpose() * Pci;
  tmp << SKEW_SYM_MATRX(Pic);
  Jdp_dR = -Rci * tmp;

  patch_size_total = sample_count;
  if(sphere_template.empty()) sphere_template=spherical::capTemplate(sample_count,patch_radius);
  patch_size_half = static_cast<int>(patch_size / 2);
  patch_buffer.resize(patch_size_total);
  warp_len = patch_size_total * patch_pyrimid_level;
  border = 5; // Shi-Tomasi border; full spherical footprint checked per sample.

  sub_feat_map.clear();
}

void VIOManager::resetGrid()
{
  // Keep the upstream per-frame reset entry point; selection is global now.
  update_flag.assign(max_patches,0);
  visible_map_positions.clear();
  total_points=0;
  candidate_count=selected_candidate_count=0;
}

// void VIOManager::resetRvizDisplay()
// {
  // sub_map_ray.clear();
  // sub_map_ray_fov.clear();
  // visual_sub_map_cur.clear();
  // visual_converged_point.clear();
  // map_cur_frame.clear();
  // sample_points.clear();
// }

void VIOManager::insertPointIntoVoxelMap(VisualPoint *pt_new)
{
  V3D pt_w(pt_new->pos_[0], pt_new->pos_[1], pt_new->pos_[2]);
  double voxel_size = 0.5;
  float loc_xyz[3];
  for (int j = 0; j < 3; j++)
  {
    loc_xyz[j] = pt_w[j] / voxel_size;
    if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
  }
  VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);
  auto iter = feat_map.find(position);
  if (iter != feat_map.end())
  {
    iter->second->voxel_points.push_back(pt_new);
    iter->second->count++;
  }
  else
  {
    VOXEL_POINTS *ot = new VOXEL_POINTS(0);
    ot->voxel_points.push_back(pt_new);
    feat_map[position] = ot;
  }
}

void VIOManager::retrieveFromVisualSparseMap(cv::Mat img, vector<pointWithVar> &pg, const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map)
{
  visual_submap->reset();
  if(feat_map.empty()) return;
  // Retain upstream local-map discovery. Pixel cells do not discard candidates.
  sub_feat_map.clear();
  for(const auto& point:pg) {
    if(!point.point_w.allFinite()) continue;
    int64_t key[3];
    for(int axis=0;axis<3;++axis) {
      double coordinate=point.point_w[axis]/0.5;
      if(coordinate<0) coordinate-=1.0;
      key[axis]=static_cast<int64_t>(coordinate);
    }
    sub_feat_map[VOXEL_LOCATION(key[0],key[1],key[2])]=0;
  }
  struct Candidate {
    std::array<double,3> position;
    float score;
    VisualPoint* point;
    Feature* reference;
    std::vector<SphericalWarp> warps;
    float error=0;
    bool accepted=false;
  };
  std::vector<Candidate> candidates;
  for(const auto& local:sub_feat_map) {
    const auto voxel=feat_map.find(local.first);
    if(voxel==feat_map.end()) continue;
    for(auto* point:voxel->second->voxel_points) {
      if(!point || !point->pos_.allFinite()) continue;
      const V2D pixel=new_frame_->w2c(point->pos_);
      if(!pixel.allFinite() || !cam->isInFrame(pixel.cast<int>(),border)) continue;
      const std::array<double,3> position{point->pos_.x(),point->pos_.y(),point->pos_.z()};
      // Exact duplicate suppression only: no one-per-grid or angular quota.
      if(!visible_map_positions.insert(position).second) continue;
      if(point->obs_.empty() || !refreshPlane(*point)) continue;
      Feature* reference=nullptr;
      if(!point->getCloseViewObs(new_frame_->pos(),reference,pixel) || !reference->sphere_image) continue;
      const float score=vk::shiTomasiScore(img,int(pixel.x()),int(pixel.y()));
      if(!std::isfinite(score)) continue;
      candidates.push_back({position,score,point,reference,{}});
    }
  }
  candidate_count=int(candidates.size());
  const int keep=std::min(max_patches,candidate_count);
  if(keep<candidate_count) {
    std::partial_sort(candidates.begin(),candidates.begin()+keep,candidates.end(),
      [](const Candidate& a,const Candidate& b) {
        return a.score!=b.score ? a.score>b.score : a.position<b.position;
      });
    candidates.resize(keep);
  }
  // Deterministic accumulation independent of unordered voxel traversal.
  std::sort(candidates.begin(),candidates.end(),
    [](const Candidate& a,const Candidate& b) { return a.position<b.position; });
  selected_candidate_count=keep;
  parallelPatches(keep,visual_threads,[&](int j) {
    auto& candidate=candidates[j];
    if(!prepareSphericalPyramid(*candidate.reference,candidate.warps)) return;
    const auto& fine=candidate.warps[0];
    for(int k=0;k<patch_size_total;++k) {
      const double residual=255*state->inv_expo_time*fine.current[k].value-
        candidate.reference->inv_expo_time_*fine.reference[k];
      candidate.error+=residual*residual;
    }
    candidate.accepted=candidate.error<=outlier_threshold*patch_size_total;
  });
  for(auto& candidate:candidates) {
    if(!candidate.accepted) continue;
    visual_submap->voxel_points.push_back(candidate.point);
    visual_submap->references.push_back(candidate.reference);
    visual_submap->sphere_warps.push_back(std::move(candidate.warps));
    visual_submap->propa_errors.push_back(candidate.error);
    visual_submap->errors.push_back(candidate.error);
    visual_submap->search_levels.push_back(0);
    visual_submap->inv_expo_list.push_back(candidate.reference->inv_expo_time_);
  }
  total_points=int(visual_submap->voxel_points.size());
  printf("[ VIO ] Retrieve %d points from visual sparse map\n",total_points);
}

void VIOManager::computeJacobianAndUpdateEKF(cv::Mat img)
{
  compute_jacobian_time = update_ekf_time = prepare_patches_time = 0;
  G.setZero();
  if (total_points == 0) return;
  
  compute_jacobian_time = update_ekf_time = 0.0;

  for (int level = patch_pyrimid_level - 1; level >= 0; level--)
  {
    updateState(img, level);
  }
  state->cov -= G * state->cov;
  updateFrameState(*state);
}

void VIOManager::generateVisualMapPoints(cv::Mat img, vector<pointWithVar> &pg)
{
  const int budget=std::max(0,max_patches-total_points);
  if(pg.size()<=10 || budget==0) return;
  struct Candidate {
    const pointWithVar* point;
    V2D pixel;
    float score;
    std::array<double,3> position;
  };
  std::vector<Candidate> candidates;
  candidates.reserve(pg.size());
  auto positions=visible_map_positions;
  const auto collect=[&](const pointWithVar& point) {
    if(!point.point_w.allFinite() || !point.normal.allFinite() || point.normal.squaredNorm()==0) return;
    const V2D pixel=new_frame_->w2c(point.point_w);
    if(!pixel.allFinite() || !cam->isInFrame(pixel.cast<int>(),border)) return;
    const std::array<double,3> position{point.point_w.x(),point.point_w.y(),point.point_w.z()};
    if(!positions.insert(position).second) return;
    const float score=vk::shiTomasiScore(img,int(pixel.x()),int(pixel.y()));
    if(std::isfinite(score) && score>0) candidates.push_back({&point,pixel,score,position});
  };
  for(const auto& point:pg) collect(point);
  for(const auto& point:visual_submap->add_from_voxel_map) collect(point);
  const int keep=std::min(budget,int(candidates.size()));
  if(keep<int(candidates.size())) {
    std::partial_sort(candidates.begin(),candidates.begin()+keep,candidates.end(),
      [](const Candidate& a,const Candidate& b) {
        return a.score!=b.score ? a.score>b.score : a.position<b.position;
      });
    candidates.resize(keep);
  }
  std::sort(candidates.begin(),candidates.end(),
    [](const Candidate& a,const Candidate& b) { return a.position<b.position; });
  std::vector<std::unique_ptr<VisualPoint>> points(keep);
  parallelPatches(keep,visual_threads,[&](int j) {
    const auto& pt_var=*candidates[j].point;
    auto point=std::make_unique<VisualPoint>(pt_var.point_w);
    point->covariance_=pt_var.var;point->is_normal_initialized_=true;
    const V3D pc3=new_frame_->w2f(pt_var.point_w);
    point->normal_=pt_var.normal;
    if(pc3.dot(new_frame_->T_f_w_.rotationMatrix()*pt_var.normal)<0) point->normal_=-point->normal_;
    point->previous_normal_=point->normal_;
    if(!refreshPlane(*point)) return; // exclusively owned candidate
    const V2D pc=candidates[j].pixel;
    auto patch=std::make_unique<float[]>(patch_size_total);
    auto feature=std::make_unique<Feature>(point.get(),patch.get(),pc,cam->cam2world(pc),new_frame_->T_f_w_,0);
    patch.release();
    feature->img_=img;feature->id_=new_frame_->id_;feature->inv_expo_time_=state->inv_expo_time;
    if(!buildSphericalReference(*feature)) return;
    point->addFrameRef(feature.get());feature.release();
    points[j]=std::move(point);
  });
  int add=0;
  for(auto& point:points) if(point) { insertPointIntoVoxelMap(point.get());point.release();++add; }
  printf("[ VIO ] Append %d new visual map points\n",add);
}

void VIOManager::updateVisualMapPoints(cv::Mat img)
{
  if(total_points==0) return;
  struct PendingReference { int index; V2D pixel; std::unique_ptr<Feature> feature; };
  std::vector<PendingReference> pending;
  const SE3 pose_cur=new_frame_->T_f_w_;
  for(int i=0;i<total_points;++i) {
    auto* pt=visual_submap->voxel_points[i];
    if(!pt) continue;
    if(pt->is_converged_) { pt->deleteNonRefPatchFeatures();continue; }
    const V2D pc=new_frame_->w2c(pt->pos_);
    auto* previous=pt->obs_.front();
    const SE3 delta_pose=previous->T_f_w_*pose_cur.inverse();
    const double trace=delta_pose.rotationMatrix().trace();
    const double angle=trace>3.0-1e-6 ? 0.0 : std::acos(std::max(-1.0,std::min(1.0,0.5*(trace-1))));
    const bool add=delta_pose.translation().norm()>0.5 || angle>0.3 ||
      previous->camera_id!=active_camera || (pc-previous->px_).norm()>40;
    if(pt->obs_.size()>=30) {
      Feature* remove=nullptr;pt->findMinScoreFeature(new_frame_->pos(),remove);pt->deleteFeatureRef(remove);
    }
    if(add) pending.push_back({i,pc,nullptr});
  }
  parallelPatches(int(pending.size()),visual_threads,[&](int j) {
    auto& item=pending[j];
    auto* pt=visual_submap->voxel_points[item.index];
    auto patch=std::make_unique<float[]>(patch_size_total);
    auto feature=std::make_unique<Feature>(pt,patch.get(),item.pixel,cam->cam2world(item.pixel),
      new_frame_->T_f_w_,visual_submap->search_levels[item.index]);
    patch.release();
    feature->img_=img;feature->id_=new_frame_->id_;feature->inv_expo_time_=state->inv_expo_time;
    if(buildSphericalReference(*feature)) item.feature=std::move(feature);
  });
  int update_num=0;
  for(auto& item:pending) if(item.feature) {
    visual_submap->voxel_points[item.index]->addFrameRef(item.feature.get());item.feature.release();
    update_flag[item.index]=1;++update_num;
  }
  printf("[ VIO ] Update %d points in visual submap\n",update_num);
}

void VIOManager::updateReferencePatch(const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &plane_map)
{
  if (total_points == 0) return;

  for (int i = 0; i < visual_submap->voxel_points.size(); i++)
  {
    VisualPoint *pt = visual_submap->voxel_points[i];

    if (!pt->is_normal_initialized_) continue;
    if (pt->is_converged_) continue;
    if (pt->obs_.size() <= 5) continue;
    if (update_flag[i] == 0) continue;

    const V3D &p_w = pt->pos_;
    float loc_xyz[3];
    for (int j = 0; j < 3; j++)
    {
      loc_xyz[j] = p_w[j] / map_voxel_size;
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
    }
    VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);
    auto iter = plane_map.find(position);
    if (iter != plane_map.end())
    {
      VoxelOctoTree *current_octo;
      current_octo = iter->second->find_correspond(p_w);
      if (current_octo && current_octo->plane_ptr_->is_plane_)
      {
        VoxelPlane &plane = *current_octo->plane_ptr_;
        float dis_to_plane = plane.normal_(0) * p_w(0) + plane.normal_(1) * p_w(1) + plane.normal_(2) * p_w(2) + plane.d_;
        float dis_to_plane_abs = fabs(dis_to_plane);
        float dis_to_center = (plane.center_(0) - p_w(0)) * (plane.center_(0) - p_w(0)) +
                              (plane.center_(1) - p_w(1)) * (plane.center_(1) - p_w(1)) + (plane.center_(2) - p_w(2)) * (plane.center_(2) - p_w(2));
        float range_dis = sqrt(dis_to_center - dis_to_plane * dis_to_plane);
        if (range_dis <= 3 * plane.radius_)
        {
          Eigen::Matrix<double, 1, 6> J_nq;
          J_nq.block<1, 3>(0, 0) = p_w - plane.center_;
          J_nq.block<1, 3>(0, 3) = -plane.normal_;
          double sigma_l = J_nq * plane.plane_var_ * J_nq.transpose();
          sigma_l += plane.normal_.transpose() * pt->covariance_ * plane.normal_;

          if (dis_to_plane_abs < 3 * sqrt(sigma_l))
          {
            // V3D norm_vec(new_frame_->T_f_w_.rotationMatrix() * plane.normal_);
            // V3D pf(new_frame_->T_f_w_ * pt->pos_);
            // V3D pf_ref(pt->ref_patch->T_f_w_ * pt->pos_);
            // V3D norm_vec_ref(pt->ref_patch->T_f_w_.rotationMatrix() *
            // plane.normal); double cos_ref = pf_ref.dot(norm_vec_ref);
            
            if (pt->previous_normal_.dot(plane.normal_) < 0) { pt->normal_ = -plane.normal_; }
            else { pt->normal_ = plane.normal_; }

            double normal_update = (pt->normal_ - pt->previous_normal_).norm();

            pt->previous_normal_ = pt->normal_;

            if (normal_update < 0.0001 && pt->obs_.size() > 10)
            {
              pt->is_converged_ = true;
              // visual_converged_point.push_back(pt);
            }
          }
        }
      }
    }

    Feature* closest=nullptr;
    if(pt->getCloseViewObs(new_frame_->pos(),closest,new_frame_->w2c(pt->pos_))) {
      pt->ref_patch=closest; pt->has_ref_patch_=true;
    }
  }
}

void VIOManager::updateState(cv::Mat img, int level)
{
  if (total_points == 0) return;
  StatesGroup old_state = (*state);
  const double prepare_start=omp_get_wtime();
  updateFrameState(*state);
  std::vector<unsigned char> valid(total_points,0);
  parallelPatches(total_points,visual_threads,[&](int i) {
    auto& warp=visual_submap->sphere_warps[i][level];
    // Retrieval and the first coarse iteration have exactly the same pose.
    // Reuse samples/Jacobian inputs only with an exact pose match; any update
    // requires freshly prepared bandwidth and current-image samples.
    if(warp.ready && (warp.pose_R-Rcw).squaredNorm()==0 && (warp.pose_t-Pcw).squaredNorm()==0) {
      valid[i]=1;return;
    }
    valid[i]=prepareSphericalWarp(*visual_submap->references[i],level,warp);
  });
  prepare_patches_time+=omp_get_wtime()-prepare_start;
  for(unsigned char ready:valid) if(!ready) return;

  VectorXd z;
  MatrixXd H_sub;
  bool EKF_end = false;
  float last_error = std::numeric_limits<float>::max();

  const int H_DIM = total_points * patch_size_total;
  z.resize(H_DIM);
  z.setZero();
  H_sub.resize(H_DIM, 7);
  H_sub.setZero();
  struct PatchResult { float error=0; int measurements=0; bool invalid=false; };
  std::vector<PatchResult> results(total_points);

  for (int iteration = 0; iteration < max_iterations; iteration++)
  {
    double t1 = omp_get_wtime();

    M3D Rwi(state->rot_end);
    V3D Pwi(state->pos_end);
    Rcw = Rci * Rwi.transpose();
    Pcw = -Rci * Rwi.transpose() * Pwi + Pci;
    Jdp_dt = Rci * Rwi.transpose();
    
    float error = 0.0;
    int n_meas = 0;
    int invalid=0;
    H_sub.setZero(); z.setZero();
    parallelPatches(total_points,visual_threads,[&](int i)
    {
      auto& result=results[i];result=PatchResult{};
      const auto& warp=visual_submap->sphere_warps[i][level];
      const double inv_ref_expo=visual_submap->inv_expo_list[i];
      spherical::Sample sampled;
      for(int k=0;k<patch_size_total;++k) {
        const V3D pi=Rwi.transpose()*(warp.world[k]-Pwi),pc=Rci*pi+Pci;
        const double range=pc.norm();
        if(range<1e-6) { result.invalid=true;break; }
        const spherical::Sample* sample=&warp.current[k];
        if(iteration!=0) {
          if(!sphere_image->sample(pc/range,warp.current_filters[k],sampled)) { result.invalid=true;break; }
          sample=&sampled;
        }
        const double cur_value=255*sample->value;
        const double res=state->inv_expo_time*cur_value-inv_ref_expo*warp.reference[k];
        const int row=i*patch_size_total+k;
        z(row)=res;
        // sample.gradient is already tangent to pc/range. Reuse pi/range and
        // Rwi/Rci instead of normalizing and transforming the same point twice.
        const Eigen::RowVector3d j=(255*state->inv_expo_time/range)*sample->gradient.transpose()*Rci;
        H_sub.block<1,3>(row,0)=j*spherical::skew(pi);
        H_sub.block<1,3>(row,3)=-j*Rwi.transpose();
        if(exposure_estimate_en) H_sub(row,6)=cur_value;
        result.error+=res*res;++result.measurements;
      }
    });
    // Fixed summation order, as in cake_slam; no worker writes shared state.
    for(int i=0;i<total_points;++i) {
      visual_submap->errors[i]=results[i].error;
      error+=results[i].error;n_meas+=results[i].measurements;
      invalid+=results[i].invalid;
    }
    // Keep the same support during an iteration; losing samples is not an
    // improvement in photometric cost.
    if(invalid || n_meas==0) { *state=old_state; break; }

    error = error / n_meas;
    
    compute_jacobian_time += omp_get_wtime() - t1;

    // printf("\nPYRAMID LEVEL %i\n---------------\n", level);
    // std::cout << "It. " << iteration
    //           << "\t last_error = " << last_error
    //           << "\t new_error = " << error
    //           << std::endl;

    double t3 = omp_get_wtime();

    if (error <= last_error)
    {
      old_state = (*state);
      last_error = error;

      // K = (H.transpose() / img_point_cov * H + state->cov.inverse()).inverse() * H.transpose() / img_point_cov; auto
      // vec = (*state_propagat) - (*state); G = K*H;
      // (*state) += (-K*z + vec - G*vec);

      auto &&H_sub_T = H_sub.transpose();
      H_T_H.setZero();
      G.setZero();
      H_T_H.block<7, 7>(0, 0) = H_sub_T * H_sub;
      MD(DIM_STATE, DIM_STATE) &&K_1 = (H_T_H + (state->cov / img_point_cov).inverse()).inverse();
      auto &&HTz = H_sub_T * z;
      // K = K_1.block<DIM_STATE,6>(0,0) * H_sub_T;
      auto vec = (*state_propagat) - (*state);
      G.block<DIM_STATE, 7>(0, 0) = K_1.block<DIM_STATE, 7>(0, 0) * H_T_H.block<7, 7>(0, 0);
      MD(DIM_STATE, 1)
      solution = -K_1.block<DIM_STATE, 7>(0, 0) * HTz + vec - G.block<DIM_STATE, 7>(0, 0) * vec.block<7, 1>(0, 0);

      if(!solution.allFinite()) { *state=old_state; break; }
      (*state) += solution;
      if(!(state->inv_expo_time>0) || !std::isfinite(state->inv_expo_time)) { *state=old_state; break; }
      auto &&rot_add = solution.block<3, 1>(0, 0);
      auto &&t_add = solution.block<3, 1>(3, 0);

      auto &&expo_add = solution.block<1, 1>(6, 0);
      // if ((rot_add.norm() * 57.3f < 0.001f) && (t_add.norm() * 100.0f < 0.001f) && (expo_add.norm() < 0.001f)) EKF_end = true;
      if ((rot_add.norm() * 57.3f < 0.001f) && (t_add.norm() * 100.0f < 0.001f))  EKF_end = true;
    }
    else
    {
      (*state) = old_state;
      EKF_end = true;
    }

    update_ekf_time += omp_get_wtime() - t3;

    if (iteration == max_iterations || EKF_end) break;
  }
  // if (state->inv_expo_time < 0.0)  {ROS_ERROR("reset expo time!!!!!!!!!!\n"); state->inv_expo_time = 0.0;}
}

void VIOManager::updateFrameState(StatesGroup state)
{
  M3D Rwi(state.rot_end);
  V3D Pwi(state.pos_end);
  Rcw = Rci * Rwi.transpose();
  Pcw = -Rci * Rwi.transpose() * Pwi + Pci;
  new_frame_->T_f_w_ = SE3(Rcw, Pcw);
  updateDepthIndex();
}

void VIOManager::plotTrackedPoints()
{
  int total_points = visual_submap->voxel_points.size();
  if (total_points == 0) return;
  // int inlier_count = 0;
  // for (int i = 0; i < img_cp.rows / grid_size; i++)
  // {
  //   cv::line(img_cp, cv::Poaint2f(0, grid_size * i), cv::Point2f(img_cp.cols, grid_size * i), cv::Scalar(255, 255, 255), 1, CV_AA);
  // }
  // for (int i = 0; i < img_cp.cols / grid_size; i++)
  // {
  //   cv::line(img_cp, cv::Point2f(grid_size * i, 0), cv::Point2f(grid_size * i, img_cp.rows), cv::Scalar(255, 255, 255), 1, CV_AA);
  // }
  // for (int i = 0; i < img_cp.rows / grid_size; i++)
  // {
  //   cv::line(img_cp, cv::Point2f(0, grid_size * i), cv::Point2f(img_cp.cols, grid_size * i), cv::Scalar(255, 255, 255), 1, CV_AA);
  // }
  // for (int i = 0; i < img_cp.cols / grid_size; i++)
  // {
  //   cv::line(img_cp, cv::Point2f(grid_size * i, 0), cv::Point2f(grid_size * i, img_cp.rows), cv::Scalar(255, 255, 255), 1, CV_AA);
  // }
  for (int i = 0; i < total_points; i++)
  {
    VisualPoint *pt = visual_submap->voxel_points[i];
    V2D pc(new_frame_->w2c(pt->pos_));
    for(const V3D& world:visual_submap->sphere_warps[i][0].world) {
      V2D uv=new_frame_->w2c(world);
      if(cam->isInFrame(uv.cast<int>(),1)) cv::circle(img_cp,cv::Point2f(uv.x(),uv.y()),1,cv::Scalar(0,255,255),-1);
    }

    if (visual_submap->errors[i] <= visual_submap->propa_errors[i])
    {
      // inlier_count++;
      cv::circle(img_cp, cv::Point2f(pc[0], pc[1]), 7, cv::Scalar(0, 255, 0), -1, 8); // Green Sparse Align tracked
    }
    else
    {
      cv::circle(img_cp, cv::Point2f(pc[0], pc[1]), 7, cv::Scalar(255, 0, 0), -1, 8); // Blue Sparse Align tracked
    }
  }
  // std::string text = std::to_string(inlier_count) + " " + std::to_string(total_points);
  // cv::Point2f origin;
  // origin.x = img_cp.cols - 110;
  // origin.y = 20;
  // cv::putText(img_cp, text, origin, cv::FONT_HERSHEY_COMPLEX, 0.7, cv::Scalar(0, 255, 0), 2, 8, 0);
}

V3F VIOManager::getInterpolatedPixel(cv::Mat img, V2D pc)
{
  const float u_ref = pc[0];
  const float v_ref = pc[1];
  const int u_ref_i = floorf(pc[0]);
  const int v_ref_i = floorf(pc[1]);
  const float subpix_u_ref = (u_ref - u_ref_i);
  const float subpix_v_ref = (v_ref - v_ref_i);
  const float w_ref_tl = (1.0 - subpix_u_ref) * (1.0 - subpix_v_ref);
  const float w_ref_tr = subpix_u_ref * (1.0 - subpix_v_ref);
  const float w_ref_bl = (1.0 - subpix_u_ref) * subpix_v_ref;
  const float w_ref_br = subpix_u_ref * subpix_v_ref;
  const size_t stride=img.step;
  const uint8_t *img_ptr = img.ptr<uint8_t>(v_ref_i) + u_ref_i*3;
  float B = w_ref_tl * img_ptr[0] + w_ref_tr * img_ptr[0 + 3] + w_ref_bl * img_ptr[stride] + w_ref_br * img_ptr[stride + 0 + 3];
  float G = w_ref_tl * img_ptr[1] + w_ref_tr * img_ptr[1 + 3] + w_ref_bl * img_ptr[1 + stride] + w_ref_br * img_ptr[stride + 1 + 3];
  float R = w_ref_tl * img_ptr[2] + w_ref_tr * img_ptr[2 + 3] + w_ref_bl * img_ptr[2 + stride] + w_ref_br * img_ptr[stride + 2 + 3];
  V3F pixel(B, G, R);
  return pixel;
}

void VIOManager::processFrame(cv::Mat &img, vector<pointWithVar> &pg, const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &feat_map, double img_time)
{
  const double frame_start=omp_get_wtime();
  if(img.empty() || width!=img.cols || height!=img.rows)
    throw std::runtime_error("raw image dimensions differ from camera calibration");
  current_points=&pg; current_planes=&feat_map; depth_index_valid=false;
  img_rgb = img;
  img_cp = img.clone();
  // img_test = img.clone();

  if (img.channels() == 3) cv::cvtColor(img, img, CV_BGR2GRAY);

  sphere_image=std::make_shared<spherical::Image>(cameras_[active_camera].camera->atlas,img);
  new_frame_.reset(new Frame(cam, img));
  updateFrameState(*state);
  
  resetGrid();

  double t1 = omp_get_wtime();

  retrieveFromVisualSparseMap(img, pg, feat_map);

  double t2 = omp_get_wtime();

  computeJacobianAndUpdateEKF(img);

  double t3 = omp_get_wtime();

  generateVisualMapPoints(img, pg);

  double t4 = omp_get_wtime();
  
  plotTrackedPoints();

  double t5 = omp_get_wtime();

  updateVisualMapPoints(img);

  double t6 = omp_get_wtime();

  updateReferencePatch(feat_map);

  double t7 = omp_get_wtime();
  
  cameras_[active_camera].img_rgb=img_rgb;
  cameras_[active_camera].img_cp=img_cp;

  frame_count++;
  ave_total = ave_total * (frame_count - 1) / frame_count + (t7 - frame_start) / frame_count;

  // printf("[ VIO ] feat_map.size(): %zu\n", feat_map.size());
  // printf("\033[1;32m[ VIO time ]: current frame: retrieveFromVisualSparseMap time: %.6lf secs.\033[0m\n", t2 - t1);
  // printf("\033[1;32m[ VIO time ]: current frame: computeJacobianAndUpdateEKF time: %.6lf secs, comp H: %.6lf secs, ekf: %.6lf secs.\033[0m\n", t3 - t2, computeH, ekf_time);
  // printf("\033[1;32m[ VIO time ]: current frame: generateVisualMapPoints time: %.6lf secs.\033[0m\n", t4 - t3);
  // printf("\033[1;32m[ VIO time ]: current frame: updateVisualMapPoints time: %.6lf secs.\033[0m\n", t6 - t5);
  // printf("\033[1;32m[ VIO time ]: current frame: updateReferencePatch time: %.6lf secs.\033[0m\n", t7 - t6);
  // printf("\033[1;32m[ VIO time ]: current total time: %.6lf, average total time: %.6lf secs.\033[0m\n", t7 - t1 - (t5 - t4), ave_total);

  // ave_build_residual_time = ave_build_residual_time * (frame_count - 1) / frame_count + (t2 - t1) / frame_count;
  // ave_ekf_time = ave_ekf_time * (frame_count - 1) / frame_count + (t3 - t2) / frame_count;
 
  // cout << BLUE << "ave_build_residual_time: " << ave_build_residual_time << RESET << endl;
  // cout << BLUE << "ave_ekf_time: " << ave_ekf_time << RESET << endl;
  
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m|                         VIO Time                            |\033[0m\n");
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m| %-29s | %-27zu |\033[0m\n", "Sparse Map Size", feat_map.size());
  printf("\033[1;34m| %-29s | %-27d |\033[0m\n", "Camera ID", active_camera);
  printf("\033[1;34m| %-29s | %-27d |\033[0m\n", "Eligible Candidates", candidate_count);
  printf("\033[1;34m| %-29s | %-27d |\033[0m\n", "Selected Candidates", selected_candidate_count);
  printf("\033[1;34m| %-29s | %-27d |\033[0m\n", "Accepted Patches", total_points);
  printf("\033[1;34m| %-29s | %-27d |\033[0m\n", "Visual Workers", visual_threads);
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m| %-29s | %-27s |\033[0m\n", "Algorithm Stage", "Time (secs)");
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "image / depth index", t1-frame_start);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "retrieveFromVisualSparseMap", t2 - t1);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "computeJacobianAndUpdateEKF", t3 - t2);
  printf("\033[1;32m| %-27s   | %-27lf |\033[0m\n", "-> preparePatches", prepare_patches_time);
  printf("\033[1;32m| %-27s   | %-27lf |\033[0m\n", "-> computeJacobian", compute_jacobian_time);
  printf("\033[1;32m| %-27s   | %-27lf |\033[0m\n", "-> updateEKF", update_ekf_time);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "generateVisualMapPoints", t4 - t3);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "draw spherical samples", t5-t4);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "updateVisualMapPoints", t6 - t5);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "updateReferencePatch", t7 - t6);
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "Current Total Time", t7 - frame_start);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "Average Total Time", ave_total);
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");

  // std::string text = std::to_string(int(1 / (t7 - t1 - (t5 - t4)))) + " HZ";
  // cv::Point2f origin;
  // origin.x = 20;
  // origin.y = 20;
  // cv::putText(img_cp, text, origin, cv::FONT_HERSHEY_COMPLEX, 0.6, cv::Scalar(255, 255, 255), 1, 8, 0);
  // cv::imwrite("/home/chunran/Desktop/raycasting/" + std::to_string(new_frame_->id_) + ".png", img_cp);
}

// Native spherical geometry integrated into the upstream VIOManager.
// The shared VisualPoint/Feature map, selection, reference lifecycle and IEKF
// remain in vio.cpp; this replaces planar patch construction and affine warping.
void VIOManager::activateCamera(int id)
{
  active_camera=id;
  auto& view=cameras_.at(id);
  cam=view.camera.get(); Rcl=view.Rcl; Pcl=view.Pcl;
  if(sample_count<4 || max_patches<1 || !(patch_radius>0 && patch_radius<0.5))
    throw std::invalid_argument("invalid spherical patch support");
  initializeVIO();
}

bool VIOManager::getColorFromCamera(int id,const V3D& world,V3F& color,double blind)
{
  const auto& view=cameras_.at(id);
  if(view.img_rgb.empty()) return false;
  const V3D pc=view.color_Rcw*world+view.color_Pcw;
  if(pc.norm()<=blind) return false;
  const V2D uv=view.camera->world2cam(pc);
  if(!view.camera->isInFrame(uv.cast<int>(),3)) return false;
  color=getInterpolatedPixel(view.img_rgb,uv);
  return true;
}

bool VIOManager::refreshPlane(VisualPoint& pt)
{
  if(!current_planes) return false;
  // Use the same voxel indexing convention as the upstream LiDAR map.
  int64_t index[3];
  for(int j=0;j<3;++j) {
    double x=pt.pos_[j]/map_voxel_size;
    if(x<0) x-=1.0;
    index[j]=static_cast<int64_t>(x);
  }
  auto it=current_planes->find(VOXEL_LOCATION(index[0],index[1],index[2]));
  if(it==current_planes->end()) return false;
  const VoxelOctoTree* cell=it->second;
  while(cell && !cell->plane_ptr_->is_plane_ && cell->octo_state_!=0) {
    const int child=(pt.pos_.x()>cell->voxel_center_[0]?4:0)+
      (pt.pos_.y()>cell->voxel_center_[1]?2:0)+(pt.pos_.z()>cell->voxel_center_[2]?1:0);
    cell=cell->leaves_[child];
  }
  if(!cell || !cell->plane_ptr_->is_plane_) return false;
  const auto& plane=*cell->plane_ptr_;
  if(plane.normal_.squaredNorm()<0.9 || plane.radius_<=0) return false;
  const V3D delta=pt.pos_-plane.center_;
  const double distance=plane.normal_.dot(delta);
  Eigen::Matrix<double,1,6> J;
  J << delta.transpose(), -plane.normal_.transpose();
  const double variance=(J*plane.plane_var_*J.transpose())(0,0)+
    plane.normal_.dot(pt.covariance_*plane.normal_)+range_noise*range_noise;
  if(distance*distance>9*std::max(variance,1e-12) ||
     (delta-distance*plane.normal_).norm()>3*plane.radius_) return false;
  pt.plane_center=plane.center_; pt.plane_radius=plane.radius_;
  pt.plane_cov=plane.plane_var_; pt.roughness=std::max(0.0,double(plane.min_eigen_value_));
  const bool flip=pt.normal_.dot(plane.normal_)<0;
  pt.normal_=plane.normal_;
  if(flip) pt.normal_=-pt.normal_;
  if(flip) { pt.plane_cov.topRightCorner<3,3>()*=-1; pt.plane_cov.bottomLeftCorner<3,3>()*=-1; }
  pt.is_normal_initialized_=true;
  return true;
}

bool VIOManager::buildSphericalReference(Feature& feature)
{
  feature.camera_id=active_camera;
  feature.sphere_image=sphere_image;
  const V3D pc=feature.T_f_w_*feature.point_->pos_;
  if(pc.squaredNorm()<1e-12) return false;
  const M3D rotation=Eigen::Quaterniond::FromTwoVectors(V3D::UnitZ(),pc.normalized()).toRotationMatrix();
  feature.rays.clear();
  feature.ray_pitch.clear();feature.ray_basis.clear();
  feature.rays.reserve(sample_count);feature.ray_pitch.reserve(sample_count);feature.ray_basis.reserve(sample_count);
  for(const auto& node:sphere_template) {
    const V3D b=rotation*node;
    V2D uv;double pitch;
    if(!sphere_image->locate(b,uv,pitch)) return false;
    feature.rays.push_back(b);feature.ray_pitch.push_back(pitch);
    feature.ray_basis.push_back(spherical::tangentBasis(b));
  }
  std::vector<SphericalWarp> warps;
  if(!prepareSphericalPyramid(feature,warps)) return false;
  std::copy(warps[0].reference.begin(),warps[0].reference.end(),feature.patch_);
  return true;
}

bool VIOManager::prepareSphericalGeometry(const Feature& feature,int level,SphericalWarp& warp) const
{
  using namespace spherical;
  warp.ready=false;warp.geometry_ready=false;warp.cover=0;
  const auto& pt=*feature.point_;
  if(!feature.sphere_image || feature.rays.size()!=size_t(sample_count) || feature.ray_pitch.size()!=feature.rays.size() || feature.ray_basis.size()!=feature.rays.size()) return false;
  warp.world.resize(sample_count);warp.covariance.resize(sample_count);
  warp.bearings.resize(sample_count);warp.reference_radii.resize(sample_count);
  const M3D Rrw=feature.T_f_w_.rotationMatrix();
  const V3D trw=feature.T_f_w_.translation();
  const V3D nr=Rrw*pt.normal_,ref_center=Rrw*pt.plane_center+trw;
  const double height=nr.dot(ref_center);
  const M3D Rcr=Rcw*Rrw.transpose();
  const V3D nc=Rcw*pt.normal_;
  for(int k=0;k<sample_count;++k) {
    const V3D& br=feature.rays[k];V3D pr;
    if(!intersectPlane(br,nr,height,pr)) return false;
    const V3D world=Rrw.transpose()*(pr-trw),pc=Rcw*world+Pcw;
    const double range=pc.norm();
    if(range<1e-6 || (world-pt.plane_center).norm()>3*pt.plane_radius) return false;
    const V3D bc=pc/range;
    const double den=nr.dot(br);
    if(std::abs(den)<1e-6 || den*nc.dot(bc)<=0) return false;
    V2D uv;double cp;
    if(!sphere_image->locate(bc,uv,cp)) return false;
    const double rp=feature.ray_pitch[k];
    const M3D dp=pr.norm()*(M3D::Identity()-br*nr.transpose()/den);
    const auto basis=tangentBasis(bc);
    const Eigen::Matrix2d A=basis.transpose()*(M3D::Identity()-bc*bc.transpose())/
      range*Rcr*dp*feature.ray_basis[k];
    const Eigen::Matrix2d metric=A*A.transpose();
    // Closed-form 2x2 singular values: identical bandwidth rule, no SVD.
    const double maximum=0.5*(metric.trace()+std::hypot(metric(0,0)-metric(1,1),2*metric(0,1)));
    const double determinant=A.determinant();
    if(!(maximum>0)) return false;
    const double minimum=std::abs(determinant)/std::sqrt(maximum);
    if(!(minimum>1e-6)) return false;
    const double radius=(1<<level)*std::max({2.5*rp,2.5*cp/minimum,patch_radius/std::sqrt(double(sample_count))});
    // Match the local second moment of the old 2-D Wendland footprint:
    // covariance per tangent axis is (5/72)*radius^2. The sampler then uses
    // the conservative major sensor-pixel axis to select prefiltered levels.
    const double bound=radius*std::sqrt(maximum);
    if(!(radius>0 && radius<=0.3 && bound>0 && bound<=0.3)) return false;
    warp.world[k]=world;
    warp.covariance[k]=(5.0/72.0)*radius*radius*basis*metric*basis.transpose();
    warp.bearings[k]=bc;warp.reference_radii[k]=radius;
  }
  warp.pose_R=Rcw;warp.pose_t=Pcw;
  return configureSphericalSampling(feature,warp);
}

bool VIOManager::configureSphericalSampling(const Feature& feature,SphericalWarp& warp) const
{
  using namespace spherical;
  warp.geometry_ready=false;warp.ready=false;warp.cover=0;
  const auto& pt=*feature.point_;
  const M3D Rrw=feature.T_f_w_.rotationMatrix();
  const V3D ref_center=feature.T_f_w_*pt.plane_center,nr=Rrw*pt.normal_;
  const V3D cur_center=Rcw*pt.plane_center+Pcw,nc=Rcw*pt.normal_;
  const PlaneSupport ref_support{nr,ref_center,nr.dot(ref_center),3*pt.plane_radius};
  const PlaneSupport cur_support{nc,cur_center,nc.dot(cur_center),3*pt.plane_radius};
  const V3D center=(Rcw*pt.pos_+Pcw).normalized();
  warp.reference_filters.resize(sample_count);warp.current_filters.resize(sample_count);
  for(int k=0;k<sample_count;++k) {
    const double radius=warp.reference_radii[k];
    const auto& basis=feature.ray_basis[k];
    const Mat covariance=(5.0/72.0)*radius*radius*basis*basis.transpose();
    double ref_cover,cur_cover;
    if(!feature.sphere_image->prepare(feature.rays[k],covariance,warp.reference_filters[k],ref_cover,&ref_support) ||
       !sphere_image->prepare(warp.bearings[k],warp.covariance[k],warp.current_filters[k],cur_cover,&cur_support)) return false;
    warp.cover=std::max(warp.cover,(warp.bearings[k]-center).norm()+cur_cover);
  }
  warp.geometry_ready=true;
  return true;
}

bool VIOManager::sampleSphericalWarp(const Feature& feature,SphericalWarp& warp) const
{
  warp.ready=false;
  if(!warp.geometry_ready) return false;
  warp.reference.resize(sample_count);warp.current.resize(sample_count);
  spherical::Sample reference;
  for(int k=0;k<sample_count;++k) {
    // Geometry already validated complete prefilter supports; intensity access
    // now needs only 3 vertices per level (6 for a fractional scale).
    if(!feature.sphere_image->sample(feature.rays[k],warp.reference_filters[k],reference,false) ||
       !sphere_image->sample(warp.bearings[k],warp.current_filters[k],warp.current[k])) return false;
    warp.reference[k]=float(255*reference.value);
  }
  warp.ready=true;
  return true;
}

bool VIOManager::prepareSphericalWarp(const Feature& feature,int level,SphericalWarp& warp) const
{
  return prepareSphericalGeometry(feature,level,warp) &&
         depthConsistent(*feature.point_,warp) && sampleSphericalWarp(feature,warp);
}

bool VIOManager::prepareSphericalPyramid(const Feature& feature,std::vector<SphericalWarp>& warps) const
{
  warps.resize(patch_pyrimid_level);
  const int coarse=patch_pyrimid_level-1;
  auto& base=warps[coarse];
  if(!prepareSphericalGeometry(feature,coarse,base)) return false;
  double cover=base.cover;
  for(int level=coarse-1;level>=0;--level) {
    scaleWarpGeometry(base,level-coarse,warps[level]);
    if(!configureSphericalSampling(feature,warps[level])) return false;
    cover=std::max(cover,warps[level].cover);
  }
  // Actual triangle/prefilter supports can change with pyramid level. Inspect
  // their union envelope, rather than assuming the old kernel nesting rule.
  base.cover=cover;
  if(!depthConsistent(*feature.point_,base)) return false;
  for(int level=coarse;level>=0;--level)
    if(!sampleSphericalWarp(feature,warps[level])) return false;
  return true;
}


void VIOManager::updateDepthIndex()
{
  if(!current_points) return;
  if(depth_index_valid && (depth_pose_R-Rcw).squaredNorm()==0 && (depth_pose_t-Pcw).squaredNorm()==0) return;
  depth_heads.fill(-1);depth_rays.clear();depth_rays.reserve(current_points->size());
  depth_origin=-Rcw.transpose()*Pcw;
  for(const auto& point:*current_points) {
    const V3D delta=point.point_w-depth_origin;const double range=delta.norm();
    if(!std::isfinite(range) || range<1e-6) continue;
    const V3D ray=delta/range,bearing=Rcw*ray;
    const int key=depthKey(depthCell(bearing.x()),depthCell(bearing.y()),depthCell(bearing.z()));
    depth_rays.push_back({ray,bearing,range,std::max(0.0,ray.dot(point.var*ray)),depth_heads[key]});
    depth_heads[key]=int(depth_rays.size())-1;
  }
  depth_pose_R=Rcw;depth_pose_t=Pcw;depth_index_valid=true;
}

bool VIOManager::depthConsistent(const VisualPoint& pt,const SphericalWarp& warp) const
{
  if(!depth_index_valid || warp.world.empty() || !warp.geometry_ready) return false;
  const V3D center=(Rcw*pt.pos_+Pcw).normalized();
  const double cover=warp.cover+1e-10,cover_squared=cover*cover;
  const int xmin=depthCell(center.x()-cover),xmax=depthCell(center.x()+cover);
  const int ymin=depthCell(center.y()-cover),ymax=depthCell(center.y()+cover);
  const int zmin=depthCell(center.z()-cover),zmax=depthCell(center.z()+cover);
  for(int z=zmin;z<=zmax;++z) for(int y=ymin;y<=ymax;++y) for(int x=xmin;x<=xmax;++x)
    for(int i=depth_heads[depthKey(x,y,z)];i>=0;i=depth_rays[i].next) {
      const auto& point=depth_rays[i];
      if((point.bearing-center).squaredNorm()>cover_squared) continue;
      const double den=pt.normal_.dot(point.ray);
      if(std::abs(den)<1e-6) return false;
      const double predicted=pt.normal_.dot(pt.plane_center-depth_origin)/den;
      if(predicted<=0) return false;
      Eigen::Matrix<double,1,6> J;
      J.head<3>()=(pt.plane_center-depth_origin-predicted*point.ray).transpose()/den;
      J.tail<3>()=pt.normal_.transpose()/den;
      const double variance=std::max(0.0,(J*pt.plane_cov*J.transpose())(0,0))+
        point.variance+range_noise*range_noise+pt.roughness/(den*den);
      const double error=point.range-predicted;
      if(error*error>9*std::max(variance,1e-12)) return false;
    }
  return true;
}
