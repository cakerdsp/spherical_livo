/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#include "vio.h"
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <omp.h>

VIOManager::VIOManager()
{
  // downSizeFilter.setLeafSize(0.2, 0.2, 0.2);
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

  if (grid_size > 10)
  {
    grid_n_width = ceil(static_cast<double>(width) / grid_size);
    grid_n_height = ceil(static_cast<double>(height) / grid_size);
  }
  else
  {
    grid_size = static_cast<int>(height / grid_n_height);
    grid_n_height = ceil(static_cast<double>(height) / grid_size);
    grid_n_width = ceil(static_cast<double>(width) / grid_size);
  }
  length = grid_n_width * grid_n_height;

  grid_num.resize(length);
  map_index.resize(length);
  map_dist.resize(length);
  update_flag.resize(length);
  scan_value.resize(length);

  patch_size_total = sample_count;
  sphere_template=spherical::capTemplate(sample_count,patch_radius);
  patch_size_half = static_cast<int>(patch_size / 2);
  patch_buffer.resize(patch_size_total);
  warp_len = patch_size_total * patch_pyrimid_level;
  border = 5; // Shi-Tomasi border; full spherical footprint checked per sample.

  retrieve_voxel_points.reserve(length);
  append_voxel_points.reserve(length);

  sub_feat_map.clear();
}

void VIOManager::resetGrid()
{
  fill(grid_num.begin(), grid_num.end(), TYPE_UNKNOWN);
  fill(map_index.begin(), map_index.end(), 0);
  fill(map_dist.begin(), map_dist.end(), 10000.0f);
  fill(update_flag.begin(), update_flag.end(), 0);
  fill(scan_value.begin(), scan_value.end(), 0.0f);

  retrieve_voxel_points.clear();
  retrieve_voxel_points.resize(length);

  append_voxel_points.clear();
  append_voxel_points.resize(length);

  total_points = 0;
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
  if (feat_map.size() <= 0) return;
  double ts0 = omp_get_wtime();

  // pg_down->reserve(feat_map.size());
  // downSizeFilter.setInputCloud(pg);
  // downSizeFilter.filter(*pg_down);

  // resetRvizDisplay();

  // Controls whether to include the visual submap from the previous frame.
  sub_feat_map.clear();

  float voxel_size = 0.5;

  double loc_xyz[3];

  // printf("A0. initial depthmap: %.6lf \n", omp_get_wtime() - ts0);
  // double ts1 = omp_get_wtime();

  // printf("pg size: %zu \n", pg.size());

  for (int i = 0; i < pg.size(); i++)
  {
    // double t0 = omp_get_wtime();

    V3D pt_w = pg[i].point_w;

    for (int j = 0; j < 3; j++)
    {
      loc_xyz[j] = pt_w[j] / voxel_size;
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
    }
    VOXEL_LOCATION position(loc_xyz[0], loc_xyz[1], loc_xyz[2]);

    // t_position += omp_get_wtime()-t0;
    // double t1 = omp_get_wtime();

    auto iter = sub_feat_map.find(position);
    if (iter == sub_feat_map.end()) { sub_feat_map[position] = 0; }
    else { iter->second = 0; }

    // t_insert += omp_get_wtime()-t1;
    // double t2 = omp_get_wtime();

    // t_depth += omp_get_wtime()-t2;
  }

  // imshow("depth_img", depth_img);
  // printf("A1: %.6lf \n", omp_get_wtime() - ts1);
  // printf("A11. calculate pt position: %.6lf \n", t_position);
  // printf("A12. sub_postion.insert(position): %.6lf \n", t_insert);
  // printf("A13. generate depth map: %.6lf \n", t_depth);
  // printf("A. projection: %.6lf \n", omp_get_wtime() - ts0);

  // double t1 = omp_get_wtime();
  vector<VOXEL_LOCATION> DeleteKeyList;

  for (auto &iter : sub_feat_map)
  {
    VOXEL_LOCATION position = iter.first;

    // double t4 = omp_get_wtime();
    auto corre_voxel = feat_map.find(position);
    // double t5 = omp_get_wtime();

    if (corre_voxel != feat_map.end())
    {
      bool voxel_in_fov = false;
      std::vector<VisualPoint *> &voxel_points = corre_voxel->second->voxel_points;
      int voxel_num = voxel_points.size();

      for (int i = 0; i < voxel_num; i++)
      {
        VisualPoint *pt = voxel_points[i];
        if (pt == nullptr) continue;
        if (pt->obs_.size() == 0) continue;

        V3D norm_vec(new_frame_->T_f_w_.rotationMatrix() * pt->normal_);
        V3D dir(new_frame_->T_f_w_ * pt->pos_);

        // dir.normalize();
        // if (dir.dot(norm_vec) <= 0.17) continue; // 0.34 70 degree  0.17 80 degree 0.08 85 degree

        V2D pc(new_frame_->w2c(pt->pos_));
        if (new_frame_->cam_->isInFrame(pc.cast<int>(), border))
        {
          // cv::circle(img_cp, cv::Point2f(pc[0], pc[1]), 3, cv::Scalar(0, 255, 255), -1, 8);
          voxel_in_fov = true;
          int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size);
          grid_num[index] = TYPE_MAP;
          Vector3d obs_vec(new_frame_->pos() - pt->pos_);
          float cur_dist = obs_vec.norm();
          if (cur_dist <= map_dist[index])
          {
            map_dist[index] = cur_dist;
            retrieve_voxel_points[index] = pt;
          }
        }
      }
      if (!voxel_in_fov) { DeleteKeyList.push_back(position); }
    }
  }

  for (auto &key : DeleteKeyList)
  {
    sub_feat_map.erase(key);
  }

  // double t2 = omp_get_wtime();

  // cout<<"B. feat_map.find: "<<t2-t1<<endl;

  // double t_2, t_3, t_4, t_5;
  // t_2=t_3=t_4=t_5=0;

  for (int i = 0; i < length; i++)
  {
    if (grid_num[i] == TYPE_MAP)
    {
      // double t_1 = omp_get_wtime();

      VisualPoint *pt = retrieve_voxel_points[i];
      // visual_sub_map_cur.push_back(pt); // before

      V2D pc(new_frame_->w2c(pt->pos_));

      // cv::circle(img_cp, cv::Point2f(pc[0], pc[1]), 3, cv::Scalar(0, 0, 255), -1, 8); // Green Sparse Align tracked

      if (total_points >= max_patches || !refreshPlane(*pt)) continue;
      Feature* ref_ftr=nullptr;
      if(!pt->getCloseViewObs(new_frame_->pos(),ref_ftr,pc) || !ref_ftr->sphere_image) continue;
      std::vector<SphericalWarp> warps(patch_pyrimid_level);
      bool valid=true;
      for(int level=0;level<patch_pyrimid_level;++level)
        if(!prepareSphericalWarp(*ref_ftr,level,warps[level]) || !depthConsistent(*pt,warps[level])) { valid=false; break; }
      if(!valid) continue;
      float error=0;
      for(int k=0;k<patch_size_total;++k) {
        spherical::Sample sample;
        const V3D bc=new_frame_->w2f(warps[0].world[k]).normalized();
        if(!sphere_image->sample(bc,warps[0].precision[k],sample)) { valid=false; break; }
        const double residual=255*state->inv_expo_time*sample.value-ref_ftr->inv_expo_time_*warps[0].reference[k];
        error+=residual*residual;
      }
      if(!valid || error > outlier_threshold*patch_size_total) continue;
      visual_submap->voxel_points.push_back(pt);
      visual_submap->references.push_back(ref_ftr);
      visual_submap->sphere_warps.push_back(std::move(warps));
      visual_submap->propa_errors.push_back(error);
      visual_submap->errors.push_back(error);
      visual_submap->search_levels.push_back(0);
      visual_submap->inv_expo_list.push_back(ref_ftr->inv_expo_time_);
      ++total_points;
      // t_5 += omp_get_wtime() - t_1;
    }
  }
  total_points = visual_submap->voxel_points.size();

  // double t3 = omp_get_wtime();
  // cout<<"C. addSubSparseMap: "<<t3-t2<<endl;
  // cout<<"depthcontinuous: C1 "<<t_2<<" C2 "<<t_3<<" C3 "<<t_4<<" C4
  // "<<t_5<<endl;
  printf("[ VIO ] Retrieve %d points from visual sparse map\n", total_points);
}

void VIOManager::computeJacobianAndUpdateEKF(cv::Mat img)
{
  compute_jacobian_time = update_ekf_time = 0;
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
  if (pg.size() <= 10) return;

  // double t0 = omp_get_wtime();
  for (int i = 0; i < pg.size(); i++)
  {
    if (pg[i].normal == V3D(0, 0, 0)) continue;

    V3D pt = pg[i].point_w;
    V2D pc(new_frame_->w2c(pt));

    if (new_frame_->cam_->isInFrame(pc.cast<int>(), border)) // 20px is the patch size in the matcher
    {
      int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size);

      if (grid_num[index] != TYPE_MAP)
      {
        float cur_value = vk::shiTomasiScore(img, pc[0], pc[1]);
        // if (cur_value < 5) continue;
        if (cur_value > scan_value[index])
        {
          scan_value[index] = cur_value;
          append_voxel_points[index] = pg[i];
          grid_num[index] = TYPE_POINTCLOUD;
        }
      }
    }
  }

  for (int j = 0; j < visual_submap->add_from_voxel_map.size(); j++)
  {
    V3D pt = visual_submap->add_from_voxel_map[j].point_w;
    V2D pc(new_frame_->w2c(pt));

    if (new_frame_->cam_->isInFrame(pc.cast<int>(), border)) // 20px is the patch size in the matcher
    {
      int index = static_cast<int>(pc[1] / grid_size) * grid_n_width + static_cast<int>(pc[0] / grid_size);

      if (grid_num[index] != TYPE_MAP)
      {
        float cur_value = vk::shiTomasiScore(img, pc[0], pc[1]);
        if (cur_value > scan_value[index])
        {
          scan_value[index] = cur_value;
          append_voxel_points[index] = visual_submap->add_from_voxel_map[j];
          grid_num[index] = TYPE_POINTCLOUD;
        }
      }
    }
  }

  // double t_b1 = omp_get_wtime() - t0;
  // t0 = omp_get_wtime();

  int add = 0;
  for (int i = 0; i < length; i++)
  {
    if (add + total_points >= max_patches) break;
    if (grid_num[i] == TYPE_POINTCLOUD) // && (scan_value[i]>=50))
    {
      pointWithVar pt_var = append_voxel_points[i];
      V3D pt = pt_var.point_w;

      V3D norm_vec(new_frame_->T_f_w_.rotationMatrix() * pt_var.normal);
      V3D dir(new_frame_->T_f_w_ * pt);
      dir.normalize();
      double cos_theta = dir.dot(norm_vec);
      // if(std::fabs(cos_theta)<0.34) continue; // 70 degree
      V2D pc(new_frame_->w2c(pt));

      float *patch = new float[patch_size_total];

      VisualPoint *pt_new = new VisualPoint(pt);

      Vector3d f = cam->cam2world(pc);
      Feature *ftr_new = new Feature(pt_new, patch, pc, f, new_frame_->T_f_w_, 0);
      ftr_new->img_ = img;
      ftr_new->id_ = new_frame_->id_;
      ftr_new->inv_expo_time_ = state->inv_expo_time;

      pt_new->addFrameRef(ftr_new);
      pt_new->covariance_ = pt_var.var;
      pt_new->is_normal_initialized_ = true;

      if (cos_theta < 0) { pt_new->normal_ = -pt_var.normal; }
      else { pt_new->normal_ = pt_var.normal; }
      
      pt_new->previous_normal_ = pt_new->normal_;

      if(!refreshPlane(*pt_new) || !buildSphericalReference(*ftr_new)) { delete pt_new; continue; }
      insertPointIntoVoxelMap(pt_new);
      add += 1;
      // map_cur_frame.push_back(pt_new);
    }
  }

  // double t_b2 = omp_get_wtime() - t0;

  printf("[ VIO ] Append %d new visual map points\n", add);
  // printf("pg.size: %d \n", pg.size());
  // printf("B1. : %.6lf \n", t_b1);
  // printf("B2. : %.6lf \n", t_b2);
}

void VIOManager::updateVisualMapPoints(cv::Mat img)
{
  if (total_points == 0) return;

  int update_num = 0;
  SE3 pose_cur = new_frame_->T_f_w_;
  for (int i = 0; i < total_points; i++)
  {
    VisualPoint *pt = visual_submap->voxel_points[i];
    if (pt == nullptr) continue;
    if (pt->is_converged_)
    { 
      pt->deleteNonRefPatchFeatures();
      continue;
    }

    V2D pc(new_frame_->w2c(pt->pos_));
    bool add_flag = false;
    
    float *patch_temp = new float[patch_size_total];

    // TODO: condition: distance and view_angle
    // Step 1: time
    Feature *last_feature = pt->obs_.front();
    // if(new_frame_->id_ >= last_feature->id_ + 10) add_flag = true; // 10

    // Step 2: delta_pose
    SE3 pose_ref = last_feature->T_f_w_;
    SE3 delta_pose = pose_ref * pose_cur.inverse();
    double delta_p = delta_pose.translation().norm();
    double delta_theta = (delta_pose.rotationMatrix().trace() > 3.0 - 1e-6) ? 0.0 : std::acos(0.5 * (delta_pose.rotationMatrix().trace() - 1));
    if (delta_p > 0.5 || delta_theta > 0.3) add_flag = true; // 0.5 || 0.3

    // Step 3: pixel distance
    Vector2d last_px = last_feature->px_;
    double pixel_dist = (pc - last_px).norm();
    if (last_feature->camera_id != active_camera || pixel_dist > 40) add_flag = true;

    // Maintain the size of 3D point observation features.
    if (pt->obs_.size() >= 30)
    {
      Feature *ref_ftr;
      pt->findMinScoreFeature(new_frame_->pos(), ref_ftr);
      pt->deleteFeatureRef(ref_ftr);
      // cout<<"pt->obs_.size() exceed 20 !!!!!!"<<endl;
    }
    if (add_flag)
    {
      update_num += 1;
      update_flag[i] = 1;
      Vector3d f = cam->cam2world(pc);
      Feature *ftr_new = new Feature(pt, patch_temp, pc, f, new_frame_->T_f_w_, visual_submap->search_levels[i]);
      ftr_new->img_ = img;
      ftr_new->id_ = new_frame_->id_;
      ftr_new->inv_expo_time_ = state->inv_expo_time;
      if(!buildSphericalReference(*ftr_new)) { delete ftr_new; continue; }
      pt->addFrameRef(ftr_new);
    }
    else delete[] patch_temp;
  }
  printf("[ VIO ] Update %d points in visual submap\n", update_num);
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
  updateFrameState(*state);
  std::vector<SphericalWarp> prepared(total_points);
  for(int i=0;i<total_points;++i) {
    if(!prepareSphericalWarp(*visual_submap->references[i],level,prepared[i]) ||
       !depthConsistent(*visual_submap->voxel_points[i],prepared[i])) return;
  }
  for(int i=0;i<total_points;++i) visual_submap->sphere_warps[i][level]=std::move(prepared[i]);

  VectorXd z;
  MatrixXd H_sub;
  bool EKF_end = false;
  float last_error = std::numeric_limits<float>::max();

  const int H_DIM = total_points * patch_size_total;
  z.resize(H_DIM);
  z.setZero();
  H_sub.resize(H_DIM, 7);
  H_sub.setZero();

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
    #pragma omp parallel for num_threads(4) reduction(+:error,n_meas,invalid)
    for(int i=0;i<total_points;++i)
    {
      const auto& warp=visual_submap->sphere_warps[i][level];
      const double inv_ref_expo=visual_submap->inv_expo_list[i];
      float patch_error=0;
      for(int k=0;k<patch_size_total;++k) {
        const V3D pc=Rcw*warp.world[k]+Pcw;
        spherical::Sample sample;
        if(pc.squaredNorm()<1e-12 || !sphere_image->sample(pc.normalized(),warp.precision[k],sample)) { ++invalid; break; }
        const double cur_value=255*sample.value;
        const double res=state->inv_expo_time*cur_value-inv_ref_expo*warp.reference[k];
        const int row=i*patch_size_total+k;
        z(row)=res;
        H_sub.block<1,6>(row,0)=255*state->inv_expo_time*sample.gradient.transpose()*
          spherical::bearingJacobian(Rwi,Pwi,Rci,Pci,warp.world[k]);
        if(exposure_estimate_en) H_sub(row,6)=cur_value;
        patch_error+=res*res; ++n_meas;
      }
      visual_submap->errors[i]=patch_error;
      error+=patch_error;
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
  uint8_t *img_ptr = (uint8_t *)img.data + ((v_ref_i)*width + (u_ref_i)) * 3;
  float B = w_ref_tl * img_ptr[0] + w_ref_tr * img_ptr[0 + 3] + w_ref_bl * img_ptr[width * 3] + w_ref_br * img_ptr[width * 3 + 0 + 3];
  float G = w_ref_tl * img_ptr[1] + w_ref_tr * img_ptr[1 + 3] + w_ref_bl * img_ptr[1 + width * 3] + w_ref_br * img_ptr[width * 3 + 1 + 3];
  float R = w_ref_tl * img_ptr[2] + w_ref_tr * img_ptr[2 + 3] + w_ref_bl * img_ptr[2 + width * 3] + w_ref_br * img_ptr[width * 3 + 2 + 3];
  V3F pixel(B, G, R);
  return pixel;
}

void VIOManager::processFrame(cv::Mat &img, vector<pointWithVar> &pg, const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &feat_map, double img_time)
{
  if(img.empty() || width!=img.cols || height!=img.rows)
    throw std::runtime_error("raw image dimensions differ from camera calibration");
  current_points=&pg; current_planes=&feat_map;
  img_rgb = img.clone();
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
  ave_total = ave_total * (frame_count - 1) / frame_count + (t7 - t1 - (t5 - t4)) / frame_count;

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
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;34m| %-29s | %-27s |\033[0m\n", "Algorithm Stage", "Time (secs)");
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "retrieveFromVisualSparseMap", t2 - t1);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "computeJacobianAndUpdateEKF", t3 - t2);
  printf("\033[1;32m| %-27s   | %-27lf |\033[0m\n", "-> computeJacobian", compute_jacobian_time);
  printf("\033[1;32m| %-27s   | %-27lf |\033[0m\n", "-> updateEKF", update_ekf_time);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "generateVisualMapPoints", t4 - t3);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "updateVisualMapPoints", t6 - t5);
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "updateReferencePatch", t7 - t6);
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "Current Total Time", t7 - t1 - (t5 - t4));
  printf("\033[1;32m| %-29s | %-27lf |\033[0m\n", "Average Total Time", ave_total);
  printf("\033[1;34m+-------------------------------------------------------------+\033[0m\n");

  // std::string text = std::to_string(int(1 / (t7 - t1 - (t5 - t4)))) + " HZ";
  // cv::Point2f origin;
  // origin.x = 20;
  // origin.y = 20;
  // cv::putText(img_cp, text, origin, cv::FONT_HERSHEY_COMPLEX, 0.6, cv::Scalar(255, 255, 255), 1, 8, 0);
  // cv::imwrite("/home/chunran/Desktop/raycasting/" + std::to_string(new_frame_->id_) + ".png", img_cp);
}// Native spherical geometry integrated into the upstream VIOManager.
// The shared VisualPoint/Feature map, selection, reference lifecycle and IEKF
// remain in vio.cpp; this replaces planar patch construction and affine warping.
void VIOManager::activateCamera(int id)
{
  active_camera=id;
  auto& view=cameras_.at(id);
  cam=view.camera.get(); Rcl=view.Rcl; Pcl=view.Pcl;
  if(sample_count<4 || max_patches<1 || !(patch_radius>0 && patch_radius<0.5))
    throw std::invalid_argument("invalid spherical patch support");
  if(grid_size<1) grid_size=std::max(10,cam->height()/17);
  initializeVIO();
}

bool VIOManager::getColorFromCamera(int id,const V3D& world,const StatesGroup& pose,V3F& color,double blind)
{
  const auto& view=cameras_.at(id);
  if(view.img_rgb.empty()) return false;
  const M3D rci=view.Rcl*Rli;
  const V3D pci=view.Rcl*Pli+view.Pcl;
  const V3D pc=rci*pose.rot_end.transpose()*(world-pose.pos_end)+pci;
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
  for(const auto& b:sphere_template) feature.rays.push_back(rotation*b);
  // Check the complete coarsest footprint before retaining a reference.
  SphericalWarp warp;
  if(!prepareSphericalWarp(feature,patch_pyrimid_level-1,warp) || !depthConsistent(*feature.point_,warp)) return false;
  if(!prepareSphericalWarp(feature,0,warp)) return false;
  std::copy(warp.reference.begin(),warp.reference.end(),feature.patch_);
  return true;
}

bool VIOManager::prepareSphericalWarp(const Feature& feature,int level,SphericalWarp& warp) const
{
  using namespace spherical;
  warp=SphericalWarp{};
  const auto& pt=*feature.point_;
  if(!feature.sphere_image || feature.rays.size()!=size_t(sample_count)) return false;
  const M3D Rrw=feature.T_f_w_.rotationMatrix();
  const V3D trw=feature.T_f_w_.translation();
  const V3D nr=Rrw*pt.normal_;
  const double height=nr.dot(Rrw*pt.plane_center+trw);
  const M3D Rcr=Rcw*Rrw.transpose();
  const V3D nc=Rcw*pt.normal_;
  const double current_height=nc.dot(Rcw*pt.plane_center+Pcw);
  const auto finiteSupport=[&](const Image& image,const Sample& sample,const M3D& rotation,
                             const V3D& translation,const V3D& normal,double h) {
    for(const auto& weight:sample.weights) {
      const int width=image.atlas().width();
      const V3D ray=image.atlas().ray(weight.first%width,weight.first/width).normalized();
      V3D surface;
      if(!intersectPlane(ray,normal,h,surface)) return false;
      const V3D world=rotation.transpose()*(surface-translation);
      if((world-pt.plane_center).norm()>3*pt.plane_radius) return false;
    }
    return true;
  };
  for(const V3D& br:feature.rays) {
    V3D pr;
    if(!intersectPlane(br,nr,height,pr)) return false;
    const V3D world=Rrw.transpose()*(pr-trw), pc=Rcw*world+Pcw;
    if(pc.squaredNorm()<1e-12 || (world-pt.plane_center).norm()>3*pt.plane_radius) return false;
    const V3D bc=pc.normalized();
    const double den=nr.dot(br);
    if(std::abs(den)<1e-6 || den*nc.dot(bc)<=0) return false;
    V2D uv; double rp,cp;
    if(!feature.sphere_image->locate(br,uv,rp) || !sphere_image->locate(bc,uv,cp)) return false;
    const M3D dp=pr.norm()*(M3D::Identity()-br*nr.transpose()/den);
    const Eigen::Matrix2d A=tangentBasis(bc).transpose()*(M3D::Identity()-bc*bc.transpose())/
      pc.norm()*Rcr*dp*tangentBasis(br);
    const double minimum=Eigen::JacobiSVD<Eigen::Matrix2d>(A).singularValues().minCoeff();
    if(!(minimum>1e-6)) return false;
    const double radius=(1<<level)*std::max({2.5*rp,2.5*cp/minimum,patch_radius/std::sqrt(double(sample_count))});
    const Eigen::Matrix2d footprint=Eigen::Matrix2d::Identity()*radius*radius;
    const M3D ar=precision(br,footprint), ac=precision(bc,A*footprint*A.transpose());
    Sample reference,current;
    if(!feature.sphere_image->sample(br,ar,reference) || !sphere_image->sample(bc,ac,current)) return false;
    if(!finiteSupport(*feature.sphere_image,reference,Rrw,trw,nr,height) ||
       !finiteSupport(*sphere_image,current,Rcw,Pcw,nc,current_height)) return false;
    warp.world.push_back(world); warp.precision.push_back(ac);
    warp.reference.push_back(float(255*reference.value));
  }
  return true;
}

bool VIOManager::depthConsistent(const VisualPoint& pt,const SphericalWarp& warp) const
{
  if(!current_points || warp.world.empty()) return false;
  const V3D origin=-Rcw.transpose()*Pcw;
  const V3D center=(Rcw*pt.pos_+Pcw).normalized();
  double cover=0;
  for(size_t k=0;k<warp.world.size();++k) {
    const V3D b=(Rcw*warp.world[k]+Pcw).normalized();
    const auto basis=spherical::tangentBasis(b);
    const Eigen::Matrix2d precision=basis.transpose()*warp.precision[k]*basis;
    const double minimum=Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>(precision).eigenvalues().minCoeff();
    if(!(minimum>0)) return false;
    cover=std::max(cover,(b-center).norm()+1/std::sqrt(minimum));
  }
  for(const auto& point:*current_points) {
    const V3D delta=point.point_w-origin; const double measured=delta.norm();
    if(measured<1e-6) continue;
    const V3D ray=delta/measured;
    if((Rcw*ray-center).norm()>cover) continue;
    const double den=pt.normal_.dot(ray);
    if(std::abs(den)<1e-6) return false;
    const double predicted=pt.normal_.dot(pt.plane_center-origin)/den;
    if(predicted<=0) return false;
    Eigen::Matrix<double,1,6> J;
    J.head<3>()=(pt.plane_center-origin-predicted*ray).transpose()/den;
    J.tail<3>()=pt.normal_.transpose()/den;
    const double variance=std::max(0.0,(J*pt.plane_cov*J.transpose())(0,0))+
      std::max(0.0,ray.dot(point.var*ray))+range_noise*range_noise+pt.roughness/(den*den);
    if(std::abs(measured-predicted)>3*std::sqrt(std::max(variance,1e-12))) return false;
  }
  return true;
}
