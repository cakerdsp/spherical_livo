#include "spherical_vio.h"
#include "timing.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace spherical {
namespace {
constexpr double huber = 1.345;
double loss(double r) {const double a=std::abs(r);return a<=huber?0.5*r*r:huber*(a-0.5*huber);}
double weight(double r) {return std::abs(r)<=huber?1.0:huber/std::abs(r);}
Vec cameraCenter(const Camera& c,const StatesGroup& s) {return s.pos_end-s.rot_end*c.Rci.transpose()*c.tci;}
Vec toCamera(const Camera& c,const StatesGroup& s,const Vec& p) {return c.Rci*s.rot_end.transpose()*(p-s.pos_end)+c.tci;}
Eigen::Matrix<double,1,6> rangePlaneJacobian(const Vec& origin,const Vec& ray,const Vec& normal,
                                           const Vec& center,double distance) {
  Eigen::Matrix<double,1,6> j;
  const double den=normal.dot(ray);
  j.head<3>()=(center-origin-distance*ray).transpose()/den;
  j.tail<3>()=normal.transpose()/den;return j;
}
}

VisualEstimator::VisualEstimator(VisualConfig config,std::vector<Camera> cameras)
 :config_(config),cameras_(std::move(cameras)),template_(capTemplate(config.samples,config.radius)) {
  if(cameras_.empty()||config_.max_patches<4) throw std::invalid_argument("invalid visual budget/camera list");
}

bool VisualEstimator::consistent(const Landmark& lm,const Camera& camera,const StatesGroup& state,
                                 const std::vector<pointWithVar>& points,double range_noise,double full_cover) const {
  const Vec origin=cameraCenter(camera,state),center=toCamera(camera,state,lm.center).normalized();
  const Mat Rcw=camera.Rci*state.rot_end.transpose();
  double cover=0;
  for(const Vec& p:lm.world) cover=std::max(cover,(toCamera(camera,state,p).normalized()-center).norm());
  // Cover the largest (4x) smoothing footprint too, rather than just sample centres.
  Eigen::Vector2d uv;double pitch;
  if(!camera.atlas->locate(center,uv,pitch)) return false;
  cover=std::max(cover+16*pitch,full_cover);
  if(lm.normal.dot(origin-lm.center)<=0) return false;
  for(const auto& point:points) {
    const Vec delta=point.point_w-origin;const double measured=delta.norm();
    if(measured<=1e-6) continue;
    const Vec ray=delta/measured;
    if((Rcw*ray-center).norm()>cover) continue;
    const double den=lm.normal.dot(ray);
    if(std::abs(den)<1e-6) return false;
    const double predicted=lm.normal.dot(lm.plane_center-origin)/den;
    if(predicted<=0) return false;
    const auto j=rangePlaneJacobian(origin,ray,lm.normal,lm.plane_center,predicted);
    const double variance=std::max(0.0,(j*lm.plane_cov*j.transpose())(0,0))+
      std::max(0.0,ray.dot(point.var*ray))+range_noise*range_noise+lm.roughness/(den*den);
    if(std::abs(measured-predicted)>3*std::sqrt(std::max(variance,1e-12))) return false;
  }
  return true;
}

bool VisualEstimator::prepare(const Landmark& lm,const Image& image,const Camera& camera,
                              const StatesGroup& state,int scale,Track& track) const {
  track.point=&lm;track.precision.clear();track.reference.clear();
  const Mat Rcr=camera.Rci*state.rot_end.transpose()*lm.Rrw.transpose();
  const Vec nr=lm.Rrw*lm.normal;
  for(size_t k=0;k<lm.world.size();++k) {
    const Vec pr=lm.Rrw*lm.world[k]+lm.trw,pc=toCamera(camera,state,lm.world[k]);
    if(pc.norm()<1e-6) return false;
    const Vec br=lm.rays[k],bc=pc.normalized();
    Eigen::Vector2d uv;double rp,cp;
    if(!lm.image->locate(br,uv,rp)||!image.locate(bc,uv,cp)) return false;
    const double den=nr.dot(br);if(std::abs(den)<1e-6) return false;
    const Mat dp=pr.norm()*(Mat::Identity()-br*nr.transpose()/den);
    const Eigen::Matrix2d warp=tangentBasis(bc).transpose()*(Mat::Identity()-bc*bc.transpose())/
      pc.norm()*Rcr*dp*tangentBasis(br);
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(warp);
    const double minimum=svd.singularValues().minCoeff();
    if(!(minimum>1e-6)) return false;
    // Native pixel footprints from BOTH cameras bound the common physical
    // bandwidth. Warp transfers the whole ellipse; no scalar scale is fitted.
    const double radius=scale*std::max({2.5*rp,2.5*cp/minimum,
      config_.radius/std::sqrt(double(config_.samples))});
    const Eigen::Matrix2d ref=Eigen::Matrix2d::Identity()*radius*radius;
    const Eigen::Matrix2d cur=warp*ref*warp.transpose();
    const Mat ar=precision(br,ref),ac=precision(bc,cur);
    Sample sample;
    if(!lm.image->sample(br,ar,sample)) return false;
    // All contributing reference rays, including the low-pass footprint, must
    // remain on the finite measured surface support.
    const Vec origin=-lm.Rrw.transpose()*lm.trw;
    const double h=nr.dot(lm.Rrw*lm.plane_center+lm.trw);
    for(const auto& weight:sample.weights){
      const int width=lm.image->atlas().width();
      const Vec ray=lm.image->atlas().ray(weight.first%width,weight.first/width).normalized();
      Vec surface;
      if(!intersectPlane(ray,nr,h,surface))return false;
      const Vec world=origin+lm.Rrw.transpose()*surface;
      if((world-lm.plane_center).norm()>3*lm.plane_radius)return false;
    }
    track.reference.push_back(std::move(sample));track.precision.push_back(ac);
    if(!image.sample(bc,ac,sample)) return false;
  }
  return true;
}

bool VisualEstimator::linearize(const std::vector<Track>& tracks,const Image& image,const Camera& camera,
                                const StatesGroup& state,double log_gain,const StatesGroup& prior,
                                const Eigen::Matrix<double,18,18>& information,Linearization& out,
                                const std::vector<Eigen::MatrixXd>* fixed) const {
  out=Linearization{};
  if(!std::isfinite(log_gain)||std::abs(log_gain)>10) return false;
  const double gain=std::exp(log_gain);
  const auto difference=state-prior;
  Eigen::Matrix<double,18,18> jp=Eigen::Matrix<double,18,18>::Identity();
  jp.topLeftCorner<3,3>()=rightJacobianInverse(difference.head<3>());
  out.H.topLeftCorner<18,18>()=jp.transpose()*information*jp;
  out.g.head<18>()=jp.transpose()*information*difference;
  out.cost=0.5*difference.dot(information*difference);
  const Mat Rcw=camera.Rci*state.rot_end.transpose();
  for(size_t ti=0;ti<tracks.size();++ti) {
    const auto& tr=tracks[ti];const auto& lm=*tr.point;const int n=lm.world.size();
    std::vector<Sample> cur(n);
    Eigen::VectorXd residual(n);
    Eigen::MatrixXd j=Eigen::MatrixXd::Zero(n,19),jg(n,6),cov(n,n);
    for(int k=0;k<n;++k) {
      const Vec pc=toCamera(camera,state,lm.world[k]),bc=pc.normalized();
      if(!image.sample(bc,tr.precision[k],cur[k])) return false;
      residual[k]=gain*cur[k].value-lm.gain*tr.reference[k].value;
      j.block<1,6>(k,0)=gain*cur[k].gradient.transpose()*
        bearingJacobian(state.rot_end,state.pos_end,camera.Rci,camera.tci,lm.world[k]);
      j(k,18)=gain*cur[k].value;
      jg.row(k)=gain*cur[k].gradient.transpose()*(Mat::Identity()-bc*bc.transpose())/
        pc.norm()*Rcw*lm.plane_jacobians[k];
    }
    Eigen::MatrixXd whitening;
    if(fixed) whitening=(*fixed)[ti];
    else {
      const double cs=std::pow(gain*image.noise(),2),rs=std::pow(lm.gain*lm.image->noise(),2);
      for(int k=0;k<n;++k)for(int l=0;l<=k;++l) {
        const double v=cs*weightOverlap(cur[k],cur[l])+rs*weightOverlap(tr.reference[k],tr.reference[l]);
        cov(k,l)=cov(l,k)=v;
      }
      cov.noalias()+=jg*lm.plane_cov*jg.transpose();
      cov.diagonal().array()+=std::max(1e-12,cov.trace()/n*1e-8);
      Eigen::LLT<Eigen::MatrixXd> llt(cov);if(llt.info()!=Eigen::Success) return false;
      whitening=llt.matrixL().solve(Eigen::MatrixXd::Identity(n,n));
    }
    const Eigen::VectorXd r=whitening*residual;const Eigen::MatrixXd h=whitening*j;
    for(int k=0;k<n;++k) {
      const double w=weight(r[k]);out.cost+=loss(r[k]);
      out.H.noalias()+=w*h.row(k).transpose()*h.row(k);
      out.g.noalias()+=w*h.row(k).transpose()*r[k];
    }
    out.squared_error+=residual.squaredNorm();out.samples+=n;
    out.whitening.push_back(std::move(whitening));
  }
  return out.samples>0&&out.H.allFinite()&&out.g.allFinite()&&std::isfinite(out.cost);
}

VisualStats VisualEstimator::process(int id,const cv::Mat& gray,StatesGroup& state,
                                     const std::vector<pointWithVar>& points,const VoxelMapManager& geometry,
                                     std::vector<VisualPatch>* display) {
  const auto start=SteadyClock::now();
  if(display)display->clear();
  std::vector<const Landmark*> final_landmarks;
  if(id<0||id>=int(cameras_.size()))throw std::out_of_range("camera id");
  Camera& camera=cameras_[id];auto image=std::make_shared<Image>(camera.atlas,gray);
  VisualStats stats;std::vector<const Landmark*> candidates;
  for(const auto& p:map_) {
    Eigen::Vector2d uv;double pitch;
    if(!image->locate(toCamera(camera,state,p->center).normalized(),uv,pitch)) continue;
    ++stats.candidates;
    if(!consistent(*p,camera,state,points,geometry.config_setting_.dept_err_)){++stats.depth_rejected;continue;}
    candidates.push_back(p.get());
  }
  std::stable_sort(candidates.begin(),candidates.end(),[&](const Landmark* a,const Landmark* b){
    return (a->center-state.pos_end).squaredNorm()<(b->center-state.pos_end).squaredNorm();});
  stats.candidates_ms=elapsedMs(start);
  const auto optimization_start=SteadyClock::now();
  const StatesGroup prior=state;
  Eigen::LLT<Eigen::Matrix<double,18,18>> pll(0.5*(prior.cov+prior.cov.transpose()));
  if(pll.info()!=Eigen::Success)throw std::runtime_error("non-positive inertial covariance");
  const Eigen::Matrix<double,18,18> information=pll.solve(Eigen::Matrix<double,18,18>::Identity());
  double log_gain=camera.last_log_gain;Linearization final_system;bool have_final=false;
  StatesGroup last_valid=prior;double last_valid_gain=log_gain;
  for(int scale: {4,2,1}) {
    const auto prepare_start=SteadyClock::now();
    std::vector<Track> tracks;std::vector<std::pair<Vec,double>> footprints;
    for(const auto* p:candidates) {
      Track tr;if(!prepare(*p,*image,camera,state,scale,tr))continue;
      const Vec centre=toCamera(camera,state,p->center).normalized();double radius=0;
      for(size_t k=0;k<p->world.size();++k){
        Eigen::SelfAdjointEigenSolver<Mat> es(tr.precision[k]);
        radius=std::max(radius,(toCamera(camera,state,p->world[k]).normalized()-centre).norm()+
                        1/std::sqrt(es.eigenvalues().minCoeff()));
      }
      if(!consistent(*p,camera,state,points,geometry.config_setting_.dept_err_,radius))continue;
      bool overlap=false;for(const auto& f:footprints)if((centre-f.first).norm()<radius+f.second){overlap=true;break;}
      if(overlap)continue;footprints.emplace_back(centre,radius);tracks.push_back(std::move(tr));
      if(tracks.size()>=size_t(config_.max_patches))break;
    }
    stats.prepare_ms+=elapsedMs(prepare_start);
    if(tracks.size()<3)continue;
    for(int iteration=0;iteration<8;++iteration) {
      Linearization lin;if(!linearize(tracks,*image,camera,state,log_gain,prior,information,lin))break;
      Eigen::LDLT<Eigen::Matrix<double,19,19>> solver(lin.H);
      if(solver.info()!=Eigen::Success||!solver.isPositive())break;
      const Eigen::Matrix<double,19,1> delta=solver.solve(-lin.g);if(!delta.allFinite())break;
      bool accepted=false;
      for(double step=1;step>=1.0/64;step*=0.5) {
        StatesGroup trial=state;trial+=Eigen::Matrix<double,18,1>(step*delta.head<18>());
        Linearization evaluated;
        if(linearize(tracks,*image,camera,trial,log_gain+step*delta[18],prior,information,evaluated,&lin.whitening)&&
           evaluated.cost<=lin.cost+1e-10) {
          state=trial;log_gain+=step*delta[18];accepted=true;break;
        }
      }
      if(!accepted||delta.head<6>().norm()<1e-6)break;
    }
    Linearization lin;
    if(linearize(tracks,*image,camera,state,log_gain,prior,information,lin)){
      final_system=std::move(lin);have_final=true;stats.patches=tracks.size();stats.cross_camera=0;
      last_valid=state;last_valid_gain=log_gain;
      final_landmarks.clear();
      for(const auto& t:tracks) {
        if(t.point->camera!=id)++stats.cross_camera;
        if(display)final_landmarks.push_back(t.point);
      }
    }else{state=last_valid;log_gain=last_valid_gain;}
  }
  state=last_valid;log_gain=last_valid_gain;
  if(have_final) {
    Eigen::LDLT<Eigen::Matrix<double,19,19>> solver(final_system.H);
    if(solver.info()==Eigen::Success&&solver.isPositive()) {
      // Marginalizing the frame's log gain retains its coupling with pose.
      const Eigen::Matrix<double,19,19> posterior=solver.solve(Eigen::Matrix<double,19,19>::Identity());
      state.cov=posterior.topLeftCorner<18,18>();state.cov=0.5*(state.cov+state.cov.transpose()).eval();
      stats.updated=true;stats.rms=std::sqrt(final_system.squared_error/final_system.samples);
      camera.last_log_gain=log_gain;
    } else state=prior;
  } else state=prior;
  stats.log_gain=camera.last_log_gain;
  stats.optimize_ms=elapsedMs(optimization_start)-stats.prepare_ms;
  // Copy final accepted supports before insertion can evict reference landmarks.
  if(display&&stats.updated)for(const auto* p:final_landmarks)
    display->push_back({p->center,p->world,p->camera,false});
  const auto insert_start=SteadyClock::now();
  insert(id,image,state,points,geometry,std::exp(camera.last_log_gain),display);
  stats.insert_ms=elapsedMs(insert_start);
  stats.map_points=map_.size();return stats;
}

void VisualEstimator::insert(int id,const std::shared_ptr<Image>& image,const StatesGroup& state,
                             const std::vector<pointWithVar>& points,const VoxelMapManager& geometry,double gain,
                             std::vector<VisualPatch>* display) {
  const Camera& camera=cameras_[id];const Mat Rcw=camera.Rci*state.rot_end.transpose();
  const Vec tcw=camera.tci-Rcw*state.pos_end,origin=cameraCenter(camera,state);
  const double voxel=geometry.config_setting_.max_voxel_size_,spacing=voxel/4;
  int inserted=0;std::vector<Vec> inserted_directions;
  for(const auto& input:points) {
    const Vec p=input.point_w;
    VOXEL_LOCATION key(std::floor(p.x()/spacing),std::floor(p.y()/spacing),std::floor(p.z()/spacing));
    if(occupied_.count(key))continue;
    const Vec pc=Rcw*p+tcw;Eigen::Vector2d uv;double pitch;
    if(pc.norm()<1e-6||!image->locate(pc.normalized(),uv,pitch))continue;
    bool duplicate_support=false;
    for(const auto& b:inserted_directions)if((b-pc.normalized()).norm()<2*config_.radius){duplicate_support=true;break;}
    if(duplicate_support)continue;
    auto it=geometry.voxel_map_.find(VOXEL_LOCATION(std::floor(p.x()/voxel),std::floor(p.y()/voxel),std::floor(p.z()/voxel)));
    if(it==geometry.voxel_map_.end())continue;
    const VoxelOctoTree* cell=it->second;
    while(cell && !cell->plane_ptr_->is_plane_ && cell->octo_state_!=0) {
      const int index=(p.x()>cell->voxel_center_[0]?4:0)+(p.y()>cell->voxel_center_[1]?2:0)+(p.z()>cell->voxel_center_[2]?1:0);
      cell=cell->leaves_[index];
    }
    if(!cell||!cell->plane_ptr_->is_plane_)continue;
    const auto& plane=*cell->plane_ptr_;
    if(plane.normal_.squaredNorm()<0.9||plane.radius_<=0)continue;
    auto lm=std::make_unique<Landmark>();lm->camera=id;lm->center=p;lm->plane_center=plane.center_;
    lm->normal=plane.normal_;lm->plane_cov=plane.plane_var_;
    if(lm->normal.dot(origin-p)<0){lm->normal=-lm->normal;
      lm->plane_cov.topRightCorner<3,3>()*=-1;lm->plane_cov.bottomLeftCorner<3,3>()*=-1;}
    lm->plane_radius=plane.radius_;lm->roughness=std::max(0.0,double(plane.min_eigen_value_));
    lm->Rrw=Rcw;lm->trw=tcw;lm->gain=gain;lm->image=image;lm->key=key;
    const Mat q=Eigen::Quaterniond::FromTwoVectors(Vec::UnitZ(),pc.normalized()).toRotationMatrix();
    const Vec nr=Rcw*lm->normal;const double h=nr.dot(Rcw*plane.center_+tcw);
    bool valid=true;double minx=uv.x(),maxx=uv.x(),miny=uv.y(),maxy=uv.y(),energy=0;
    for(const Vec& t:template_) {
      const Vec br=q*t;Vec pr;
      if(!intersectPlane(br,nr,h,pr)){valid=false;break;}
      const Vec pw=Rcw.transpose()*(pr-tcw),dw=(pw-origin).normalized();
      if((pw-plane.center_).norm()+10*pitch*pr.norm()>3*plane.radius_){valid=false;break;}
      Eigen::Vector2d pixel;double local_pitch;
      if(!image->locate(br,pixel,local_pitch)){valid=false;break;}
      Sample sample;
      if(!image->sample(br,Mat::Identity()/std::pow(2.5*local_pitch,2),sample)){valid=false;break;}
      energy+=sample.gradient.squaredNorm()*local_pitch*local_pitch;
      lm->rays.push_back(br);lm->world.push_back(pw);
      lm->plane_jacobians.push_back(dw*rangePlaneJacobian(origin,dw,lm->normal,plane.center_,pr.norm()));
      minx=std::min(minx,pixel.x());maxx=std::max(maxx,pixel.x());miny=std::min(miny,pixel.y());maxy=std::max(maxy,pixel.y());
    }
    if(!valid||energy/config_.samples<image->noise()*image->noise())continue;
    if(!consistent(*lm,camera,state,points,geometry.config_setting_.dept_err_))continue;
    // Store only the native sensor tile needed by this reference; LUT is shared.
    constexpr int margin=64;
    const cv::Rect roi(int(std::floor(minx))-margin,int(std::floor(miny))-margin,
      int(std::ceil(maxx-minx))+2*margin+2,int(std::ceil(maxy-miny))+2*margin+2);
    lm->image=image->crop(roi);
    if(display)display->push_back({lm->center,lm->world,id,true});
    occupied_.insert(key);map_.push_back(std::move(lm));
    inserted_directions.push_back(pc.normalized());
    if(++inserted>=config_.max_patches)break;
  }
  // Capacity ownership only: no information/lifecycle/directional management.
  const size_t capacity=size_t(config_.max_patches)*32;
  while(map_.size()>capacity){occupied_.erase(map_.front()->key);map_.pop_front();}
}
} // namespace spherical
