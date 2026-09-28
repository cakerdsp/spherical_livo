#include "spherical_vio.h"
#include <gtest/gtest.h>
#include <ros/ros.h>

namespace {
using spherical::Vec;using spherical::Mat;
spherical::CameraModel model(){spherical::CameraModel m;m.width=160;m.height=128;m.fx=m.fy=120;m.cx=79.5;m.cy=63.5;return m;}
cv::Mat render(const spherical::CameraModel& m,const Vec& camera_origin,double gain=1){
  cv::Mat raw(m.height,m.width,CV_8UC1);
  for(int y=0;y<m.height;++y)for(int x=0;x<m.width;++x){const Vec b=m.ray(x,y);const Vec p=camera_origin+b*((4-camera_origin.z())/b.z());
    const double value=(125+40*std::sin(17*p.x())+35*std::cos(13*p.y())+15*std::sin(9*(p.x()+p.y())))/gain;
    raw.at<uchar>(y,x)=cv::saturate_cast<uchar>(value);}
  return raw;
}
void release(VoxelOctoTree* c){delete c;}
}

TEST(SphericalIntegration,SharedMapReusesOtherCameraReferenceAndReducesPoseError){
  if(!ros::isInitialized()){ros::M_string remappings;ros::init(remappings,"spherical_integration_test",ros::init_options::NoSigintHandler);}
  auto m=model();auto atlas=std::make_shared<spherical::RayAtlas>(m);
  spherical::Camera a,b;a.atlas=b.atlas=atlas;b.tci=Vec(-0.2,0,0);
  spherical::VisualConfig visual;visual.radius=0.025;visual.samples=24;visual.max_patches=32;
  spherical::VisualEstimator estimator(visual,{a,b});
  VoxelMapConfig config{};config.max_voxel_size_=1;config.max_layer_=0;config.max_iterations_=5;
  config.layer_init_num_={5};config.max_points_num_=100;config.planner_threshold_=0.001;
  config.beam_err_=0.02;config.dept_err_=0.01;config.sigma_num_=3;
  std::unordered_map<VOXEL_LOCATION,VoxelOctoTree*> empty;VoxelMapManager geometry(config,empty);
  geometry.extR_=Mat::Identity();geometry.extT_=Vec::Zero();geometry.state_=StatesGroup();
  for(double y=-2;y<=2;y+=0.1)for(double x=-2.5;x<=2.5;x+=0.1){
    PointType p{};p.x=x;p.y=y;p.z=4;geometry.feats_down_body_->push_back(p);geometry.feats_down_world_->push_back(p);
  }
  geometry.BuildVoxelMap();std::vector<pointWithVar> points;
  for(const auto& p:geometry.feats_down_world_->points){pointWithVar v;v.point_w=Vec(p.x,p.y,p.z);v.var=Mat::Identity()*1e-5;points.push_back(v);}
  StatesGroup state;state.cov=Eigen::Matrix<double,18,18>::Identity()*0.01;
  const auto seed=estimator.process(0,render(m,Vec::Zero()),state,points,geometry);
  ASSERT_GT(seed.map_points,0u);
  state.pos_end=Vec(0.012,-0.008,0.006);const double before=state.pos_end.norm();
  const auto tracked=estimator.process(1,render(m,Vec(0.2,0,0)),state,points,geometry);
  EXPECT_TRUE(tracked.updated);EXPECT_GT(tracked.cross_camera,0);
  EXPECT_LT(state.pos_end.norm(),before);
  Eigen::LLT<Eigen::Matrix<double,18,18>> llt(state.cov);EXPECT_EQ(llt.info(),Eigen::Success);
  for(auto& kv:geometry.voxel_map_)release(kv.second);
}
