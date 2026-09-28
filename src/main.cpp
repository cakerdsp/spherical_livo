#include "IMU_Processing.h"
#include "preprocess.h"
#include "spherical_vio.h"
#include <cv_bridge/cv_bridge.h>
#include <nav_msgs/Path.h>
#include <sensor_msgs/CompressedImage.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl_conversions/pcl_conversions.h>
#include <Eigen/Cholesky>
#include <boost/make_shared.hpp>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace {
using spherical::Mat;
using spherical::Vec;
template<class T>T parameter(const ros::NodeHandle& nh,const std::string& key,T value){nh.param<T>(key,value,value);return value;}
std::vector<double> required(const ros::NodeHandle& nh,const std::string& key,size_t count) {
  std::vector<double> values;
  if(!nh.getParam(key,values)||values.size()!=count)throw std::invalid_argument(key+": wrong/missing array size");
  for(double x:values)if(!std::isfinite(x))throw std::invalid_argument(key+": non-finite value");
  return values;
}
Mat rotation(const ros::NodeHandle& nh,const std::string& key) {
  const auto a=required(nh,key,9);Mat r;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)r(i,j)=a[3*i+j];
  if((r.transpose()*r-Mat::Identity()).norm()>1e-3||std::abs(r.determinant()-1)>1e-3)
    throw std::invalid_argument(key+": invalid rotation");
  return Eigen::Quaterniond(r).normalized().toRotationMatrix();
}
Vec translation(const ros::NodeHandle& nh,const std::string& key) {
  const auto a=required(nh,key,3);return Vec(a[0],a[1],a[2]);
}
spherical::CameraModel cameraModel(const ros::NodeHandle& nh,const std::string& ns) {
  spherical::CameraModel c;const auto model=parameter<std::string>(nh,ns+"/model","");
  if(model=="Pinhole")c.type=spherical::CameraModel::Type::Pinhole;
  else if(model=="EquidistantCamera")c.type=spherical::CameraModel::Type::KannalaBrandt;
  else if(model=="MEI")c.type=spherical::CameraModel::Type::Mei;
  else throw std::invalid_argument(ns+": supported models are Pinhole, EquidistantCamera, MEI");
  c.width=parameter(nh,ns+"/width",0);c.height=parameter(nh,ns+"/height",0);
  c.fx=parameter(nh,ns+"/fx",0.0);c.fy=parameter(nh,ns+"/fy",0.0);
  c.cx=parameter(nh,ns+"/cx",0.0);c.cy=parameter(nh,ns+"/cy",0.0);c.xi=parameter(nh,ns+"/xi",0.0);
  if(parameter(nh,ns+"/scale",1.0)!=1.0)throw std::invalid_argument("use native image dimensions; camera scale must be 1");
  if(c.type==spherical::CameraModel::Type::KannalaBrandt)
    for(int i=0;i<4;++i)c.distortion[i]=parameter(nh,ns+"/k"+std::to_string(i+1),0.0);
  else if(c.type==spherical::CameraModel::Type::Mei){
    c.distortion[0]=parameter(nh,ns+"/k1",0.0);c.distortion[1]=parameter(nh,ns+"/k2",0.0);
    c.distortion[2]=parameter(nh,ns+"/p1",0.0);c.distortion[3]=parameter(nh,ns+"/p2",0.0);
  }else for(int i=0;i<5;++i)c.distortion[i]=parameter(nh,ns+"/d"+std::to_string(i),0.0);
  c.validate();return c;
}

struct Frame {double time;int camera;cv::Mat gray;};
struct TimedPoint {double time;PointType point;};

class Node {
 public:
  explicit Node(ros::NodeHandle nh):nh_(nh) {
    Ril_=rotation(nh_,"extrin_calib/extrinsic_R");til_=translation(nh_,"extrin_calib/extrinsic_T");
    VoxelMapConfig config;loadVoxelConfig(nh_,config);
    if(config.max_layer_<0||config.max_layer_>=int(config.layer_init_num_.size())||config.max_voxel_size_<=0||
       config.dept_err_<=0||config.beam_err_<=0)throw std::invalid_argument("invalid LiDAR map/noise configuration");
    std::unordered_map<VOXEL_LOCATION,VoxelOctoTree*> empty;
    geometry_=std::make_unique<VoxelMapManager>(config,empty);geometry_->extR_=Ril_;geometry_->extT_=til_;
    imu_.set_extrinsic(til_,Ril_);
    imu_.set_gyr_cov_scale(Vec::Constant(parameter(nh_,"imu/gyr_cov",0.1)));
    imu_.set_acc_cov_scale(Vec::Constant(parameter(nh_,"imu/acc_cov",0.1)));
    imu_.set_gyr_bias_cov(Vec::Constant(parameter(nh_,"imu/b_gyr_cov",0.0001)));
    imu_.set_acc_bias_cov(Vec::Constant(parameter(nh_,"imu/b_acc_cov",0.0001)));
    imu_.set_imu_init_frame_num(parameter(nh_,"imu/imu_int_frame",100));
    imu_offset_=parameter(nh_,"time_offset/imu_time_offset",0.0);
    lidar_offset_=parameter(nh_,"time_offset/lidar_time_offset",0.0);
    pre_.set(false,parameter(nh_,"preprocess/lidar_type",2),parameter(nh_,"preprocess/blind",0.5),
             parameter(nh_,"preprocess/point_filter_num",1));
    pre_.N_SCANS=parameter(nh_,"preprocess/scan_line",32);pre_.blind_sqr=pre_.blind*pre_.blind;
    imu_.lidar_type=pre_.lidar_type;
    const double leaf=parameter(nh_,"preprocess/filter_size_surf",0.1);
    if(leaf<=0||pre_.N_SCANS<1||pre_.N_SCANS>128||pre_.point_filter_num<1||pre_.lidar_type<1||pre_.lidar_type>8)
      throw std::invalid_argument("invalid preprocessing configuration");
    downsample_.setLeafSize(leaf,leaf,leaf);
    spherical::VisualConfig visual;
    visual.radius=parameter(nh_,"visual/patch_radius_deg",0.7)*M_PI/180;
    visual.samples=parameter(nh_,"visual/patch_samples",48);
    visual.max_patches=parameter(nh_,"visual/max_patches",120);
    const int count=parameter(nh_,"common/num_cameras",1);
    if(count<1)throw std::invalid_argument("at least one camera is required");
    offsets_.resize(count);last_camera_.assign(count,-1);
    std::vector<spherical::Camera> cameras;
    for(int i=0;i<count;++i) {
      const std::string key="cameras/camera"+std::to_string(i);
      const auto ns=parameter<std::string>(nh_,key+"/camera_ns","camera"+std::to_string(i));
      spherical::Camera camera;camera.atlas=std::make_shared<spherical::RayAtlas>(cameraModel(nh_,ns));
      const Mat Rcl=rotation(nh_,key+"/Rcl");const Vec tcl=translation(nh_,key+"/Pcl");
      camera.Rci=Rcl*Ril_.transpose();camera.tci=tcl-camera.Rci*til_;
      cameras.push_back(camera);offsets_[i]=parameter(nh_,key+"/time_offset",0.0);
      const auto topic=parameter<std::string>(nh_,key+"/img_topic","");
      if(topic.empty())throw std::invalid_argument(key+"/img_topic is required");
      if(topic.size()>=11&&topic.substr(topic.size()-11)=="/compressed") {
        images_.push_back(nh_.subscribe<sensor_msgs::CompressedImage>(topic,8,[this,i](const sensor_msgs::CompressedImage::ConstPtr& msg){
          receiveImage(i,msg->header.stamp.toSec(),cv_bridge::toCvCopy(msg,"mono8")->image);
        }));
      }else{
        images_.push_back(nh_.subscribe<sensor_msgs::Image>(topic,8,[this,i](const sensor_msgs::Image::ConstPtr& msg){
          receiveImage(i,msg->header.stamp.toSec(),cv_bridge::toCvCopy(msg,"mono8")->image);
        }));
      }
    }
    visual_=std::make_unique<spherical::VisualEstimator>(visual,std::move(cameras));
    const auto lidar_topic=parameter<std::string>(nh_,"common/lid_topic","/velodyne_points");
    if(pre_.lidar_type==AVIA)lidar_=nh_.subscribe<livox_ros_driver2::CustomMsg>(lidar_topic,16,[this](const livox_ros_driver2::CustomMsg::ConstPtr& msg){
      PointCloudXYZI::Ptr cloud(new PointCloudXYZI);pre_.process(msg,cloud);receiveCloud(msg->header.stamp.toSec(),cloud);
    });
    else lidar_=nh_.subscribe<sensor_msgs::PointCloud2>(lidar_topic,16,[this](const sensor_msgs::PointCloud2::ConstPtr& msg){
      PointCloudXYZI::Ptr cloud(new PointCloudXYZI);pre_.process(msg,cloud);receiveCloud(msg->header.stamp.toSec(),cloud);
    });
    imu_sub_=nh_.subscribe<sensor_msgs::Imu>(parameter<std::string>(nh_,"common/imu_topic","/imu"),4000,[this](const sensor_msgs::Imu::ConstPtr& msg){
      sensor_msgs::Imu::Ptr corrected(new sensor_msgs::Imu(*msg));corrected->header.stamp.fromSec(msg->header.stamp.toSec()+imu_offset_);
      const double t=corrected->header.stamp.toSec();
      if(t<=last_imu_time_)throw std::runtime_error("IMU timestamps moved backwards or repeated; restart for a new bag");
      last_imu_time_=t;imu_queue_.push_back(corrected);
    });
    odom_=nh_.advertise<nav_msgs::Odometry>("odometry",10);
    path_pub_=nh_.advertise<nav_msgs::Path>("path",1);
    cloud_pub_=nh_.advertise<sensor_msgs::PointCloud2>("cloud",1);
    geometry_->voxel_map_pub_=nh_.advertise<visualization_msgs::MarkerArray>("planes",1);
    const auto dir=parameter<std::string>(nh_,"output_dir","");
    if(!dir.empty()){
      std::filesystem::create_directories(dir);trajectory_.open(dir+"/trajectory.txt");metrics_.open(dir+"/visual.csv");
      if(!trajectory_||!metrics_)throw std::runtime_error("cannot open output files");
      metrics_<<"time,camera,candidates,depth_rejected,patches,cross_camera,map_points,updated,rms,log_gain,seconds\n";
    }
  }

  ~Node(){for(auto& kv:geometry_->voxel_map_)release(kv.second);}

  void processReady() {
    if(imu_queue_.size()>20000||points_.size()>2000000||frames_.size()>256)
      throw std::runtime_error("sensor buffers exceeded capacity; check timestamps/topics/input rates");
    if(points_.empty()||imu_queue_.size()<2||scan_ends_.empty())return;
    if(time_<0){time_=std::max(points_.front().time,imu_queue_.front()->header.stamp.toSec());imu_.first_lidar_time=time_;}
    while(!scan_ends_.empty()&&scan_ends_.front()<=time_+1e-9)scan_ends_.pop_front();
    while(!frames_.empty()&&frames_.front().time<time_-1e-9){ROS_WARN_THROTTLE(1,"Dropping late image; no timestamp substitution");frames_.pop_front();}
    for(int work=0;work<32;++work) {
      double target=frames_.empty()?std::numeric_limits<double>::infinity():frames_.front().time;
      // One native scan of reorder buffering when no image is queued. A missing
      // camera never blocks LIO or the other cameras.
      if(!scan_ends_.empty()&&scan_ends_.front()<=latest_scan_begin_)target=std::min(target,scan_ends_.front());
      if(!std::isfinite(target)||target>latest_lidar_end_||target>last_imu_time_||target<time_-1e-9)break;
      if(target>time_+1e-9)advance(target);
      while(!frames_.empty()&&std::abs(frames_.front().time-target)<=1e-9) {
        Frame frame=std::move(frames_.front());frames_.pop_front();
        if(imu_.imu_need_init||!map_ready_)continue;
        const auto start=std::chrono::steady_clock::now();
        const auto stats=visual_->process(frame.camera,frame.gray,state_,recent_points_,*geometry_);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(metrics_)metrics_<<std::setprecision(12)<<target<<','<<frame.camera<<','<<stats.candidates<<','<<stats.depth_rejected<<','
          <<stats.patches<<','<<stats.cross_camera<<','<<stats.map_points<<','<<stats.updated<<','<<stats.rms<<','<<stats.log_gain<<','<<seconds<<'\n';
        ROS_INFO_THROTTLE(1,"Spherical VIO: camera=%d patches=%d cross_camera=%d map=%zu %.3fs",frame.camera,stats.patches,stats.cross_camera,stats.map_points,seconds);
      }
      while(!scan_ends_.empty()&&scan_ends_.front()<=time_+1e-9)scan_ends_.pop_front();
      if(!imu_.imu_need_init)publish();
      if(points_.empty()&&frames_.empty())break;
    }
    if(imu_queue_.size()>20000||points_.size()>2000000||frames_.size()>256)
      throw std::runtime_error("sensor buffers exceeded capacity; check timestamps/topics/input rates");
  }

 private:
  static void release(VoxelOctoTree* cell){delete cell;}
  void receiveImage(int camera,double stamp,const cv::Mat& gray){
    const double t=stamp+offsets_[camera];
    if(t<=last_camera_[camera])throw std::runtime_error("camera timestamps repeated or moved backwards");
    last_camera_[camera]=t;
    auto pos=std::upper_bound(frames_.begin(),frames_.end(),t,[](double a,const Frame& b){return a<b.time;});
    frames_.insert(pos,Frame{t,camera,gray.clone()});
  }
  void receiveCloud(double stamp,const PointCloudXYZI::Ptr& cloud){
    stamp+=lidar_offset_;if(stamp<=latest_scan_begin_)throw std::runtime_error("LiDAR timestamps repeated or moved backwards");
    latest_scan_begin_=stamp;
    double end=stamp;
    for(const auto& p:cloud->points){
      if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||!std::isfinite(p.curvature)||p.curvature<0)continue;
      const double t=stamp+p.curvature*0.001;if(t>time_)points_.push_back({t,p});end=std::max(end,t);
    }
    std::stable_sort(points_.begin(),points_.end(),[](const TimedPoint& a,const TimedPoint& b){return a.time<b.time;});
    if(end>stamp){scan_ends_.push_back(end);latest_lidar_end_=std::max(latest_lidar_end_,end);}
  }
  sensor_msgs::Imu::Ptr interpolated(double target,const sensor_msgs::Imu& a,const sensor_msgs::Imu& b){
    const double span=b.header.stamp.toSec()-a.header.stamp.toSec();
    const double w=span>0?(target-a.header.stamp.toSec())/span:0;
    auto out=boost::make_shared<sensor_msgs::Imu>(a);out->header.stamp.fromSec(target);
    auto mix=[w](double x,double y){return x+(y-x)*w;};
    out->angular_velocity.x=mix(a.angular_velocity.x,b.angular_velocity.x);out->angular_velocity.y=mix(a.angular_velocity.y,b.angular_velocity.y);
    out->angular_velocity.z=mix(a.angular_velocity.z,b.angular_velocity.z);
    out->linear_acceleration.x=mix(a.linear_acceleration.x,b.linear_acceleration.x);out->linear_acceleration.y=mix(a.linear_acceleration.y,b.linear_acceleration.y);
    out->linear_acceleration.z=mix(a.linear_acceleration.z,b.linear_acceleration.z);return out;
  }
  void advance(double target){
    LidarMeasureGroup measurement;measurement.lio_vio_flg=LIO;measurement.last_lio_update_time=time_;
    MeasureGroup group;group.lio_time=target;
    while(imu_queue_.size()>1&&imu_queue_[1]->header.stamp.toSec()<=time_)imu_queue_.pop_front();
    if(imu_queue_.size()<2||imu_queue_.front()->header.stamp.toSec()>time_)
      throw std::runtime_error("no IMU bracket at propagation start");
    group.imu.push_back(interpolated(time_,*imu_queue_[0],*imu_queue_[1]));
    while(imu_queue_.size()>1&&imu_queue_[1]->header.stamp.toSec()<target){
      imu_queue_.pop_front();group.imu.push_back(imu_queue_.front());
    }
    if(imu_queue_.size()<2)throw std::runtime_error("no IMU bracket at propagation end");
    group.imu.push_back(interpolated(target,*imu_queue_[0],*imu_queue_[1]));
    const auto &gyro = group.imu.back()->angular_velocity;
    measured_gyro_ = Vec(gyro.x, gyro.y, gyro.z);
    while(!points_.empty()&&points_.front().time<=target){
      auto point=points_.front();points_.pop_front();
      if(point.time<time_)continue;point.point.curvature=(point.time-time_)*1000;
      measurement.pcl_proc_cur->push_back(point.point);
    }
    measurement.measures.push_back(group);PointCloudXYZI::Ptr undistorted(new PointCloudXYZI);
    imu_.Process2(measurement,state_,undistorted);time_=target;
    if(imu_.imu_need_init||undistorted->size()<6)return;
    PointCloudXYZI::Ptr down(new PointCloudXYZI),world(new PointCloudXYZI);
    downsample_.setInputCloud(undistorted);downsample_.filter(*down);if(down->size()<6)return;
    geometry_->state_=state_;geometry_->feats_down_body_=down;geometry_->feats_down_size_=down->size();
    for(const auto& p:down->points){PointType q=p;const Vec w=state_.rot_end*(Ril_*Vec(p.x,p.y,p.z)+til_)+state_.pos_end;
      q.x=w.x();q.y=w.y();q.z=w.z();world->push_back(q);}
    geometry_->feats_down_world_=world;
    if(!map_ready_){geometry_->BuildVoxelMap();map_ready_=true;}
    StatesGroup prior=state_;geometry_->StateEstimation(prior);state_=geometry_->state_;
    state_.cov=0.5*(state_.cov+state_.cov.transpose()).eval();
    world->clear();
    for(size_t i=0;i<down->size();++i){
      auto& pv=geometry_->pv_list_[i];const auto& p=down->points[i];const Vec pi=Ril_*Vec(p.x,p.y,p.z)+til_;
      pv.point_w=state_.rot_end*pi+state_.pos_end;
      Eigen::Matrix<double,3,18> j=Eigen::Matrix<double,3,18>::Zero();
      j.leftCols<3>()=-state_.rot_end*spherical::skew(pi);j.block<3,3>(0,3)=Mat::Identity();
      pv.var_nostate=state_.rot_end*Ril_*geometry_->body_cov_list_[i]*Ril_.transpose()*state_.rot_end.transpose();
      pv.var=pv.var_nostate+j*state_.cov*j.transpose();
      PointType q=p;q.x=pv.point_w.x();q.y=pv.point_w.y();q.z=pv.point_w.z();world->push_back(q);
    }
    geometry_->UpdateVoxelMap(geometry_->pv_list_);recent_points_=geometry_->pv_list_;
    if(geometry_->config_setting_.map_sliding_en)geometry_->mapSliding();
    if(geometry_->config_setting_.is_pub_plane_map_)geometry_->pubVoxelMap();
    sensor_msgs::PointCloud2 msg;pcl::toROSMsg(*world,msg);msg.header.frame_id="world";msg.header.stamp.fromSec(time_);cloud_pub_.publish(msg);
  }
  void publish(){
    nav_msgs::Odometry odom;odom.header.frame_id="world";odom.child_frame_id="imu";odom.header.stamp.fromSec(time_);
    odom.pose.pose.position.x=state_.pos_end.x();odom.pose.pose.position.y=state_.pos_end.y();odom.pose.pose.position.z=state_.pos_end.z();
    const Eigen::Quaterniond q(state_.rot_end);odom.pose.pose.orientation.x=q.x();odom.pose.pose.orientation.y=q.y();
    odom.pose.pose.orientation.z=q.z();odom.pose.pose.orientation.w=q.w();
    // ROS pose covariance uses world translation then fixed-axis orientation.
    Eigen::Matrix<double,6,18> j=Eigen::Matrix<double,6,18>::Zero();j.block<3,3>(0,3)=Mat::Identity();j.block<3,3>(3,0)=state_.rot_end;
    const Eigen::Matrix<double,6,6> covariance=j*state_.cov*j.transpose();
    for(int r=0;r<6;++r)for(int c=0;c<6;++c)odom.pose.covariance[6*r+c]=covariance(r,c);
    const Vec vi=state_.rot_end.transpose()*state_.vel_end;
    odom.twist.twist.linear.x=vi.x();odom.twist.twist.linear.y=vi.y();odom.twist.twist.linear.z=vi.z();
    const Vec omega=measured_gyro_-state_.bias_g;
    odom.twist.twist.angular.x=omega.x();odom.twist.twist.angular.y=omega.y();odom.twist.twist.angular.z=omega.z();
    Eigen::Matrix<double,6,18> jt=Eigen::Matrix<double,6,18>::Zero();
    jt.block<3,3>(0,0)=spherical::skew(vi);jt.block<3,3>(0,6)=state_.rot_end.transpose();
    jt.block<3,3>(3,9)=-Mat::Identity();
    Eigen::Matrix<double,6,6> ct=jt*state_.cov*jt.transpose();
    ct.bottomRightCorner<3,3>()+=imu_.cov_gyr.asDiagonal();
    for(int r=0;r<6;++r)for(int c=0;c<6;++c)odom.twist.covariance[6*r+c]=ct(r,c);
    odom_.publish(odom);
    geometry_msgs::PoseStamped pose;pose.header=odom.header;pose.pose=odom.pose.pose;
    path_.header=odom.header;path_.poses.push_back(pose);if(path_.poses.size()>20000)path_.poses.erase(path_.poses.begin(),path_.poses.begin()+1000);
    path_pub_.publish(path_);
    if(trajectory_)trajectory_<<std::fixed<<std::setprecision(9)<<time_<<' '<<state_.pos_end.transpose()<<' '<<q.x()<<' '<<q.y()<<' '<<q.z()<<' '<<q.w()<<'\n';
  }
  ros::NodeHandle nh_;ros::Subscriber lidar_,imu_sub_;std::vector<ros::Subscriber> images_;
  ros::Publisher odom_,cloud_pub_,path_pub_;nav_msgs::Path path_;
  Preprocess pre_;ImuProcess imu_;StatesGroup state_;Mat Ril_;Vec til_;Vec measured_gyro_=Vec::Zero();
  std::unique_ptr<VoxelMapManager> geometry_;std::unique_ptr<spherical::VisualEstimator> visual_;
  pcl::VoxelGrid<PointType> downsample_;std::deque<Frame> frames_;std::deque<TimedPoint> points_;
  std::deque<sensor_msgs::Imu::ConstPtr> imu_queue_;std::deque<double> scan_ends_;
  std::vector<pointWithVar> recent_points_;std::vector<double> offsets_,last_camera_;
  double time_=-1,last_imu_time_=-1,latest_scan_begin_=-1,latest_lidar_end_=-1,imu_offset_=0,lidar_offset_=0;
  bool map_ready_=false;std::ofstream trajectory_,metrics_;
};
}

int main(int argc,char** argv){
  ros::init(argc,argv,"spherical_livo");ros::NodeHandle nh("~");
  try{Node node(nh);ROS_INFO("Spherical LIVO ready: camera LUTs built, sensor subscribers ready");
    ros::Rate rate(500);while(ros::ok()){ros::spinOnce();node.processReady();rate.sleep();}}
  catch(const std::exception& e){ROS_FATAL("Spherical LIVO stopped: %s",e.what());return 1;}
  return 0;
}
