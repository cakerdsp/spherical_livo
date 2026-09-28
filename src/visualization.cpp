#include "visualization.h"
#include <cv_bridge/cv_bridge.h>
#include <opencv2/imgproc.hpp>
#include <pcl_conversions/pcl_conversions.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>

namespace spherical {
Visualization::Visualization(ros::NodeHandle& nh,const std::vector<Camera>& cameras):cameras_(cameras) {
  cloud_=nh.advertise<sensor_msgs::PointCloud2>("cloud_rgb",1);
  for(size_t i=0;i<cameras.size();++i)
    images_.push_back(nh.advertise<sensor_msgs::Image>("image/camera_"+std::to_string(i),1));
}

bool Visualization::wantsImage(int camera) const {return images_.at(camera).getNumSubscribers()>0;}

void Visualization::publish(int id,double stamp,const cv::Mat& bgr,const StatesGroup& state,
                            const std::vector<pointWithVar>& points,const VisualStats& stats,
                            const std::vector<VisualPatch>& patches,bool initialized) {
  const Camera& camera=cameras_.at(id);
  const Mat Rcw=camera.Rci*state.rot_end.transpose();
  const Vec tcw=camera.tci-Rcw*state.pos_end;
  // No positive-z gate: a calibrated fisheye may see beyond a hemisphere.
  const auto project=[&](const Vec& world,Eigen::Vector2d& uv) {
    const Vec pc=Rcw*world+tcw;double pitch;
    return pc.allFinite()&&pc.norm()>1e-6&&camera.atlas->locate(pc.normalized(),uv,pitch)&&
      uv.allFinite()&&uv.x()>=0&&uv.y()>=0&&uv.x()<bgr.cols-1&&uv.y()<bgr.rows-1;
  };
  if(initialized&&cloud_.getNumSubscribers()>0) {
    std::vector<pcl::PointXYZRGB> rgb(points.size());
    // Unlike vector<bool>, separate bytes are independent memory locations.
    std::vector<uint8_t> valid(points.size(),0);
    const int workers=std::max(1,std::min(4,omp_get_max_threads()));
    const long long count=static_cast<long long>(points.size());
    // cake_slam's pattern: independent slots, then deterministic serial compaction.
    // All image, atlas, pose and point data stay read-only inside this loop.
#pragma omp parallel for schedule(static) num_threads(workers) if(count>=1024)
    for(long long i=0;i<count;++i) {
      Eigen::Vector2d uv;const Vec& world=points[i].point_w;
      if(!project(world,uv))continue;
      const int x=static_cast<int>(uv.x()),y=static_cast<int>(uv.y());
      const double dx=uv.x()-x,dy=uv.y()-y;
      cv::Vec3d color(0,0,0);
      for(int v=0;v<2;++v)for(int u=0;u<2;++u) {
        const double w=(u?dx:1-dx)*(v?dy:1-dy);
        const cv::Vec3b pixel=bgr.at<cv::Vec3b>(y+v,x+u);
        for(int c=0;c<3;++c)color[c]+=w*pixel[c];
      }
      auto& p=rgb[i];p.x=world.x();p.y=world.y();p.z=world.z();
      p.r=cv::saturate_cast<uint8_t>(color[2]);p.g=cv::saturate_cast<uint8_t>(color[1]);
      p.b=cv::saturate_cast<uint8_t>(color[0]);p.a=255;valid[i]=1;
    }
    pcl::PointCloud<pcl::PointXYZRGB> colored;colored.reserve(points.size());
    for(size_t i=0;i<rgb.size();++i)if(valid[i])colored.push_back(rgb[i]);
    sensor_msgs::PointCloud2 msg;pcl::toROSMsg(colored,msg);
    msg.header.frame_id="world";msg.header.stamp.fromSec(stamp);cloud_.publish(msg);
  }
  if(!wantsImage(id))return;
  cv::Mat overlay=bgr.clone();
  for(const auto& patch:patches) {
    const cv::Scalar color=patch.inserted?cv::Scalar(0,220,255):
      (patch.reference_camera==id?cv::Scalar(0,255,0):cv::Scalar(255,0,255));
    Eigen::Vector2d uv;
    for(const Vec& world:patch.samples)if(project(world,uv))
      cv::circle(overlay,cv::Point(cvRound(uv.x()),cvRound(uv.y())),1,color,-1,cv::LINE_AA);
    if(project(patch.center,uv))
      cv::circle(overlay,cv::Point(cvRound(uv.x()),cvRound(uv.y())),4,color,1,cv::LINE_AA);
  }
  std::ostringstream label;label<<"camera "<<id<<" patches "<<stats.patches<<" cross "<<stats.cross_camera
    <<(initialized?(stats.updated?" | updated":" | no visual update"):" | initializing");
  const auto text=[&](const std::string& s,int y) {
    cv::putText(overlay,s,{12,y},cv::FONT_HERSHEY_SIMPLEX,0.5,{0,0,0},3,cv::LINE_AA);
    cv::putText(overlay,s,{12,y},cv::FONT_HERSHEY_SIMPLEX,0.5,{255,255,255},1,cv::LINE_AA);
  };
  text(label.str(),24);text("green: tracked  magenta: cross-camera  yellow: new reference",46);
  std_msgs::Header header;header.stamp.fromSec(stamp);header.frame_id="camera_"+std::to_string(id);
  images_[id].publish(cv_bridge::CvImage(header,"bgr8",overlay).toImageMsg());
}
} // namespace spherical
