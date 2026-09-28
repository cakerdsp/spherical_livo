#pragma once

#include "spherical_vio.h"
#include <ros/ros.h>

namespace spherical {
// Display-only projections use the same native-camera ray atlas as the solver.
class Visualization {
 public:
  Visualization(ros::NodeHandle& nh,const std::vector<Camera>& cameras);
  bool wantsImage(int camera) const;
  void publish(int camera,double stamp,const cv::Mat& bgr,const StatesGroup& state,
               const std::vector<pointWithVar>& points,const VisualStats& stats,
               const std::vector<VisualPatch>& patches,bool initialized=true);
 private:
  std::vector<Camera> cameras_;
  std::vector<ros::Publisher> images_;
  ros::Publisher cloud_;
};
} // namespace spherical
