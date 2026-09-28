#pragma once
#include "spherical.h"
#include <vikit/abstract_camera.h>

class SphericalCamera final : public vk::AbstractCamera
{
public:
  spherical::CameraModel model;
  std::shared_ptr<const spherical::RayAtlas> atlas;
  explicit SphericalCamera(const spherical::CameraModel& m)
    : vk::AbstractCamera(m.width,m.height,1.0), model(m),
      atlas(std::make_shared<spherical::RayAtlas>(m)) {}
  Eigen::Vector3d cam2world(const double& u,const double& v) const override { return model.ray(u,v).normalized(); }
  Eigen::Vector3d cam2world(const Eigen::Vector2d& uv) const override { return cam2world(uv.x(),uv.y()); }
  Eigen::Vector2d world2cam(const Eigen::Vector3d& p) const override {
    Eigen::Vector2d uv(-10000,-10000); double pitch;
    if(p.allFinite() && p.squaredNorm()>1e-12) atlas->locate(p.normalized(),uv,pitch);
    return uv;
  }
  Eigen::Vector2d world2cam(const Eigen::Vector2d& uv) const override { return world2cam(Eigen::Vector3d(uv.x(),uv.y(),1)); }
  double fx() const override { return model.fx; }
  double fy() const override { return model.fy; }
  double cx() const override { return model.cx; }
  double cy() const override { return model.cy; }
  double errorMultiplier() const override { return model.fx; }
  double errorMultiplier2() const override { return model.fx*model.fy; }
};
