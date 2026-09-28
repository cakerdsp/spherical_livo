#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <opencv2/core.hpp>
#include <array>
#include <memory>
#include <vector>

namespace spherical {
using Vec = Eigen::Vector3d;
using Mat = Eigen::Matrix3d;
using Basis = Eigen::Matrix<double, 3, 2>;

Mat skew(const Vec& x);
Basis tangentBasis(const Vec& b);
Mat rightJacobianInverse(const Vec& x);
std::vector<Vec> capTemplate(int count, double radius);
bool intersectPlane(const Vec& b, const Vec& n, double h, Vec& p);
Eigen::Matrix<double, 3, 6> bearingJacobian(const Mat& Rwi, const Vec& pwi,
                                          const Mat& Rci, const Vec& tci,
                                          const Vec& pw);

// Calibration builds the ray LUT and seeds its cell lookup.
// Runtime photometric sampling has no virtual camera or pixel warp.
struct CameraModel {
  enum class Type { Pinhole, KannalaBrandt, Mei };
  Type type = Type::Pinhole;
  int width = 0, height = 0;
  double fx = 0, fy = 0, cx = 0, cy = 0, xi = 0;
  std::array<double, 5> distortion{};
  Vec ray(double u, double v) const;
  bool project(const Vec& bearing, Eigen::Vector2d& pixel) const;
  void validate() const;
};

// Barycentric coordinates and their derivative with respect to a 3-D ray.
// Vertices are indices in the selected decimated sensor lattice.
struct Triangle {
  std::array<cv::Point,3> vertices;
  Vec weights;
  Mat derivative;
  Eigen::Vector2d pixel;
};

class RayAtlas {
 public:
  explicit RayAtlas(const CameraModel& model);
  int width() const { return width_; }
  int height() const { return height_; }
  Vec ray(int x, int y) const;
  float area(int x, int y) const;
  // Returns an original sensor cell only as an indexing aid. Containment is
  // checked with its three unit rays, including directions with z <= 0.
  bool locate(const Vec& b, Eigen::Vector2d& pixel, double& pitch) const;
  bool triangle(const Vec& b,int level,Triangle& result) const;
  int levels() const { return int(masses_.size()); }
  const cv::Mat& mass(int level) const { return masses_.at(level); }
  double supportRadius(int level,const cv::Point& pixel) const;
  cv::Rect supportBox(int level,const cv::Point& pixel) const;
 private:
  struct Seed { Eigen::Vector3f b; int pixel; };
  struct Node { int seed, left = -1, right = -1, axis = 0; };
  int build(std::vector<int>& ids, int begin, int end, int depth);
  void nearest(int node, const Vec& b, int& best, double& distance) const;
  CameraModel model_;
  bool locateCell(const Vec& b,int x,int y,Eigen::Vector2d& pixel) const;
  bool triangleCell(const Vec& b,int x,int y,int level,Triangle& result) const;
  int width_, height_, root_ = -1;
  std::vector<Eigen::Vector3f> rays_;
  std::vector<float> areas_;
  std::vector<Seed> seeds_;
  std::vector<Node> nodes_;
  std::vector<cv::Mat> masses_, radii_;
};

struct Sample {
  double value = 0;
  Vec gradient = Vec::Zero(); // derivative of the same interpolant on S2
};

// Frozen during each level's IEKF iterations. Blend adjacent prefiltered levels.
struct SamplingFilter {
  int level = 0;
  double blend = 0;
};

struct PlaneSupport {
  Vec normal, center;
  double height = 0, radius = 0;
  bool contains(const Vec& ray) const;
  bool containsKernel(const Vec& bearing,double kernel_radius) const;
};

class Image {
 public:
  Image(std::shared_ptr<const RayAtlas> atlas, const cv::Mat& gray);
  bool locate(const Vec& b, Eigen::Vector2d& pixel, double& pitch) const;
  // Select a conservative sensor-pyramid bandwidth from angular covariance.
  // Also check all intensity-independent support conditions before sampling.
  bool prepare(const Vec& b,const Mat& covariance,SamplingFilter& filter,
               double& cover,const PlaneSupport* support = nullptr) const;
  bool sample(const Vec& b,const SamplingFilter& filter,Sample& result,
              bool gradient = true,const PlaneSupport* support = nullptr) const;
  const RayAtlas& atlas() const { return *atlas_; }
 private:
  struct Level { cv::Mat values, valid; };
  bool inspect(const Vec& b,int level,Triangle& triangle,double& cover,
               const PlaneSupport* support) const;
  std::shared_ptr<const RayAtlas> atlas_;
  std::vector<Level> levels_;
};
} // namespace spherical
