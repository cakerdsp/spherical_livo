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

// Calibration is used only to build the original sensor's ray LUT.
// Runtime photometric sampling has no virtual camera or pixel warp.
struct CameraModel {
  enum class Type { Pinhole, KannalaBrandt, Mei };
  Type type = Type::Pinhole;
  int width = 0, height = 0;
  double fx = 0, fy = 0, cx = 0, cy = 0, xi = 0;
  std::array<double, 5> distortion{};
  Vec ray(double u, double v) const;
  void validate() const;
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
 private:
  struct Seed { Eigen::Vector3f b; int pixel; };
  struct Node { int seed, left = -1, right = -1, axis = 0; };
  int build(std::vector<int>& ids, int begin, int end, int depth);
  void nearest(int node, const Vec& b, int& best, double& distance) const;
  int width_, height_, root_ = -1;
  std::vector<Eigen::Vector3f> rays_;
  std::vector<float> areas_;
  std::vector<Seed> seeds_;
  std::vector<Node> nodes_;
};

struct Sample {
  double value = 0;
  Vec gradient = Vec::Zero(); // Riemannian gradient on the unit sphere.
  std::vector<std::pair<int, double>> weights;
};

class Image {
 public:
  Image(std::shared_ptr<const RayAtlas> atlas, const cv::Mat& gray);
  std::shared_ptr<Image> crop(const cv::Rect& roi) const;
  // precision is frozen for a linearization. The compact Wendland kernel has
  // zero value and derivative at its boundary; changing support is continuous.
  bool sample(const Vec& b, const Mat& precision, Sample& result) const;
  bool locate(const Vec& b, Eigen::Vector2d& pixel, double& pitch) const;
  const RayAtlas& atlas() const { return *atlas_; }
  double noise() const { return noise_; }
  cv::Rect bounds() const { return bounds_; }
 private:
  Image() = default;
  std::shared_ptr<const RayAtlas> atlas_;
  cv::Mat gray_;
  cv::Rect bounds_;
  double noise_ = 1.0 / 255.0;
};

double weightOverlap(const Sample& a, const Sample& b);
Mat precision(const Vec& b, const Eigen::Matrix2d& footprint);
} // namespace spherical
