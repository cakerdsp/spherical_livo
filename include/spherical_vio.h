#pragma once

#include "spherical.h"
#include "voxel_map.h"
#include <deque>
#include <memory>
#include <unordered_set>

namespace spherical {
struct Camera {
  std::shared_ptr<RayAtlas> atlas;
  Mat Rci = Mat::Identity();
  Vec tci = Vec::Zero();
  double last_log_gain = 0;
};

struct VisualConfig {
  double radius = 0.012;
  int samples = 48;
  int max_patches = 120;
};

struct VisualPatch {
  Vec center;
  std::vector<Vec> samples;
  int reference_camera = 0;
  bool inserted = false;
};

struct VisualStats {
  size_t map_points = 0;
  int candidates = 0, patches = 0, cross_camera = 0, depth_rejected = 0;
  bool updated = false;
  double rms = 0, log_gain = 0;
  double candidates_ms = 0, prepare_ms = 0, optimize_ms = 0, insert_ms = 0;
};

class VisualEstimator {
 public:
  VisualEstimator(VisualConfig config, std::vector<Camera> cameras);
  VisualStats process(int camera, const cv::Mat& gray, StatesGroup& state,
                      const std::vector<pointWithVar>& points, const VoxelMapManager& geometry,
                      std::vector<VisualPatch>* display = nullptr);
 private:
  struct Landmark {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    int camera;
    Vec center, normal, plane_center;
    Eigen::Matrix<double,6,6> plane_cov;
    double plane_radius, roughness, gain;
    Mat Rrw;
    Vec trw;
    std::shared_ptr<Image> image;
    std::vector<Vec> rays, world;
    std::vector<Eigen::Matrix<double,3,6>> plane_jacobians;
    VOXEL_LOCATION key{0,0,0};
  };
  struct Track {
    const Landmark* point;
    std::vector<Mat> precision;
    std::vector<Sample> reference;
  };
  struct Linearization {
    Eigen::Matrix<double,19,19> H = Eigen::Matrix<double,19,19>::Zero();
    Eigen::Matrix<double,19,1> g = Eigen::Matrix<double,19,1>::Zero();
    std::vector<Eigen::MatrixXd> whitening;
    double cost = 0, squared_error = 0;
    int samples = 0;
  };
  bool consistent(const Landmark& lm, const Camera& camera, const StatesGroup& state,
                  const std::vector<pointWithVar>& points, double range_noise, double full_cover = 0) const;
  bool prepare(const Landmark& lm, const Image& image, const Camera& camera,
               const StatesGroup& state, int scale, Track& track) const;
  bool linearize(const std::vector<Track>& tracks, const Image& image, const Camera& camera,
                 const StatesGroup& state, double log_gain, const StatesGroup& prior,
                 const Eigen::Matrix<double,18,18>& information, Linearization& result,
                 const std::vector<Eigen::MatrixXd>* fixed_whitening = nullptr) const;
  void insert(int camera_id, const std::shared_ptr<Image>& image, const StatesGroup& state,
              const std::vector<pointWithVar>& points, const VoxelMapManager& geometry, double gain,
              std::vector<VisualPatch>* display);
  VisualConfig config_;
  std::vector<Camera> cameras_;
  std::vector<Vec> template_;
  std::deque<std::unique_ptr<Landmark>> map_;
  std::unordered_set<VOXEL_LOCATION> occupied_;
};
} // namespace spherical
