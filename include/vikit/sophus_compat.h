#ifndef VIKIT_SOPHUS_COMPAT_H
#define VIKIT_SOPHUS_COMPAT_H

#if __has_include(<sophus/se3.hpp>)
#include <sophus/se3.hpp>
#define VIKIT_SOPHUS_HAS_SE3_HPP 1
#elif __has_include(<sophus/se3.h>)
#include <sophus/se3.h>
#define VIKIT_SOPHUS_HAS_SE3_HPP 0
#else
#error "Cannot find Sophus SE3 header. Expected sophus/se3.hpp or sophus/se3.h."
#endif

namespace fast_livo
{
#if VIKIT_SOPHUS_HAS_SE3_HPP
using SE3d = Sophus::SE3<double>;
#else
using SE3d = Sophus::SE3;
#endif
}

using fast_livo::SE3d;

#if !VIKIT_SOPHUS_HAS_SE3_HPP
#define rotationMatrix rotation_matrix
#endif

#endif
