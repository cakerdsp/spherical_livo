#include "spherical.h"
#include <gtest/gtest.h>
#include <Eigen/Cholesky>
#include <Eigen/LU>

using namespace spherical;
namespace {
CameraModel pinhole(){CameraModel m;m.width=160;m.height=128;m.fx=m.fy=120;m.cx=79.5;m.cy=63.5;return m;}
Mat expRotation(const Vec& x){if(x.norm()<1e-15)return Mat::Identity();return Eigen::AngleAxisd(x.norm(),x.normalized()).toRotationMatrix();}
}

TEST(SphericalGeometry,PlaneWarpIsExactUnderTranslationAndRotation){
  const Vec n=Vec(0.1,-0.2,1).normalized();const double h=4;
  const Mat r=expRotation(Vec(0.03,0.1,-0.04));const Vec t(0.4,0.1,-0.2);
  for(const Vec& b:capTemplate(48,0.03)){
    Vec p;ASSERT_TRUE(intersectPlane(b,n,h,p));EXPECT_NEAR(n.dot(p),h,1e-12);
    const Vec exact=(r*p+t).normalized();const Vec homography=((r+t*n.transpose()/h)*b).normalized();
    EXPECT_LT((exact-homography).norm(),1e-12);
  }
}

TEST(SphericalGeometry,BearingDerivativeMatchesRightPerturbedState){
  const Mat rw=expRotation(Vec(0.13,-0.23,0.05)),rc=expRotation(Vec(-0.12,0.3,0.1));
  const Vec p(0.2,-0.1,1),tc(0.3,0.02,-0.05),pw(1,0.5,6);
  const auto j=bearingJacobian(rw,p,rc,tc,pw);const double eps=1e-6;
  for(int axis=0;axis<6;++axis){
    Vec d=Vec::Zero();d[axis%3]=eps;
    const Vec plus=axis<3?(rc*(rw*expRotation(d)).transpose()*(pw-p)+tc).normalized():
      (rc*rw.transpose()*(pw-p-d)+tc).normalized();
    const Vec minus=axis<3?(rc*(rw*expRotation(-d)).transpose()*(pw-p)+tc).normalized():
      (rc*rw.transpose()*(pw-p+d)+tc).normalized();
    EXPECT_LT(((plus-minus)/(2*eps)-j.col(axis)).norm(),1e-8);
  }
}

TEST(SphericalGeometry,TriangleGradientMatchesItsFrozenScaleInterpolant){
  auto model=pinhole();auto atlas=std::make_shared<RayAtlas>(model);cv::Mat raw(128,160,CV_8UC1);
  for(int y=0;y<raw.rows;++y)for(int x=0;x<raw.cols;++x)raw.at<uchar>(y,x)=40+(3*x+2*y)%150;
  Image image(atlas,raw);const Vec b=model.ray(83.17,59.39);
  const Basis e=tangentBasis(b);const double eps=1e-6;
  for(const SamplingFilter filter:{SamplingFilter{0,0},SamplingFilter{1,0.35}}) {
    Sample base;ASSERT_TRUE(image.sample(b,filter,base));
    for(int k=0;k<2;++k){Sample plus,minus;
      ASSERT_TRUE(image.sample((b+eps*e.col(k)).normalized(),filter,plus));
      ASSERT_TRUE(image.sample((b-eps*e.col(k)).normalized(),filter,minus));
      EXPECT_NEAR((plus.value-minus.value)/(2*eps),base.gradient.dot(e.col(k)),1e-5);
    }
    EXPECT_NEAR(base.gradient.dot(b),0,1e-10);
  }
}

TEST(SphericalGeometry,ClippingPropagatesThroughPrefilterSupport){
  auto atlas=std::make_shared<RayAtlas>(pinhole());cv::Mat raw(128,160,CV_8UC1,cv::Scalar(100));
  const Vec b=(0.4*atlas->ray(80,64)+0.3*atlas->ray(81,64)+0.3*atlas->ray(81,65)).normalized();
  // Outside the native triangle, but inside the next level's filter footprint.
  raw.at<uchar>(63,78)=255;Image image(atlas,raw);Sample sample;
  EXPECT_TRUE(image.sample(b,SamplingFilter{0,0},sample));
  EXPECT_FALSE(image.sample(b,SamplingFilter{1,0},sample));
}

TEST(SphericalGeometry,SharedTriangleEdgeIsContinuous){
  auto atlas=std::make_shared<RayAtlas>(pinhole());cv::Mat raw(128,160,CV_8UC1);
  for(int y=0;y<128;++y)for(int x=0;x<160;++x)raw.at<uchar>(y,x)=40+(3*x+2*y)%150;
  Image image(atlas,raw);
  const Vec a=atlas->ray(80,60),c=atlas->ray(81,61);
  const Vec edge=(0.37*a+0.63*c).normalized(),normal=a.cross(c).normalized();
  Sample minus,plus,center;
  ASSERT_TRUE(image.sample(edge,SamplingFilter{},center));
  ASSERT_TRUE(image.sample((edge-1e-8*normal).normalized(),SamplingFilter{},minus));
  ASSERT_TRUE(image.sample((edge+1e-8*normal).normalized(),SamplingFilter{},plus));
  EXPECT_NEAR(plus.value,minus.value,1e-5);
  const double expected=(0.37*raw.at<uchar>(60,80)+0.63*raw.at<uchar>(61,81))/255.0;
  EXPECT_NEAR(center.value,expected,1e-7);
}

TEST(SphericalGeometry,PrefilterPreservesConstantBrightnessAndSuppressesCheckerboard){
  auto model=pinhole();auto atlas=std::make_shared<RayAtlas>(model);
  cv::Mat constant(128,160,CV_8UC1,cv::Scalar(130)),checker(128,160,CV_8UC1);
  for(int y=0;y<128;++y)for(int x=0;x<160;++x)checker.at<uchar>(y,x)=(x+y)%2?80:180;
  Image flat(atlas,constant),alternating(atlas,checker);
  const Vec b=atlas->ray(80,60).normalized();Sample smooth,original;
  ASSERT_TRUE(flat.sample(b,SamplingFilter{2,0.4},smooth));
  EXPECT_NEAR(smooth.value,130.0/255,1e-6);EXPECT_LT(smooth.gradient.norm(),1e-4);
  ASSERT_TRUE(alternating.sample(b,SamplingFilter{},original));
  ASSERT_TRUE(alternating.sample(b,SamplingFilter{2,0},smooth));
  EXPECT_GT(std::abs(original.value-130.0/255),0.1);
  EXPECT_NEAR(smooth.value,130.0/255,0.005);
}

TEST(SphericalGeometry,LookupUsesFullCameraSphereNotPositiveZGate){
  auto m=pinhole();m.type=CameraModel::Type::KannalaBrandt;m.fx=m.fy=45;
  RayAtlas atlas(m);const Vec b=m.ray(155,64);ASSERT_LT(b.z(),0);
  Eigen::Vector2d uv;double pitch;ASSERT_TRUE(atlas.locate(b,uv,pitch));
  EXPECT_NEAR(uv.x(),155,0.01);EXPECT_NEAR(uv.y(),64,0.01);
}

TEST(SphericalGeometry,ExposureEliminationRetainsOnlyIdentifiableInformation){
  // A spatial ramp under dilation is indistinguishable from gain. With no
  // independent structure, eliminating gain must remove this visual constraint.
  Eigen::Vector3d j(1,2,3);const Eigen::Matrix2d h=(Eigen::Vector2d(1,1)*Eigen::Vector2d(1,1).transpose())*j.squaredNorm();
  EXPECT_NEAR(h(0,0)-h(0,1)*h(1,0)/h(1,1),0,1e-12);
}


TEST(SphericalGeometry,NativeProjectionSeedsStillRequireSphericalContainment){
  for(auto type:{CameraModel::Type::Pinhole,CameraModel::Type::KannalaBrandt,CameraModel::Type::Mei}) {
    auto model=pinhole();model.type=type;
    if(type==CameraModel::Type::KannalaBrandt) model.fx=model.fy=45;
    if(type==CameraModel::Type::Mei) model.xi=3.177956;
    model.distortion[0]=-0.015;model.distortion[1]=0.006;
    RayAtlas atlas(model);int checked=0;
    for(int y=30;y<100;y+=17) for(int x=25;x<140;x+=21) {
      const Vec b=model.ray(x+0.2,y+0.3);
      if(b.squaredNorm()<0.9) continue;
      Eigen::Vector2d seed,located;double pitch;
      ASSERT_TRUE(model.project(b,seed));
      EXPECT_NEAR(seed.x(),x+0.2,1e-6);EXPECT_NEAR(seed.y(),y+0.3,1e-6);
      ASSERT_TRUE(atlas.locate(b,located,pitch));
      EXPECT_LT((located-seed).norm(),0.1);EXPECT_GT(pitch,0);++checked;
    }
    EXPECT_GT(checked,0);
  }
}

TEST(SphericalGeometry,FinitePlaneIncludesAllPrefilterContributors){
  auto model=pinhole();auto atlas=std::make_shared<RayAtlas>(model);
  cv::Mat raw(128,160,CV_8UC1,cv::Scalar(100));Image image(atlas,raw);
  const Vec b=model.ray(80.3,64.2),normal=Vec::UnitZ();const double height=4;
  const Vec center=b*(height/normal.dot(b));
  Sample sample;
  PlaneSupport large{normal,center,height,3.0},small{normal,center,height,0.05};
  EXPECT_TRUE(image.sample(b,SamplingFilter{2,0},sample,true,&large));
  EXPECT_FALSE(image.sample(b,SamplingFilter{2,0},sample,true,&small));
  SamplingFilter filter;double cover;
  const Basis basis=tangentBasis(b);
  ASSERT_TRUE(image.prepare(b,0.0001*basis*basis.transpose(),filter,cover,&large));
  EXPECT_GT(cover,0);EXPECT_GE(filter.level,0);
  EXPECT_TRUE(image.sample(b,filter,sample,true,&large));
}

TEST(SphericalGeometry,TriangleSamplingSupportsDirectionsBehindTheCameraPlane){
  auto model=pinhole();model.type=CameraModel::Type::KannalaBrandt;model.fx=model.fy=45;
  auto atlas=std::make_shared<RayAtlas>(model);cv::Mat raw(128,160,CV_8UC1,cv::Scalar(100));
  Image image(atlas,raw);Sample sample;const Vec b=model.ray(155,64);
  ASSERT_LT(b.z(),0);ASSERT_TRUE(image.sample(b,SamplingFilter{},sample));
  EXPECT_NEAR(sample.value,100.0/255,1e-6);
}
