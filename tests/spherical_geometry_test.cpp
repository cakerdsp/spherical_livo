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

TEST(SphericalGeometry,KernelGradientIncludesNormalizedWeights){
  auto atlas=std::make_shared<RayAtlas>(pinhole());cv::Mat raw(128,160,CV_8UC1);
  for(int y=0;y<raw.rows;++y)for(int x=0;x<raw.cols;++x)raw.at<uchar>(y,x)=40+(3*x+2*y)%150;
  Image image(atlas,raw);const Vec b=Vec(0.04,-0.03,1).normalized();
  const Mat a=Mat::Identity()/(0.04*0.04);Sample base;ASSERT_TRUE(image.sample(b,a,base));
  const Basis e=tangentBasis(b);const double eps=1e-6;
  for(int k=0;k<2;++k){Sample plus,minus;
    ASSERT_TRUE(image.sample((b+eps*e.col(k)).normalized(),a,plus));
    ASSERT_TRUE(image.sample((b-eps*e.col(k)).normalized(),a,minus));
    EXPECT_NEAR((plus.value-minus.value)/(2*eps),base.gradient.dot(e.col(k)),1e-4);
  }
  EXPECT_NEAR(base.gradient.dot(b),0,1e-10);
  EXPECT_GT(weightOverlap(base,base),0);
  EXPECT_LE(weightOverlap(base,base),1);
}

TEST(SphericalGeometry,ClippedNativePixelInvalidatesItsWholeSupport){
  auto atlas=std::make_shared<RayAtlas>(pinhole());cv::Mat raw(128,160,CV_8UC1,cv::Scalar(100));
  raw.at<uchar>(64,80)=255;Image image(atlas,raw);Sample sample;
  EXPECT_FALSE(image.sample(atlas->ray(80,64).normalized(),Mat::Identity()/0.001,sample));
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

TEST(SphericalGeometry,FastKernelChecksTheSameFinitePlanePixels){
  auto atlas=std::make_shared<RayAtlas>(pinhole());cv::Mat raw(128,160,CV_8UC1);
  for(int y=0;y<128;++y) for(int x=0;x<160;++x) raw.at<uchar>(y,x)=40+(3*x+2*y)%150;
  Image image(atlas,raw);const Vec b=Vec(0.03,-0.02,1).normalized();
  const Mat a=precision(b,(Eigen::Vector2d(0.04*0.04,0.025*0.025)).asDiagonal());
  const auto kernel=SamplingKernel::fromMatrix(a);
  Sample full;ASSERT_TRUE(image.sample(b,a,full));
  const Vec normal=Vec(0.2,-0.1,1).normalized();const double height=4;
  const Vec center=b*(height/normal.dot(b));
  for(double radius:{0.02,0.12,2.0}) {
    bool expected=true;
    for(const auto& weight:full.weights) {
      const Vec ray=atlas->ray(weight.first%160,weight.first/160);
      const Vec surface=ray*(height/normal.dot(ray));
      if((surface-center).norm()>radius) expected=false;
    }
    PlaneSupport support{normal,center,height,radius};Sample fast;
    EXPECT_EQ(image.sample(b,kernel,fast,true,&support),expected);
    if(expected) { EXPECT_NEAR(full.value,fast.value,1e-12);EXPECT_LT((full.gradient-fast.gradient).norm(),1e-10); }
    EXPECT_TRUE(fast.weights.empty());
  }
}
