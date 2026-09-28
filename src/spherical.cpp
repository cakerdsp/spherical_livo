#include "spherical.h"
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace spherical {
Mat skew(const Vec& x) {
  Mat m;
  m << 0, -x.z(), x.y(), x.z(), 0, -x.x(), -x.y(), x.x(), 0;
  return m;
}

Basis tangentBasis(const Vec& b) {
  Basis e;
  const Vec axis = std::abs(b.z()) < 0.8 ? Vec::UnitZ() : Vec::UnitX();
  e.col(0) = axis.cross(b).normalized();
  e.col(1) = b.cross(e.col(0));
  return e;
}

Mat rightJacobianInverse(const Vec& x) {
  const double angle = x.norm();
  const Mat a = skew(x);
  const double c = angle < 1e-4 ? 1.0/12.0 + angle*angle/720.0 :
      (1.0 - 0.5*angle/std::tan(0.5*angle))/(angle*angle);
  return Mat::Identity() + 0.5*a + c*a*a;
}

std::vector<Vec> capTemplate(int count, double radius) {
  if (count < 8 || count > 256 || !(radius > 0 && radius < 0.2))
    throw std::invalid_argument("patch samples must be in [8,256], radius in (0,0.2) rad");
  std::vector<Vec> out;
  out.reserve(count);
  out.push_back(Vec::UnitZ());
  const double golden = M_PI*(3.0-std::sqrt(5.0));
  for (int k=1; k<count; ++k) {
    const double z = 1.0-(1.0-std::cos(radius))*(k-0.5)/(count-1);
    const double r = std::sqrt(std::max(0.0, 1-z*z));
    out.emplace_back(r*std::cos(golden*k), r*std::sin(golden*k), z);
  }
  return out;
}

bool intersectPlane(const Vec& b, const Vec& n, double h, Vec& p) {
  const double den = n.dot(b);
  if (std::abs(den) < 1e-6) return false;
  const double d = h/den;
  if (!(d > 0) || !std::isfinite(d)) return false;
  p = d*b;
  return true;
}

Eigen::Matrix<double,3,6> bearingJacobian(const Mat& Rwi, const Vec& pwi,
                                        const Mat& Rci, const Vec& tci, const Vec& pw) {
  const Vec pi = Rwi.transpose()*(pw-pwi);
  const Vec pc = Rci*pi+tci;
  const Vec b = pc.normalized();
  Eigen::Matrix<double,3,6> j;
  j.leftCols<3>() = Rci*skew(pi);
  j.rightCols<3>() = -Rci*Rwi.transpose();
  return (Mat::Identity()-b*b.transpose())/pc.norm()*j;
}

void CameraModel::validate() const {
  if (width < 16 || height < 16 || !std::isfinite(fx) || !std::isfinite(fy) ||
      fx <= 0 || fy <= 0 || !std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(xi))
    throw std::invalid_argument("invalid camera dimensions/intrinsics");
  for (double x: distortion) if (!std::isfinite(x))
    throw std::invalid_argument("non-finite camera distortion");
}

Vec CameraModel::ray(double u, double v) const {
  const double xd = (u-cx)/fx, yd = (v-cy)/fy;
  double x=xd, y=yd;
  if (type == Type::KannalaBrandt) {
    const double rd = std::hypot(xd, yd);
    if (rd < 1e-12) return Vec::UnitZ();
    double theta=rd;
    for (int i=0; i<20; ++i) {
      const double t2=theta*theta;
      double value=theta, derivative=1, power=theta*t2, dp=t2;
      for (int k=0;k<4;++k) {
        value += distortion[k]*power;
        derivative += (2*k+3)*distortion[k]*dp;
        power*=t2; dp*=t2;
      }
      if (!(derivative>1e-10)) return Vec::Zero();
      const double step=(value-rd)/derivative;
      theta-=step;
      if (std::abs(step)<1e-12) break;
    }
    if (!(theta>=0 && theta<M_PI-1e-4)) return Vec::Zero();
    double check=theta, power=theta*theta*theta;
    for(int k=0;k<4;++k){ check+=distortion[k]*power; power*=theta*theta; }
    if(std::abs(check-rd)>1e-8) return Vec::Zero();
    return Vec(std::sin(theta)*xd/rd,std::sin(theta)*yd/rd,std::cos(theta));
  }
  for (int i=0;i<25;++i) {
    const double r2=x*x+y*y;
    const double radial=1+distortion[0]*r2+distortion[1]*r2*r2+distortion[4]*r2*r2*r2;
    if (radial<=0 || !std::isfinite(radial)) return Vec::Zero();
    const double dx=2*distortion[2]*x*y+distortion[3]*(r2+2*x*x);
    const double dy=distortion[2]*(r2+2*y*y)+2*distortion[3]*x*y;
    const double nx=(xd-dx)/radial, ny=(yd-dy)/radial;
    const double error=std::hypot(nx-x,ny-y);
    x=nx; y=ny;
    if(error<1e-12) break;
    if(i==24) return Vec::Zero();
  }
  if(type == Type::Mei) {
    const double r2=x*x+y*y, rad=1+(1-xi*xi)*r2;
    if(rad<0) return Vec::Zero();
    const double factor=(xi+std::sqrt(rad))/(1+r2);
    return Vec(factor*x,factor*y,factor-xi).normalized();
  }
  return Vec(x,y,1).normalized();
}

RayAtlas::RayAtlas(const CameraModel& model):width_(model.width),height_(model.height) {
  model.validate();
  rays_.resize(width_*height_);
  areas_.resize(width_*height_,0);
  for(int y=0;y<height_;++y) for(int x=0;x<width_;++x)
    rays_[y*width_+x]=model.ray(x,y).cast<float>();
  for(int y=1;y<height_-1;++y) for(int x=1;x<width_-1;++x) {
    const Vec b=ray(x,y), dx=0.5*(ray(x+1,y)-ray(x-1,y)),dy=0.5*(ray(x,y+1)-ray(x,y-1));
    if(b.squaredNorm()<0.9 || ray(x+1,y).squaredNorm()<0.9 || ray(x-1,y).squaredNorm()<0.9 ||
       ray(x,y+1).squaredNorm()<0.9 || ray(x,y-1).squaredNorm()<0.9) continue;
    areas_[y*width_+x]=std::abs(b.dot(dx.cross(dy)));
    if(x%8==0 && y%8==0) seeds_.push_back({b.cast<float>(),y*width_+x});
  }
  if(seeds_.empty()) throw std::invalid_argument("camera has no valid ray support");
  std::vector<int> ids(seeds_.size()); std::iota(ids.begin(),ids.end(),0);
  nodes_.reserve(ids.size()); root_=build(ids,0,ids.size(),0);
}

int RayAtlas::build(std::vector<int>& ids,int begin,int end,int depth) {
  if(begin>=end) return -1;
  const int mid=(begin+end)/2,axis=depth%3;
  std::nth_element(ids.begin()+begin,ids.begin()+mid,ids.begin()+end,
    [&](int a,int b){return seeds_[a].b[axis]<seeds_[b].b[axis];});
  const int node=nodes_.size(); nodes_.push_back({ids[mid],-1,-1,axis});
  const int left=build(ids,begin,mid,depth+1),right=build(ids,mid+1,end,depth+1);
  nodes_[node].left=left;nodes_[node].right=right;return node;
}

void RayAtlas::nearest(int node,const Vec& b,int& best,double& distance) const {
  if(node<0) return;
  const Node& n=nodes_[node];const Vec d=b-seeds_[n.seed].b.cast<double>();
  const double sq=d.squaredNorm(); if(sq<distance){distance=sq;best=n.seed;}
  const bool low=d[n.axis]<0;
  nearest(low?n.left:n.right,b,best,distance);
  if(d[n.axis]*d[n.axis]<distance) nearest(low?n.right:n.left,b,best,distance);
}

Vec RayAtlas::ray(int x,int y) const {return rays_[y*width_+x].cast<double>();}
float RayAtlas::area(int x,int y) const {return areas_[y*width_+x];}

bool RayAtlas::locate(const Vec& b,Eigen::Vector2d& pixel,double& pitch) const {
  if(!b.allFinite() || std::abs(b.squaredNorm()-1)>1e-5) return false;
  int seed=-1;double best=std::numeric_limits<double>::infinity();nearest(root_,b,seed,best);
  int x=seeds_[seed].pixel%width_,y=seeds_[seed].pixel/width_;
  for(int step: {8,4,2,1}) {
    for(int pass=0;pass<4;++pass) {
      int nx=x,ny=y;double dist=(ray(x,y)-b).squaredNorm();
      for(int dy=-step;dy<=step;dy+=step) for(int dx=-step;dx<=step;dx+=step) {
        const int xx=x+dx,yy=y+dy;
        if(xx<1||yy<1||xx>=width_-1||yy>=height_-1||area(xx,yy)<=0) continue;
        const double d=(ray(xx,yy)-b).squaredNorm();if(d<dist){dist=d;nx=xx;ny=yy;}
      }
      if(nx==x&&ny==y) break;x=nx;y=ny;
    }
  }
  for(int v=std::max(0,y-2);v<=std::min(height_-2,y+1);++v)
    for(int u=std::max(0,x-2);u<=std::min(width_-2,x+1);++u) for(int tri=0;tri<2;++tri) {
      const std::array<Eigen::Vector2i,3> p = tri==0 ?
        std::array<Eigen::Vector2i,3>{Eigen::Vector2i(u,v),{u+1,v},{u+1,v+1}} :
        std::array<Eigen::Vector2i,3>{Eigen::Vector2i(u,v),{u+1,v+1},{u,v+1}};
      Mat basis;for(int j=0;j<3;++j) basis.col(j)=ray(p[j].x(),p[j].y());
      if(std::abs(basis.determinant())<1e-14) continue;
      Vec a=basis.fullPivLu().solve(b);
      if(a.minCoeff() < -1e-7 || a.sum()<=0) continue;a/=a.sum();pixel.setZero();
      for(int j=0;j<3;++j) pixel+=a[j]*p[j].cast<double>();
      pitch=std::sqrt(std::max(1e-16,double(area(x,y))));return true;
    }
  return false;
}

Image::Image(std::shared_ptr<const RayAtlas> atlas,const cv::Mat& gray):atlas_(std::move(atlas)) {
  if(gray.type()!=CV_8UC1||gray.cols!=atlas_->width()||gray.rows!=atlas_->height())
    throw std::invalid_argument("image must be mono8 with calibrated dimensions");
  gray_=gray.clone();bounds_=cv::Rect(0,0,gray.cols,gray.rows);
  std::vector<double> highpass;
  for(int y=1;y<gray.rows-1;y+=4) for(int x=1;x<gray.cols-1;x+=4)
    highpass.push_back(std::abs(double(gray.at<uchar>(y,x))-gray.at<uchar>(y+1,x)-
      gray.at<uchar>(y,x+1)+gray.at<uchar>(y+1,x+1))/510.0);
  if(!highpass.empty()){
    auto m=highpass.begin()+highpass.size()/2;std::nth_element(highpass.begin(),m,highpass.end());
    noise_=std::max(noise_,*m/0.6744897501960817);
  }
}

std::shared_ptr<Image> Image::crop(const cv::Rect& roi) const {
  const cv::Rect box=roi&bounds_;if(box.empty()) throw std::invalid_argument("empty reference crop");
  auto out=std::shared_ptr<Image>(new Image);out->atlas_=atlas_;out->bounds_=box;out->noise_=noise_;
  out->gray_=gray_(cv::Rect(box.x-bounds_.x,box.y-bounds_.y,box.width,box.height)).clone();return out;
}

bool Image::locate(const Vec& b,Eigen::Vector2d& pixel,double& pitch) const {
  if(!atlas_->locate(b,pixel,pitch)) return false;
  return pixel.x()>=bounds_.x && pixel.y()>=bounds_.y &&
    pixel.x()<bounds_.x+bounds_.width-1 && pixel.y()<bounds_.y+bounds_.height-1;
}

bool Image::sample(const Vec& b,const Mat& a,Sample& out) const {
  out=Sample{};Eigen::Vector2d uv;double pitch;
  if(!locate(b,uv,pitch)) return false;
  Eigen::SelfAdjointEigenSolver<Mat> eig(a);
  if(eig.info()!=Eigen::Success||eig.eigenvalues().minCoeff()<=0) return false;
  const double radius=1/std::sqrt(eig.eigenvalues().minCoeff());
  if(radius>0.3) return false;
  const int cx=std::lround(uv.x()),cy=std::lround(uv.y());
  int extent=std::max(2,int(std::ceil(1.5*radius/pitch)));
  if(extent>128) return false;
  for(int attempt=0;attempt<4;++attempt) {
    out=Sample{};double mass=0;Vec dm=Vec::Zero(),dv=Vec::Zero();bool expand=false;
    for(int y=cy-extent;y<=cy+extent;++y) for(int x=cx-extent;x<=cx+extent;++x) {
      if(x<0||y<0||x>=atlas_->width()||y>=atlas_->height()) continue;
      const Vec d=b-atlas_->ray(x,y);const double q=d.dot(a*d);
      if(q>=1||q<0) continue;
      if(x<=bounds_.x||y<=bounds_.y||x>=bounds_.x+bounds_.width-1||y>=bounds_.y+bounds_.height-1)
        return false;
      if(x==cx-extent||x==cx+extent||y==cy-extent||y==cy+extent) expand=true;
      const double area=atlas_->area(x,y);if(area<=0) return false;
      const int value=gray_.at<uchar>(y-bounds_.y,x-bounds_.x);
      if(value<=1||value>=254) return false;
      const double r=std::sqrt(q),t=1-r,w=area*t*t*t*t*(4*r+1);
      const Vec dw=-20*area*t*t*t*a*d;
      mass+=w;out.value+=w*value/255.0;dm+=dw;dv+=dw*value/255.0;
      out.weights.emplace_back(y*atlas_->width()+x,w);
    }
    if(expand){extent*=2;if(extent>128) return false;continue;}
    if(mass<=1e-20||out.weights.size()<3) return false;
    out.value/=mass;out.gradient=(Mat::Identity()-b*b.transpose())*(dv-out.value*dm)/mass;
    for(auto& w:out.weights)w.second/=mass;
    return out.gradient.allFinite();
  }
  return false;
}

double weightOverlap(const Sample& a,const Sample& b) {
  size_t i=0,j=0;double value=0;
  while(i<a.weights.size()&&j<b.weights.size()){
    if(a.weights[i].first<b.weights[j].first)++i;
    else if(a.weights[i].first>b.weights[j].first)++j;
    else {value+=a.weights[i++].second*b.weights[j++].second;}
  }
  return value;
}

Mat precision(const Vec& b,const Eigen::Matrix2d& footprint) {
  const Basis e=tangentBasis(b);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(footprint);
  return e*footprint.inverse()*e.transpose()+b*b.transpose()/eig.eigenvalues().maxCoeff();
}
} // namespace spherical
