#include "spherical.h"
#include <opencv2/imgproc.hpp>
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

bool CameraModel::project(const Vec& p,Eigen::Vector2d& pixel) const {
  if(!p.allFinite() || p.squaredNorm()<1e-20) return false;
  double x,y;
  if(type==Type::KannalaBrandt) {
    const double r=std::hypot(p.x(),p.y());
    if(r<1e-12) { if(p.z()<=0) return false; pixel=Eigen::Vector2d(cx,cy); return true; }
    const double theta=std::atan2(r,p.z()),t2=theta*theta;
    double rd=theta,power=theta*t2;
    for(int k=0;k<4;++k) { rd+=distortion[k]*power;power*=t2; }
    pixel=Eigen::Vector2d(fx*rd*p.x()/r+cx,fy*rd*p.y()/r+cy);
  } else {
    const double den=type==Type::Mei ? p.z()+xi*p.norm() : p.z();
    if(den<=1e-12) return false;
    x=p.x()/den;y=p.y()/den;
    const double r2=x*x+y*y,radial=1+distortion[0]*r2+distortion[1]*r2*r2+distortion[4]*r2*r2*r2;
    pixel=Eigen::Vector2d(fx*(x*radial+2*distortion[2]*x*y+distortion[3]*(r2+2*x*x))+cx,
                         fy*(y*radial+distortion[2]*(r2+2*y*y)+2*distortion[3]*x*y)+cy);
  }
  return pixel.allFinite();
}

bool RayAtlas::triangleCell(const Vec& b,int x,int y,int level,Triangle& out) const {
  if(level<0 || level>=levels()) return false;
  const auto size=masses_[level].size();
  if(x<0 || y<0 || x>=size.width-1 || y>=size.height-1) return false;
  const int step=1<<level;
  for(int side=0;side<2;++side) {
    const std::array<cv::Point,3> vertices{{{x,y},{x+1,side==0?y:y+1},{side==0?x+1:x,y+1}}};
    const Vec r0=ray(vertices[0].x*step,vertices[0].y*step),
              r1=ray(vertices[1].x*step,vertices[1].y*step),
              r2=ray(vertices[2].x*step,vertices[2].y*step);
    if(r0.squaredNorm()<0.9 || r1.squaredNorm()<0.9 || r2.squaredNorm()<0.9) continue;
    Mat inverse;
    inverse.row(0)=r1.cross(r2).transpose();
    inverse.row(1)=r2.cross(r0).transpose();
    inverse.row(2)=r0.cross(r1).transpose();
    const double determinant=r0.dot(inverse.row(0));
    if(std::abs(determinant)<1e-14) continue;
    inverse/=determinant;
    const Vec a=inverse*b;const double sum=a.sum();
    if(a.minCoeff() < -1e-9 || !(sum>0)) continue;
    out.vertices=vertices;out.weights=a/sum;
    out.derivative=(inverse-out.weights*inverse.colwise().sum())/sum;
    out.pixel.setZero();
    for(int i=0;i<3;++i) out.pixel+=out.weights[i]*Eigen::Vector2d(vertices[i].x,vertices[i].y);
    return true;
  }
  return false;
}

bool RayAtlas::locateCell(const Vec& b,int x,int y,Eigen::Vector2d& pixel) const {
  Triangle found;
  if(!triangleCell(b,x,y,0,found)) return false;
  pixel=found.pixel;return true;
}

RayAtlas::RayAtlas(const CameraModel& model):model_(model),width_(model.width),height_(model.height) {
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
  // Binomial 5-tap pyramid centers coincide with raw pixel (2^level*x,2^level*y).
  // Area masses and support envelopes depend only on calibration, not frames.
  masses_.emplace_back(height_,width_,CV_32FC1,areas_.data());
  radii_.emplace_back(); // native pixels have zero prefilter support radius
  while(masses_.back().rows>=16 && masses_.back().cols>=16) {
    const int level=int(masses_.size()),step=1<<level;
    cv::Mat mass;cv::pyrDown(masses_.back(),mass);
    cv::Mat radius(mass.size(),CV_32FC1,cv::Scalar(2));
    const auto& previous=masses_.back();
    for(int y=0;y<mass.rows;++y) for(int x=0;x<mass.cols;++x) {
      if(2*x<2 || 2*y<2 || 2*x+2>=previous.cols || 2*y+2>=previous.rows) continue;
      const Vec center=ray(x*step,y*step);
      if(center.squaredNorm()<0.9) continue;
      double bound=0;
      for(int dy=-2;dy<=2;++dy) for(int dx=-2;dx<=2;++dx) {
        const cv::Point child(2*x+dx,2*y+dy);
        const Vec direction=ray(child.x*(step/2),child.y*(step/2));
        const double child_radius=level==1?0:radii_.back().at<float>(child);
        bound=std::max(bound,(direction-center).norm()+child_radius);
      }
      // Inflate for float LUT/envelope storage; the exact check remains fallback.
      radius.at<float>(y,x)=float(std::min(2.0,bound+1e-6));
    }
    masses_.push_back(mass);radii_.push_back(radius);
  }
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
  // Native projection is only a seed: accept only after the same spherical
  // LUT-triangle containment test. Irregular/ambiguous regions use the old KD path.
  Eigen::Vector2d seed_pixel;
  if(model_.project(b,seed_pixel) && seed_pixel.x()>-2 && seed_pixel.y()>-2 &&
     seed_pixel.x()<width_+1 && seed_pixel.y()<height_+1) {
    const int sx=std::lround(seed_pixel.x()),sy=std::lround(seed_pixel.y());
    int nx=-1,ny=-1;double distance=std::numeric_limits<double>::infinity();
    for(int y=std::max(1,sy-1);y<=std::min(height_-2,sy+1);++y)
      for(int x=std::max(1,sx-1);x<=std::min(width_-2,sx+1);++x) {
        if(area(x,y)<=0) continue;
        const double d=(ray(x,y)-b).squaredNorm();
        if(d<distance) { distance=d;nx=x;ny=y; }
      }
    if(nx>=0) {
      const int cell_x=std::floor(seed_pixel.x()),cell_y=std::floor(seed_pixel.y());
      if(locateCell(b,cell_x,cell_y,pixel)) { pitch=std::sqrt(double(area(nx,ny)));return true; }
      for(int y=cell_y-1;y<=cell_y+1;++y) for(int x=cell_x-1;x<=cell_x+1;++x)
        if(locateCell(b,x,y,pixel)) { pitch=std::sqrt(double(area(nx,ny)));return true; }
    }
  }
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
    for(int u=std::max(0,x-2);u<=std::min(width_-2,x+1);++u)
      if(locateCell(b,u,v,pixel)) { pitch=std::sqrt(std::max(1e-16,double(area(x,y))));return true; }
  return false;
}

bool RayAtlas::triangle(const Vec& b,int level,Triangle& out) const {
  if(level<0 || level>=levels() || !b.allFinite() || std::abs(b.squaredNorm()-1)>1e-5) return false;
  Eigen::Vector2d pixel;
  const double step=double(1<<level);
  const auto nearby=[&](const Eigen::Vector2d& raw) {
    const Eigen::Vector2d coarse=raw/step;
    const auto size=masses_[level].size();
    if(coarse.x() < -2 || coarse.y() < -2 || coarse.x()>size.width+1 || coarse.y()>size.height+1) return false;
    const int x=int(std::floor(coarse.x())),y=int(std::floor(coarse.y()));
    if(triangleCell(b,x,y,level,out)) return true;
    for(int dy=-1;dy<=1;++dy) for(int dx=-1;dx<=1;++dx)
      if((dx || dy) && triangleCell(b,x+dx,y+dy,level,out)) return true;
    return false;
  };
  if(model_.project(b,pixel) && nearby(pixel)) return true;
  // Preserve the LUT fallback for ambiguous/nonlinear native projection regions.
  double pitch;
  return locate(b,pixel,pitch) && nearby(pixel);
}

double RayAtlas::supportRadius(int level,const cv::Point& pixel) const {
  return level==0 ? 0.0 : double(radii_.at(level).at<float>(pixel));
}

cv::Rect RayAtlas::supportBox(int level,const cv::Point& pixel) const {
  const int step=1<<level,radius=2*(step-1);
  return cv::Rect(pixel.x*step-radius,pixel.y*step-radius,2*radius+1,2*radius+1);
}

bool PlaneSupport::contains(const Vec& ray) const {
  const double den=normal.dot(ray);
  if(std::abs(den)<1e-6*ray.norm() || !(height/den>0)) return false;
  return (ray*(height/den)-center).squaredNorm()<=radius*radius;
}

bool PlaneSupport::containsKernel(const Vec& bearing,double radius_bound) const {
  const double length=bearing.norm();
  if(!(length>0)) return false;
  const Vec b=bearing/length;
  const double r=radius_bound+std::abs(length-1.0);
  const double den=normal.dot(b),n=normal.norm(),lower=std::abs(den)-n*r;
  if(lower<=1e-6*(1+r) || !(height/den>0)) return false;
  const double bound=std::abs(height)*r*n/(lower*std::abs(den));
  return (b*(height/den)-center).norm()+bound+1e-10<=radius;
}

Image::Image(std::shared_ptr<const RayAtlas> atlas,const cv::Mat& gray):atlas_(std::move(atlas)) {
  if(gray.type()!=CV_8UC1 || gray.cols!=atlas_->width() || gray.rows!=atlas_->height())
    throw std::invalid_argument("image must be mono8 with calibrated dimensions");
  levels_.resize(atlas_->levels());
  gray.convertTo(levels_[0].values,CV_32F,1.0/255.0);
  cv::inRange(gray,cv::Scalar(2),cv::Scalar(253),levels_[0].valid);
  levels_[0].valid.setTo(0,atlas_->mass(0)<=0);
  levels_[0].values.setTo(0,levels_[0].valid==0);
  const cv::Mat kernel=cv::Mat::ones(5,5,CV_8U);
  for(int level=1;level<atlas_->levels();++level) {
    const auto& previous=levels_[level-1];auto& next=levels_[level];
    cv::Mat numerator=previous.values.mul(atlas_->mass(level-1)),filtered,valid;
    cv::pyrDown(numerator,filtered,atlas_->mass(level).size());
    cv::divide(filtered,atlas_->mass(level),next.values);
    // Every positive-weight contributor must be valid, including image edges.
    cv::erode(previous.valid,valid,kernel,cv::Point(-1,-1),1,cv::BORDER_CONSTANT,cv::Scalar(0));
    next.valid.create(next.values.size(),CV_8UC1);
    for(int y=0;y<next.valid.rows;++y) {
      auto* dst=next.valid.ptr<uchar>(y);const auto* src=valid.ptr<uchar>(2*y);
      for(int x=0;x<next.valid.cols;++x) dst[x]=src[2*x];
    }
    next.valid.setTo(0,atlas_->mass(level)<=0);
    next.values.setTo(0,next.valid==0);
  }
}

bool Image::locate(const Vec& b,Eigen::Vector2d& pixel,double& pitch) const {
  return atlas_->locate(b,pixel,pitch);
}

bool Image::inspect(const Vec& b,int level,Triangle& triangle,double& cover,const PlaneSupport* support) const {
  if(!atlas_->triangle(b,level,triangle)) return false;
  const int step=1<<level;
  cover=0;
  // Validate the entire chosen triangle, including vertices with zero weight:
  // an invalid vertex must not become a hidden derivative contribution.
  for(const auto& vertex:triangle.vertices) {
    if(!levels_[level].valid.at<uchar>(vertex)) return false;
    const Vec ray=atlas_->ray(vertex.x*step,vertex.y*step);
    const double radius=atlas_->supportRadius(level,vertex);
    cover=std::max(cover,(ray.normalized()-b).norm()+radius+1e-6);
    if(!support) continue;
    if(level==0) { if(!support->contains(ray)) return false;continue; }
    if(support->containsKernel(ray,radius)) continue;
    // Rare finite-plane boundary case: retain the exact contributor check.
    // This is only a geometry fallback, never a weighted intensity search.
    const cv::Rect box=atlas_->supportBox(level,vertex);
    if(box.x<0 || box.y<0 || box.x+box.width>atlas_->width() || box.y+box.height>atlas_->height()) return false;
    for(int y=box.y;y<box.y+box.height;++y) for(int x=box.x;x<box.x+box.width;++x)
      if(!support->contains(atlas_->ray(x,y))) return false;
  }
  return true;
}

bool Image::prepare(const Vec& b,const Mat& covariance,SamplingFilter& filter,double& cover,
                    const PlaneSupport* support) const {
  if(!covariance.allFinite()) return false;
  Triangle native;
  if(!atlas_->triangle(b,0,native)) return false;
  Eigen::Matrix<double,2,3> pixel_derivative=Eigen::Matrix<double,2,3>::Zero();
  for(int i=0;i<3;++i) {
    pixel_derivative.row(0)+=native.vertices[i].x*native.derivative.row(i);
    pixel_derivative.row(1)+=native.vertices[i].y*native.derivative.row(i);
  }
  const Eigen::Matrix2d footprint=pixel_derivative*covariance*pixel_derivative.transpose();
  const double variance=0.5*(footprint.trace()+std::hypot(footprint(0,0)-footprint(1,1),2*footprint(0,1)));
  const double minimum=0.5*(footprint.trace()-std::hypot(footprint(0,0)-footprint(1,1),2*footprint(0,1)));
  if(!std::isfinite(variance) || minimum < -1e-10*std::max(1.0,variance) || variance<0) return false;
  // A recursive [1 4 6 4 1]/16 pyramid has sensor-pixel variance (4^L-1)/3.
  const int last=int(levels_.size())-1;
  if(variance>(std::ldexp(1.0,2*last)-1)/3.0) return false;
  filter.level=std::min(last,int(std::floor(0.5*std::log2(1+3*variance))));
  const double lower=(std::ldexp(1.0,2*filter.level)-1)/3.0;
  filter.blend=filter.level==last ? 0 : std::max(0.0,std::min(1.0,(variance-lower)/std::ldexp(1.0,2*filter.level)));
  if(filter.blend==1) { ++filter.level;filter.blend=0; }
  Triangle triangle;double radius;
  if(!inspect(b,filter.level,triangle,cover,support)) return false;
  if(filter.blend>0) {
    if(!inspect(b,filter.level+1,triangle,radius,support)) return false;
    cover=std::max(cover,radius);
  }
  return true;
}

bool Image::sample(const Vec& b,const SamplingFilter& filter,Sample& out,bool gradient,
                   const PlaneSupport* support) const {
  out.value=0;out.gradient.setZero();
  if(filter.level<0 || filter.level>=int(levels_.size()) ||
     !(filter.blend>=0 && filter.blend<=1) ||
     (filter.blend>0 && filter.level+1>=int(levels_.size()))) return false;
  for(int offset=0;offset<(filter.blend>0?2:1);++offset) {
    const int level=filter.level+offset;
    const double blend=offset==0?1-filter.blend:filter.blend;
    Triangle triangle;double cover;
    if(!inspect(b,level,triangle,cover,support)) return false;
    Vec values;
    for(int i=0;i<3;++i) values[i]=levels_[level].values.at<float>(triangle.vertices[i]);
    out.value+=blend*values.dot(triangle.weights);
    if(gradient) out.gradient+=blend*triangle.derivative.transpose()*values;
  }
  if(gradient) out.gradient-=b*b.dot(out.gradient);
  return std::isfinite(out.value) && out.gradient.allFinite();
}
} // namespace spherical
