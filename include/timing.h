#pragma once

#include <ros/ros.h>
#include <array>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

namespace spherical {
using SteadyClock = std::chrono::steady_clock;
inline double elapsedMs(SteadyClock::time_point start) {
  return std::chrono::duration<double, std::milli>(SteadyClock::now()-start).count();
}

// Aggregate every completed event, print at most once per wall second.
// Sensor /clock, bag pauses and timestamp epochs do not affect the timer.
template<size_t N> class TimingTable {
 public:
  explicit TimingTable(std::array<const char*,N> labels):labels_(labels) {}
  void record(const std::string& title,double stamp,const std::array<double,N>& values) {
    ++count_;
    for(size_t i=0;i<N;++i)sum_[i]+=values[i];
    const auto now=SteadyClock::now();
    if(count_>1&&std::chrono::duration<double>(now-last_print_).count()<1) return;
    last_print_=now;
    std::ostringstream out;
    out<<"\n[ "<<title<<" timing ] events="<<count_<<" stamp="<<std::fixed<<std::setprecision(6)<<stamp;
    out<<"\n"<<std::left<<std::setw(25)<<"Stage"<<std::right<<std::setw(13)<<"Current ms"<<std::setw(13)<<"Mean ms";
    out<<std::setprecision(3);
    for(size_t i=0;i<N;++i)
      out<<"\n"<<std::left<<std::setw(25)<<labels_[i]<<std::right<<std::setw(13)<<values[i]<<std::setw(13)<<sum_[i]/count_;
    ROS_INFO_STREAM(out.str());
  }
 private:
  std::array<const char*,N> labels_;
  std::array<double,N> sum_{};
  size_t count_=0;
  SteadyClock::time_point last_print_{};
};
} // namespace spherical
