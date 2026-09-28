#pragma once

#include <array>
#include <chrono>
#include <cstdio>
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
  explicit TimingTable(std::array<const char*,N> labels,size_t total_index=N-1,
                       const char* color="\033[1;36m")
    :labels_(labels),total_index_(total_index),color_(color) {}
  void record(const std::string& title,double stamp,const std::array<double,N>& values,
              const std::string& details="") {
    ++count_;
    for(size_t i=0;i<N;++i)sum_[i]+=values[i];
    const auto now=SteadyClock::now();
    if(count_>1&&std::chrono::duration<double>(now-last_print_).count()<1) return;
    last_print_=now;
    std::ostringstream out;
    constexpr const char* blue="\033[1;34m";
    constexpr const char* reset="\033[0m";
    const auto border=[&] {out<<blue<<"+-------------------------------------------------------------+"<<reset<<'\n';};
    const auto row=[&](const char* color,const std::string& label,const auto& value) {
      out<<color<<"| "<<std::left<<std::setw(29)<<label<<" | "<<std::setw(27)<<value<<" |"<<reset<<'\n';
    };
    out<<'\n';
    if(!details.empty())out<<color_<<"[ "<<title<<" ] "<<details<<reset<<'\n';
    border();
    const std::string heading=title.substr(0,61);
    const size_t padding=(61-heading.size())/2;
    out<<blue<<'|'<<std::string(padding,' ')<<heading<<std::string(61-heading.size()-padding,' ')<<'|'<<reset<<'\n';
    border();
    out<<std::fixed<<std::setprecision(6);
    row(blue,"Sensor Timestamp",stamp);
    row(blue,"Completed Events",count_);
    border();
    row(blue,"Algorithm Stage","Time (secs)");
    border();
    for(size_t i=0;i<N;++i)if(i!=total_index_)row(color_,labels_[i],values[i]*0.001);
    border();
    row(color_,"Current Total Time",values[total_index_]*0.001);
    row(color_,"Average Total Time",sum_[total_index_]/count_*0.001);
    border();
    // Match FAST-LIVO2/cake_slam's ANSI stdout table without a ROS log prefix.
    const std::string table=out.str();
    std::fwrite(table.data(),1,table.size(),stdout);
    std::fflush(stdout);
  }
 private:
  std::array<const char*,N> labels_;
  size_t total_index_;
  const char* color_;
  std::array<double,N> sum_{};
  size_t count_=0;
  SteadyClock::time_point last_print_{};
};
} // namespace spherical
