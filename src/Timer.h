#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>

class Timer {
 public:
  void Reset();
  double Elapsed() const;

 private:
  using Clock = std::chrono::high_resolution_clock;
  using Second = std::chrono::duration<double, std::ratio<1>>;
  std::chrono::time_point<Clock> mStart{Clock::now()};
};

#endif  // TIMER_H_
