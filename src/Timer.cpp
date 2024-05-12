#include "Timer.h"

void Timer::Reset() { mStart = Clock::now(); }

double Timer::Elapsed() const {
  return std::chrono::duration_cast<Second>(Clock::now() - mStart).count();
}
