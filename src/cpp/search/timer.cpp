#include <cstddef>
#include <sys/time.h>

#ifndef TIMER_HPP
#define TIMER_HPP
using namespace std;

class Timer {
public:
  double prevtime, curtime;
  double update_time(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    prevtime = curtime;
    curtime = tv.tv_sec + 1e-6 * tv.tv_usec;
    return curtime - prevtime;
  }
};
#endif
