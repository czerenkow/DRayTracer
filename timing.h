#pragma once
#include <time.h>
#include <errno.h>
#include <chrono>

struct TimeMeasure {
    double full_duration = 0.0; // in seconds
    std::chrono::high_resolution_clock::time_point start_time;

    inline void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    /* Returns the number of seconds elapsed since the last start() call. */
    inline double stop() {
        const auto endv = std::chrono::high_resolution_clock::now();
        if (start_time.time_since_epoch() == std::chrono::high_resolution_clock::time_point::duration::zero()) {
            start_time = std::chrono::high_resolution_clock::time_point(); // zeroing
            return 0.0;
        }
        const auto duration = endv - start_time;
        const double result = duration.count() * 1e-9;
        full_duration += result;
        return result;
    }
};




/* msleep(): Sleep for the requested number of milliseconds. */
inline int msleep(long msec)
{
    struct timespec ts;
    int res;

    if (msec < 0)
    {
        errno = EINVAL;
        return -1;
    }

    ts.tv_sec = msec / 1000;
    ts.tv_nsec = (msec % 1000) * 1000000;

    do {
        res = nanosleep(&ts, &ts);
    } while (res && errno == EINTR);

    return res;
}
