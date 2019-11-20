#ifndef TIME_KEEP_HH_INCLUDED
#define TIME_KEEP_HH_INCLUDED
#include "omp_opt.h"

namespace time_keep
{ 

/**
 * A struct to keep track of wall clock time that has elapsed.
 **/
struct timer
{
    double begin_t;
    double end_t;

    void start() { begin_t = omp_get_wtime(); };
    void stop()   { end_t = omp_get_wtime(); };
    double get_elapsed() const { return end_t - begin_t; };

};

 double get_elapsed( const struct timer& time );
}

/**
 * Get the elapsed time (in seconds) of the timer. 
 * @param time The timer from which to get elapsed time.
 **/
inline double time_keep::get_elapsed( const struct time_keep::timer& time )
{
    return time.get_elapsed();
}

#endif // TIME_KEEP_HH_INCLUDED
