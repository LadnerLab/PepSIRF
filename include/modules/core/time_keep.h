#ifndef TIME_KEEP_HH_INCLUDED
#define TIME_KEEP_HH_INCLUDED

namespace time_keep
{ 

/**
 * A struct to keep track of time that has elapsed.
 **/
struct timer
{
    double start;
    double end;
};

 double get_elapsed( const struct timer& time );
}

inline double time_keep::get_elapsed( const struct time_keep::timer& time )
{
    return time.end - time.start; 
}

#endif // TIME_KEEP_HH_INCLUDED
