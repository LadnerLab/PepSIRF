#ifndef OMP_OPT_HH_INCLUDED
#define OMP_OPT_HH_INCLUDED

#ifdef ENABLE_OPENMP
#include <omp.h>

#else
#include <time.h>
#include <sys/time.h>

#define omp_get_num_threads() 0
#define omp_set_num_threads

#define omp_get_wtime() [] () -> double { \
    struct timeval time; \
    if( gettimeofday( &time, NULL ) ) { \
    return 0; \
    } \
    return (double) time.tv_sec + (double) time.tv_usec * 0.000001; \
    }()

#endif // ifdef ENABLE_OPENMP

#endif // OMP_OPT_HH_INCLUDED
