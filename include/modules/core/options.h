#ifndef OPTIONS_HH_INCLUDED
#define OPTIONS_HH_INCLUDED
#include <string>

/*! Data class to contain and handle 
 * arguments that will be passed in from the 
 * command-line.
*/
class options
{
public:
    options(): DEFAULT_NUM_THREADS( 2 ){} //!< Default constructor.
    int num_threads; //!< The number of threads to use for computation.
    const int DEFAULT_NUM_THREADS; //!< The default number of threads to use
};

#endif //OPTIONS_HH_INCLUDED
