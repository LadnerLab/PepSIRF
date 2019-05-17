#ifndef OPTIONS_HH_INCLUDED
#define OPTIONS_HH_INCLUDED
#include <string>

class options
{
public:
    options();

    std::string input_r1_fname;
    std::string input_r2_fname;
    std::string library_fname;

    int num_threads;
};

#endif //OPTIONS_HH_INCLUDED
