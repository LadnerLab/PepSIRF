#include "options_bin.h"
#include <sstream>

options_bin::options_bin() = default;

std::string options_bin::get_arguments()
{
    std::ostringstream str_stream;

    return str_stream.str();
};


