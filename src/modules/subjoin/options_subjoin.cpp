#include <sstream>

#include "options_subjoin.h"

std::string options_subjoin::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << "hello world";

    return str_stream.str();
}
