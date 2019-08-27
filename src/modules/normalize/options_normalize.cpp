#include <sstream>
#include "options_normalize.h"


options_normalize::options_normalize() = default;

std::string options_normalize::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << "\n";

    return str_stream.str();
}
