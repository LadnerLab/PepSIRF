#include "options_info.h"
#include <sstream>

options_info::options_info() = default;
options_info::~options_info() = default;

std::string options_info::get_arguments()
{
    std::ostringstream stream;

    stream << "Hello world!\n";

    return stream.str();
}
