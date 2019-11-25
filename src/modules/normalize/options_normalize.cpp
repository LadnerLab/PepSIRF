#include <sstream>
#include "options_normalize.h"


options_normalize::options_normalize() = default;

std::string options_normalize::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << "--peptide_scores " << peptide_scores_fname << "\n " <<
                  "--output         " << output_fname << "\n" <<
                  " --threads        " << num_threads << "\n" <<
        "\n";

    return str_stream.str();
}
