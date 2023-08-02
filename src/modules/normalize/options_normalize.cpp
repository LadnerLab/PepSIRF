#include <sstream>
#include "options_normalize.h"


options_normalize::options_normalize() = default;

std::string options_normalize::get_arguments()
{
    std::ostringstream str_stream;
    str_stream << "--peptide_scores     " << peptide_scores_fname << "\n " <<
                  "--normalize_approach " << approach << "\n " <<
                  "--negative_control   " << neg_control << "\n " <<
                  "--negative_id        " << neg_id << "\n " <<
                  "--negative_names     " << neg_names << "\n " <<
                  "--precision   " << precision_digits << "\n " <<
                  "--output             " << output_fname << "\n" <<
                  "--logfile            " << logfile << "\n" <<
        "\n";

    return str_stream.str();
}

