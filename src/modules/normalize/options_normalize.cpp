#include <sstream>
#include "options_normalize.h"


options_normalize::options_normalize() = default;

std::string options_normalize::get_arguments()
{

    auto bool_to_str = []( bool val ) -> std::string
        {
            if( val )
                {
                    return "true";
                }
            return "false";
        };

    std::ostringstream str_stream;
    str_stream << "--peptide_scores     " << peptide_scores_fname << "\n " <<
                  "--normalize_apporach " << approach << "\n " <<
                  "--precision_digits   " << precision_digits << "\n " <<
                  "--output             " << output_fname << "\n" <<
        "\n";

    return str_stream.str();
}
