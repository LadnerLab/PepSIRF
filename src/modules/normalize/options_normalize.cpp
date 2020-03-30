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

    std::string col_sum_str = bool_to_str( !( col_sum_norm && size_factors_norm ) );
    std::string size_factors_str = bool_to_str( size_factors_norm );

    str_stream << "--peptide_scores     " << peptide_scores_fname << "\n " <<
                  "--col_sum            "  <<  col_sum_str << "\n " << 
                  "--size_factors       "  << size_factors_str << "\n " << 
                  "--precision_digits   "  << precision_digits << "\n " << 
                  "--output             " << output_fname << "\n" <<
        "\n";

    return str_stream.str();
}
