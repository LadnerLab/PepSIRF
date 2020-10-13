#include <sstream>
#include <algorithm>

#include "options_subjoin.h"


namespace std
{
    // template specialization of std::hash
    // for evaluation_strategy::duplicate_resolution_strategy

    using strategy = evaluation_strategy::duplicate_resolution_strategy;

    template<>
    struct hash<strategy>
    {
        std::size_t operator()( const strategy strat ) const
        {
            return static_cast<std::size_t>( strat );
        }
    };

}; // namespace std



namespace evaluation_strategy
{
    duplicate_resolution_strategy from_string( const std::string& str )
    {
        return drs_string_map.at( str );
    }

    bool is_valid( const std::string& check )
    {
        return !( drs_string_map.find( check ) == drs_string_map.end() );
    }

    std::unordered_map<std::string,duplicate_resolution_strategy> drs_string_map
    {
        { "combine", duplicate_resolution_strategy::COMBINE },
        { "include", duplicate_resolution_strategy::INCLUDE },
        { "ignore",  duplicate_resolution_strategy::IGNORE }
    };

    std::unordered_map<duplicate_resolution_strategy,std::string> string_drs_map
    {
        { duplicate_resolution_strategy::COMBINE, "combine" },
        { duplicate_resolution_strategy::INCLUDE, "include" },
        { duplicate_resolution_strategy::IGNORE,  "ignore" }
    };


    std::string to_string( duplicate_resolution_strategy strategy )
    {
        return string_drs_map.find( strategy )->second;
    }

};

options_subjoin::options_subjoin() = default;

std::string options_subjoin::get_arguments()
{
    std::ostringstream str_stream;

    std::for_each( multi_matrix_name_pairs.begin(),
                   multi_matrix_name_pairs.end(),
                   [&]( const std::pair<std::string,std::string>& str )
                   {
                       str_stream << "--multi_file        "
                                  << str.first << "," << str.second
                                  << "\n ";
                   }
                 );

    std::for_each( input_matrix_name_pairs.begin(),
                   input_matrix_name_pairs.end(),
                   [&]( const std::pair<std::string,std::string>& str )
                   {
                       str_stream << "--filter_scores        "
                                  << str.first << "," << str.second
                                  << "\n ";
                   }
                 );
    str_stream << "--filter_peptide_names " << std::boolalpha << !use_sample_names << "\n ";
    str_stream << "--duplicate_evaluation " <<
        evaluation_strategy::string_drs_map[ duplicate_resolution_strategy ]
               << "\n ";
    str_stream << "--output               " << out_matrix_fname << "\n " <<
                  "\n";

    return str_stream.str();
}

