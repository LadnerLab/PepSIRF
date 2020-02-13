#include <unordered_set>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include "module_subjoin.h"
#include "matrix.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include "time_keep.h" 

module_subjoin::module_subjoin() = default;

module_subjoin::name_replacement_list
module_subjoin::parse_namelist( std::vector<std::string>& dest,
                                std::istream& file
                              )
{
    std::string line;
    std::vector<std::string> split_line;
    name_replacement_list ret_val;

    while( std::getline( file, line ) )
        {
            boost::trim( line );
            if( line.size() )
                {
                    boost::split( split_line, line, boost::is_any_of( "\t" ) );

                    // replace the old name with the new
                    if( split_line.size() == 2 && split_line[ 1 ].length() > 0 )
                        {
                            dest.emplace_back( split_line[ 1 ] );
                            ret_val.emplace_back( std::make_pair( split_line[ 0 ],
                                                                  split_line[ 1 ]
                                                                )
                                                );
                        }
                    else
                        {
                            dest.emplace_back( line );
                        }
                }
        }
    return ret_val;
}

void module_subjoin::run( options *opts )
{
    options_subjoin *s_opts = (options_subjoin*) opts;
    std::vector<peptide_score_data_sample_major>
        parsed_score_data;
    parsed_score_data.resize( s_opts->matrix_name_pairs.size() );
    std::uint32_t idx = 0;

    bool use_peptide_names = !s_opts->use_sample_names;

    time_keep::timer time;

    time.start();

    // for all of the ( score_matrix, names_to_filter ) pairs:
    #pragma omp parallel for num_threads( 2 ) private( idx ) schedule( dynamic ) \
            shared( s_opts, parsed_score_data ) 
    for( idx = 0; idx < s_opts->matrix_name_pairs.size(); ++idx )
        {
            auto &score_name_pair = s_opts->matrix_name_pairs[ idx ];
            peptide_score_data_sample_major&
                my_data = parsed_score_data[ idx ];
                
            std::string& matrix_name_list  = score_name_pair.first;
            std::ifstream names_list( score_name_pair.second,
                                      std::ios_base::in
                                     );
            // parse the matrix
            peptide_scoring
                ::parse_peptide_scores( my_data,
                                        matrix_name_list
                                      );

        // parse the peptide scores
            std::vector<std::string> peptide_name_list;
            parse_namelist( peptide_name_list, names_list );
        // filter the data, assign to the scores and peptide name list
            if( use_peptide_names )
                {
                    my_data.scores = my_data.scores.transpose();
                }
            
            my_data.scores = my_data.scores.filter_rows( peptide_name_list );
            my_data.pep_names = peptide_name_list;
        }

    std::ofstream output( s_opts->out_matrix_fname, std::ios_base::out );


    peptide_score_data_sample_major& joined_data =
        parsed_score_data[ 0 ];

    // join all of the matrices
    for( uint idx = 1; idx < parsed_score_data.size(); ++idx )
        {
            joined_data.scores = join_with_resolve_strategy( joined_data,
                                                             parsed_score_data[ idx ],
                                                             s_opts->
                                                             duplicate_resolution_strategy
                                                             );

        }

    if( use_peptide_names )
        {
            output << joined_data.scores;
        }
    else
        {
            output << joined_data.scores.transpose();
        }

    time.stop();

    std::cout << "Took " << time.get_elapsed() << " seconds.\n";
}

labeled_matrix<double,std::string>
    module_subjoin::join_with_resolve_strategy( peptide_score_data_sample_major first,
                                                peptide_score_data_sample_major second,
                                                evaluation_strategy::duplicate_resolution_strategy
                                                resolution_strategy
                                                ) const
{
    std::unordered_set<std::string> row_intersection;
    std::unordered_set<std::string> col_intersection;

    setops::set_intersection( row_intersection,
                              first.scores.get_row_labels(),
                              second.scores.get_row_labels(),
                              setops
                              ::get_key<std::string,
                              std::uint32_t
                              >()
                            );

    if( !row_intersection.empty() )
        {
            std::cout << "WARNING: Duplicate names have been encountered, duplicates will "
                         "be resolved with the '"
                      << evaluation_strategy::to_string( resolution_strategy )
                      << "' duplicate resolution strategy.\n";
            setops::set_intersection( col_intersection,
                                      first.scores.get_col_labels(),
                                      second.scores.get_col_labels(),
                                      setops
                                      ::get_key<
                                      std::string,std::uint32_t
                                      >()
                                    );

            if( evaluation_strategy::duplicate_resolution_strategy::COMBINE
                == resolution_strategy
              )
                {
                    for( const auto& row_duplicate : row_intersection )
                        {
                            for( const auto& col_duplicate : col_intersection )
                                {
                                    double first_val = first.scores( row_duplicate, col_duplicate );
                                    second.scores( row_duplicate, col_duplicate ) += first_val;
                                }
                        }

                }
            else if( evaluation_strategy::duplicate_resolution_strategy::INCLUDE
                     == resolution_strategy
                   )
                {
                    for( const auto& row_duplicate : row_intersection )
                        {
                            first.scores.set_row_label( row_duplicate,
                                                        row_duplicate + "_" + first.file_name
                                                      );
                            second.scores.set_row_label( row_duplicate,
                                                         row_duplicate + "_" + second.file_name
                                                      );

                        }
                }
        }

    return first.scores.full_outer_join( second.scores );
}

std::string module_subjoin::get_name()
{
    return "Subjoin";
}
