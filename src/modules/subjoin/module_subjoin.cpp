#include <unordered_set>
#include "module_subjoin.h"
#include "matrix.h"
#include <iostream>
#include "time_keep.h" 

module_subjoin::module_subjoin() = default;

void module_subjoin::parse_namelist( std::vector<std::string>& dest,
                                     std::istream& file
                                   )
{
    std::string line;
    while( std::getline( file, line ).good() )
        {
            dest.push_back( line );
        }
}

void module_subjoin::run( options *opts )
{
    options_subjoin *s_opts = (options_subjoin*) opts;
    std::vector<peptide_score_data_sample_major>
        parsed_score_data;
    parsed_score_data.resize( s_opts->matrix_name_pairs.size() );
    std::uint32_t idx = 0;

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
            if( s_opts->use_sample_names )
                {
                    my_data.scores = my_data.scores.transpose();
                }
            
            my_data.scores = my_data.scores.filter_rows( peptide_name_list );
            my_data.pep_names = peptide_name_list;
        }

    std::ofstream output( s_opts->out_matrix_fname, std::ios_base::out );

    if( parsed_score_data.size() > 1 )
        {

            labeled_matrix<double,std::string> joined_matrix;
            // join all of the matrices
            for( uint idx = 0; idx < parsed_score_data.size(); ++idx )
                {
                    joined_matrix = joined_matrix
                        .full_outer_join( parsed_score_data[ idx ].scores );
                }

            if( s_opts->use_sample_names )
                {
                    output << joined_matrix.transpose();
                }
            else
                {
                    output << joined_matrix;
                }
        }
    else
        {
        if( s_opts->use_sample_names )
                {
                    output << parsed_score_data[ 0 ].scores.transpose();
                }
        else
            {
                output << parsed_score_data[ 0 ].scores;
            }
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

