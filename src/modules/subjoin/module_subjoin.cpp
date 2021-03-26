#include <unordered_set>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
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
                    if( split_line.size() == 2 && split_line[ 1 ].length() > 0 
                    
                    )
                        {
                            dest.emplace_back( split_line[ 0 ] );
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
    // join matrix name pairs provided by potentially both/either a multifile input and command input.
    std::vector<std::pair<std::string,std::string>> matrix_name_pairs;
    if( s_opts->input_matrix_name_pairs.size() + s_opts->multi_matrix_name_pairs.size() == 0 )
        {
            throw std::runtime_error( "No -i (--input) or -m (--multi_file) input has been provided. Manipulation cannot proceed.\n" );
        }
    matrix_name_pairs.reserve( s_opts->input_matrix_name_pairs.size() + s_opts->multi_matrix_name_pairs.size() );
    if( !s_opts->multi_matrix_name_pairs.empty() )
        {
            matrix_name_pairs.insert( matrix_name_pairs.end(),
                                      s_opts->multi_matrix_name_pairs.begin(),
                                      s_opts->multi_matrix_name_pairs.end()
                                    );
        }
    if( !s_opts->input_matrix_name_pairs.empty() )
        {
            matrix_name_pairs.insert( matrix_name_pairs.end(),
                                      s_opts->input_matrix_name_pairs.begin(),
                                      s_opts->input_matrix_name_pairs.end()
                                    );
        }
    parsed_score_data.resize( matrix_name_pairs.size() );
    std::uint32_t idx = 0;

    bool use_peptide_names = !s_opts->use_sample_names;

    time_keep::timer time;

    time.start();

    #pragma omp parallel for num_threads( 2 ) private( idx ) schedule( dynamic ) \
            shared( matrix_name_pairs, parsed_score_data )
    for( idx = 0; idx < matrix_name_pairs.size(); ++idx )
        {
            auto &score_name_pair = matrix_name_pairs[ idx ];

            peptide_score_data_sample_major&
                my_data = parsed_score_data[ idx ];

            std::string& matrix_name_list  = score_name_pair.first;
            // parse the matrix
            peptide_scoring::parse_peptide_scores( my_data,
                                                matrix_name_list
                                                );
            if( use_peptide_names )
                {
                    my_data.scores = my_data.scores.transpose();
                }

            // parse the peptide scores
            std::vector<std::string> peptide_name_list;
            name_replacement_list replacement_names;
            if( score_name_pair.second.empty() )
                {
                    std::string line;
                    std::ifstream matrix_names( matrix_name_list,
                                                    std::ios_base::in );
                    std::getline( matrix_names, line );

                    boost::split( peptide_name_list, line, boost::is_any_of( "\t" ) );

                }
            else
                {
                    std::ifstream names_list( score_name_pair.second,
                                            std::ios_base::in
                                            );
                    replacement_names = parse_namelist( peptide_name_list, names_list );
                }
            auto replace_begin = use_peptide_names ? my_data.pep_names.begin()
                                      : my_data.sample_names.begin();

            auto replace_end = use_peptide_names ? my_data.sample_names.end()
                                    : my_data.sample_names.end();

            std::unordered_set<std::string> names;
            std::size_t curr_name_idx;
            if( use_peptide_names )
                {
                    names.insert( my_data.pep_names.begin(),
                                    my_data.pep_names.end()
                                    );
                }
            else
                {
                    names.insert( my_data.sample_names.begin(),
                                    my_data.sample_names.end()
                                );
                }
            // verify given filter names exist
            for( curr_name_idx = 0; curr_name_idx < peptide_name_list.size(); curr_name_idx++ )
                {
                    if( names.find( peptide_name_list[ curr_name_idx ] )  == names.end()
                        && boost::to_lower_copy( peptide_name_list[ curr_name_idx ] ) != "sequence name" )
                        {
                            std::cout << "WARNING: The sample "
                                    << peptide_name_list[ curr_name_idx ]
                                    << " was not found in "
                                    << "the input matrix, "
                                    << "and will not be included in the output.\n";
                            peptide_name_list.erase( peptide_name_list.begin() + curr_name_idx );
                        }
                }
            // filter out unused rows
            my_data.scores = my_data.scores.filter_rows( peptide_name_list );
            my_data.pep_names = peptide_name_list;
            std::vector<std::string> filter_list;
            std::size_t curr_replacement_idx = 0;
            for( curr_name_idx = 0; curr_name_idx < peptide_name_list.size(); curr_name_idx++ )
                {
                    if( names.find( peptide_name_list[ curr_name_idx ] ) != names.end() )
                        {
                            // if the name is found and has a replacement for output, then use the replacement name as the filter.
                            if( curr_replacement_idx < replacement_names.size()
                                && replacement_names[ curr_replacement_idx ].first == peptide_name_list[ curr_name_idx ] )
                                {
                                    filter_list.emplace_back( replacement_names[ curr_replacement_idx ].second );
                                    my_data.scores.set_row_label( replacement_names[ curr_replacement_idx ].first,
                                                            replacement_names[ curr_replacement_idx ].second
                                                        );
                                    std::replace( replace_begin,
                                                replace_end,
                                                replacement_names[ curr_replacement_idx ].first,
                                                replacement_names[ curr_replacement_idx ].second
                                                );
                                    curr_replacement_idx++;
                                }
                            else
                                {
                                    filter_list.emplace_back( peptide_name_list[ curr_name_idx ] );
                                }
                        }
                }
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
