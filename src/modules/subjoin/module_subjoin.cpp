#include <unordered_set>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include "logger.h"
#include "module_subjoin.h"
#include "matrix.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include "time_keep.h"

module_subjoin::module_subjoin()
{
    name = "Subjoin";
}

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
                            ret_val.insert( {split_line[ 0 ], split_line[ 1 ]} );
                        }
                    else
                        {
                            dest.emplace_back( line );
                            ret_val.insert( {split_line[ 0 ], split_line[ 0 ]} );
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
            Log::error(
                "No -i (--input) or -m (--multi_file) input has"
                " been provided."
                " Manipulation cannot proceed.\n"
            );
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


    std::vector<std::string> orig_names;
    name_replacement_list replacement_names;

    #pragma omp parallel for num_threads( 2 ) private( idx ) schedule( dynamic ) \
            shared( matrix_name_pairs, parsed_score_data )
    for( idx = 0; idx < matrix_name_pairs.size(); ++idx )
        {
            auto &score_name_pair = matrix_name_pairs[ idx ];

            peptide_score_data_sample_major&
                my_data = parsed_score_data[ idx ];

            std::string& matrix_name_list = score_name_pair.first;
            // parse the matrix
            peptide_scoring::parse_peptide_scores( my_data,
                                                matrix_name_list
                                                );
            if( use_peptide_names )
                {
                    my_data.scores = my_data.scores.transpose();
                }

            // parse the peptide scores
            if( !score_name_pair.second.empty() )
                {
                    std::ifstream names_list( score_name_pair.second,
                                            std::ios_base::in
                                            );
                    if( names_list.fail() )
                        {
                            Log::error(
                                "Unable to open name list '"
                                + score_name_pair.second + "'.\n"
                            );
                        }
                    replacement_names = parse_namelist( orig_names, names_list );
                }

            std::unordered_set<std::string> names;
            
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

            bool exclude = s_opts->exclude_names;
            std::unordered_set<std::string> nonexcluded_names;

            // if exclude, create a names set for the names that are not excluded
            if( exclude )
               {
                nonexcluded_names = names;
               }

            std::size_t curr_name_idx;
            // verify given names from namelist exist and set up nonexluded names
            for( curr_name_idx = 0; curr_name_idx < orig_names.size(); curr_name_idx++ )
                {
                    if( !exclude && names.find( orig_names[ curr_name_idx ] ) == names.end()
                        && boost::to_lower_copy( orig_names[ curr_name_idx ] ) != "sequence name" )
                        {
                            Log::warn(
                                "The sample " + orig_names[curr_name_idx]
                                + " was not found in the input matrix, and"
                                " will not be included in the output.\n"
                            );
                            orig_names.erase( orig_names.begin() + curr_name_idx );
                        }
                    else
                        {
                            // remove the name from the non excluded names
                            nonexcluded_names.erase( orig_names[ curr_name_idx ] );
                        }
                }

            // test for exlude option, update names
            if( exclude )
                {  
                    orig_names.clear();
                    replacement_names.clear();
                    for( std::string name : nonexcluded_names )
                        {
                            orig_names.emplace_back( name );
                            replacement_names.insert( {name, name} );
                        }
                }

            if( !score_name_pair.second.empty() )
                {
                    // filter out unused rows using namelist file.
                    my_data.scores = my_data.scores.filter_rows( orig_names );
                    std::vector<std::string> new_rows;
                    new_rows.resize(orig_names.size());
                    // verify no duplicate names
                    std::unordered_set<std::string> output_names;
                    for( const auto& output_name : replacement_names )
                        {
                            output_names.emplace(output_name.second);
                        }
                    if( output_names.size() != replacement_names.size() )   
                        {
                            Log::error(
                                "Duplicate name found in output names provided"
                                " by '--input' or '--multi_file'. Verify"
                                " name list sample/peptide names do not"
                                " include duplicate names.\n"
                            );
                        }
                    // replace original names with updates for output
                    for( auto& orig_name : my_data.scores.get_row_labels() )
                        {
                            std::pair<std::string,std::uint32_t> new_name
                                = std::make_pair(
                                    replacement_names
                                        .find(orig_name.first)->second,
                                    orig_name.second
                                );
                            new_rows[new_name.second] = new_name.first;
                        }
                    my_data.scores.update_row_labels( new_rows );
                    my_data.sample_names = new_rows;          
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

    Log::info("Took " + std::to_string(time.get_elapsed()) + " seconds.\n");
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
            Log::warn(
                "Duplicate names have been encountered, duplicates will be"
                " resolved with the '"
                + evaluation_strategy::to_string(resolution_strategy)
                + "' duplicate resolution strategy.\n"
            );
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