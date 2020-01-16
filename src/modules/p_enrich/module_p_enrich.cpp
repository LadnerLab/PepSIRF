#include <boost/algorithm/string.hpp>
#include <stdexcept>
#include <numeric>
#include <fstream>
#include <iostream>

#include "module_p_enrich.h"
#include "time_keep.h"
#include "omp_opt.h"
#include "fs_tools.h"
#include "predicate.h"


void module_p_enrich::run( options *opts )
{
    options_p_enrich *p_opts = (options_p_enrich*) opts;
    time_keep::timer timer;
    timer.start();


    peptide_score_data_sample_major zscores;
    peptide_score_data_sample_major norm_scores;
    peptide_score_data_sample_major raw_scores;

    peptide_score_data_sample_major *zscores_ptr = &zscores;
    peptide_score_data_sample_major *norm_scores_ptr = &zscores;
    peptide_score_data_sample_major *raw_scores_ptr = nullptr;

    peptide_scoring::parse_peptide_scores( zscores, p_opts->in_zscore_fname );
    peptide_scoring::parse_peptide_scores( norm_scores, p_opts->in_norm_scores_fname );

    std::ifstream pairs_file;

    pairs_file.exceptions( std::ios::badbit
                           | std::ios::failbit
                         );

    try
        {
            pairs_file.open( p_opts->in_samples_fname );
        }

    catch( const std::ifstream::failure& e )
        {
            throw std::runtime_error( "Unable to open the provided samples file" );
        }
    
    auto sample_pairs = parse_samples( pairs_file );

    pairs_file.close();

    bool raw_counts_included = !p_opts->in_raw_scores_fname.empty();

    auto output_path = fs_tools::path( p_opts->out_dirname );
    bool dir_exists = !fs_tools::create_directories( output_path );

    if( dir_exists )
        {
            std::cout << "WARNING: the directory '" << p_opts->out_dirname
                      << "' exists, any files with " 
                      << "colliding filenames will be overwritten!\n";
        }

    using sample_name_set = std::unordered_set<std::string>;

    sample_name_set zscore_sample_names{ zscores.sample_names.begin(),
            zscores.sample_names.end()
            };

    sample_name_set norm_score_sample_names{ norm_scores.sample_names.begin(),
            norm_scores.sample_names.end()
            };

    sample_name_set raw_score_sample_names{ raw_scores.sample_names.begin(),
            raw_scores.sample_names.end()
            };

    // ensure zscore_sample_names and norm_score_sample_names are
    // equal. If included, raw_count_sample_names must also equal
    // zscore equal names, otherwise we do not care.
    if( zscore_sample_names != norm_score_sample_names
        && ( !raw_counts_included
             || ( zscore_sample_names == raw_score_sample_names
                  )
             )
        )
        {
            throw std::runtime_error( "The samplenames provided in each input file are "
                                      "not the same"
                                    );
        }


}


std::vector<module_p_enrich::sample_type>
module_p_enrich::parse_samples( std::istream& file )
{
    std::vector<sample_type> return_val;

    std::string current_line;
    std::vector<std::string> values;
    // reserve the expected amount
    values.reserve( 2 );

    while( std::getline( file, current_line ).good() )
        {
            boost::split( values,
                          current_line,
                          boost::is_any_of( "\t" )
                       );

            if( values.size() != 2 )
                {
                    throw std::runtime_error( "The input samples file must "
                                              "contain two samples per line "
                                            );
                }

            return_val.emplace_back( std::make_pair( values[ 0 ],
                                                     values[ 1 ]
                                                   )
                                   );
        }

    return return_val;
}
