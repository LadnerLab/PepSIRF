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
#include "file_io.h"


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

    pairs_file.exceptions( std::ios::badbit );

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

    auto zscore_enriched = [=]( const paired_score& val ) -> bool
        {
            return pair_threshold_met( val.zscore,
                                       p_opts->zscore_params
                                     );
        };

    auto norm_score_enriched = [=]( const paired_score& val ) -> bool
        {
            return pair_threshold_met( val.norm_score,
                                       p_opts->norm_scores_params
                                     );
        };

    auto enriched = [=]( const paired_score& val ) -> bool
        {
            return predicate::value_constrained_by( val,
                                                    zscore_enriched,
                                                    norm_score_enriched
                                                   );
        };

    for( std::size_t sample_idx = 0;
         sample_idx < sample_pairs.size();
         ++sample_idx
       )
        {
            auto enrichment_candidates =
                get_enrichment_candidates( zscores_ptr,
                                           norm_scores_ptr,
                                           raw_scores_ptr,
                                           sample_pairs[ sample_idx ]
                                         );

            std::vector<paired_score> enriched_probes;
            enriched_probes.reserve( enrichment_candidates.size() );

            // only attempt to get the column sum if raw counts
            // have been specified
            std::pair<double,double>
                col_sums{ 0.0, 0.0 };

            if( raw_counts_included )
                {
                    col_sums = get_raw_sums( enrichment_candidates.begin(),
                                             enrichment_candidates.end()
                                           );
                }

            bool raw_count_enriched = raw_counts_included
                ? pair_threshold_met( col_sums, p_opts->raw_scores_params )
                : true;

            if( raw_count_enriched )
                {
                    predicate::valid_for( enrichment_candidates.begin(),
                                          enrichment_candidates.end(),
                                          std::back_inserter( enriched_probes ),
                                          enriched
                                        );

                    if( !enriched_probes.empty() )
                        {
                            std::string outf_name =
                                p_opts->out_dirname + '/'
                                + sample_pairs[ sample_idx ].first
                                + p_opts->out_fname_join
                                + sample_pairs[ sample_idx ].second
                                + p_opts->out_suffix;
                            std::ofstream out_file{ outf_name, std::ios_base::out };

                            pepsirf_io::write_file( out_file,
                                                    enriched_probes.begin(),
                                                    enriched_probes.end(),
                                                    "\n"
                                                   );
                            
                        }
                }
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

std::vector<paired_score>
module_p_enrich::get_enrichment_candidates( const peptide_score_data_sample_major *zscore_data,
                                            const peptide_score_data_sample_major *norm_score_data,
                                            const peptide_score_data_sample_major *raw_score_data,
                                            const std::pair<std::string,std::string> sample_names
                                          )
{
    std::vector<paired_score> candidates;
    std::size_t pep_idx = 0;

    using pair = std::pair<double,double>;

    for( pep_idx = 0; pep_idx < zscore_data->pep_names.size(); ++pep_idx )
        {
            std::string pep_name = zscore_data->pep_names[ pep_idx ];
            
            // get zscore in each of the sample
            pair zscores = { zscore_data->scores( sample_names.first, pep_name ),
                             zscore_data->scores( sample_names.second, pep_name )
                           };

            pair norm_scores = { norm_score_data
                                 ->scores( sample_names.first, pep_name ),
                                 norm_score_data
                                 ->scores( sample_names.second, pep_name )
                               };

            pair raw_scores = { 0.0, 0.0 };

            if( raw_score_data != nullptr)
                {
                    raw_scores = { raw_score_data
                                   ->scores( sample_names.first, pep_name ),
                                   raw_score_data
                                   ->scores( sample_names.second, pep_name )
                                 };
                }

            paired_score candidate{ pep_name,
                                    zscores,
                                    norm_scores,
                                    raw_scores
                                  };

            candidates.emplace_back( candidate );
        }

    return candidates;
}
