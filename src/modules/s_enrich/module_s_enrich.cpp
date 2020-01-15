#include <iostream>
#include <unordered_set>
#include <stdexcept>
#include <numeric>

#include "module_s_enrich.h"
#include "time_keep.h"
#include "omp_opt.h"
#include "fs_tools.h"
#include "predicate.h"

void module_s_enrich::run( options *opts )
{
    options_s_enrich *e_opts = (options_s_enrich*) opts;
    time_keep::timer timer;
    timer.start();


    peptide_score_data_sample_major zscores;
    peptide_score_data_sample_major norm_scores;
    peptide_score_data_sample_major raw_scores;

    peptide_score_data_sample_major *zscores_ptr = &zscores;
    peptide_score_data_sample_major *norm_scores_ptr = &zscores;
    peptide_score_data_sample_major *raw_scores_ptr = nullptr;

    peptide_scoring::parse_peptide_scores( zscores, e_opts->in_zscore_fname );
    peptide_scoring::parse_peptide_scores( norm_scores, e_opts->in_norm_score_fname );

    bool raw_counts_included = !e_opts->in_raw_count_fname.empty();

    if( raw_counts_included )
        {
            peptide_scoring::parse_peptide_scores( raw_scores,
                                                   e_opts->in_raw_count_fname
                                                 );
            raw_scores_ptr = &raw_scores;
        }

    auto output_path = fs_tools::path( e_opts->out_dirname );
    bool dir_exists = !fs_tools::create_directories( output_path );

    if( dir_exists )
        {
            std::cout << "WARNING: the directory '" << e_opts->out_dirname
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

    auto zscore_enriched = [=]( const peptide_score<std::string> val ) -> bool
        {
            return val.zscore >= e_opts->min_zscore;
        };

    auto norm_score_enriched = [=]( const peptide_score<std::string> val ) -> bool
        {
            return val.norm_score >= e_opts->min_norm_score;
        };

    auto enriched = [=]( const peptide_score<std::string> val ) -> bool
    {
        return predicate::value_constrained_by( val,
                                                zscore_enriched,
                                                norm_score_enriched
                                              );
    };

    for( std::size_t sample_idx = 0;
         sample_idx < zscores.scores.nrows();
         ++sample_idx
       )
        {
            auto sample_id = norm_scores.sample_names[ sample_idx ];

            auto enrichment_candidates = get_enrichment_candidates( zscores_ptr,
                                                                    norm_scores_ptr,
                                                                    raw_scores_ptr,
                                                                    sample_id
                                                                  );

            std::vector<peptide_score<std::string>> enriched_probes;

            // only compute the raw_sum if necessary
            double raw_sum = raw_counts_included
                              ? std::accumulate( enrichment_candidates.begin(),
                                                 enrichment_candidates.end(),
                                                 static_cast<double>( 0 ),
                                                 []( const double a,
                                                     const peptide_score<std::string>& b
                                                     ) -> double
                                                 { return a + b.raw_score; }
                                               )
                              : static_cast<double>( 0.0 );

            // if raw_counts is not included, we are automatically enriched
            // otherwise, we must check that the column sum is at least the
            // specified threshold
            bool colsum_enriched = !raw_counts_included
                                   || ( raw_counts_included
                                         && raw_sum >= e_opts->min_raw_count
                                      );

            // always true if raw_counts not included
            if( colsum_enriched )
                {


                    predicate::valid_for( enrichment_candidates.begin(),
                                          enrichment_candidates.end(),
                                          std::back_inserter( enriched_probes ),
                                          enriched
                                        );

                    if( !enriched_probes.empty() )
                        {
                            std::string outf_name = e_opts->out_dirname + '/'
                                + zscores.sample_names[ sample_idx ]
                                + e_opts->out_suffix;

                            std::ofstream out_file{ outf_name, std::ios_base::out };

                            write_probe_names( out_file,
                                               enriched_probes
                                               );
                        }
                }
        }
    
    timer.stop();
    std::cout << "Took " << timer.get_elapsed() << " seconds.\n";
}

std::vector<peptide_score<std::string>>
module_s_enrich::get_enrichment_candidates( const peptide_score_data_sample_major *zscore_data,
                                            const peptide_score_data_sample_major *norm_score_data,
                                            const peptide_score_data_sample_major *raw_score_data,
                                            const std::string sample_id
                                          )
{
    std::vector<peptide_score<std::string>> ret_val;
    std::size_t pep_idx = 0;

    for( pep_idx = 0; pep_idx < zscore_data->pep_names.size(); ++pep_idx )
        {
            std::string pep_name = zscore_data->pep_names[ pep_idx ];

            double pep_zscore = zscore_data->scores( sample_id, pep_name );
            double pep_norm_score = norm_score_data->scores( sample_id, pep_name );
            double pep_raw_score = 0;

            if( raw_score_data != nullptr )
                {
                    pep_raw_score = raw_score_data->scores( sample_id, pep_name );
                }

            ret_val.emplace_back( peptide_score<std::string>( pep_name,
                                                              pep_zscore,
                                                              pep_norm_score,
                                                              pep_raw_score
                                                            )
                                );
        }

    return ret_val;
}

void module_s_enrich::write_probe_names( std::ostream &stream,
                                         const std::vector<peptide_score<std::string>>& probes
                                       )
{
    for( auto probe = probes.begin();
         probe != probes.end();
         ++probe
       )
        {
            stream << probe->peptide;
            stream << "\n";
        }
}
