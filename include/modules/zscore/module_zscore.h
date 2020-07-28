#ifndef MODULE_ZSCORE_HH_INCLUDED
#define MODULE_ZSCORE_HH_INCLUDED
#include "module.h"
#include "options_zscore.h"
#include "peptide_scoring.h"
#include "peptide_bin.h"
#include "nan_report.h"
#include "matrix.h"
#include "stats.h"

class module_zscore : public module
{
 public:
    module_zscore();
    void run( options *opts );

    /**
     * Calculate the zscores for each peptide in a bin.
     * @tparam OutputIt the type of iterator to write nan_reports to
     * @param peptide_bins The peptide bins against which zscores should
     *        be calculated.
     * @param trim_percent Percent [0.0,100.0] of peptides with the
     *        highest and lowest scores to trim.
     * @param sample_name The name of the sample to calculate the zscores for.
     * @param nan_report_location The place to write nan reports to.
     **/
    template<typename OutputIt>
    void
        calculate_zscores( const bin_collection& peptide_bins,
                           options_zscore* z_opts,
                           labeled_matrix<double,std::string>& scores,
                           const std::string sample_name,
                           OutputIt nan_report_location
                         )
        {
            std::vector<double> current_row_counts;

            std::size_t bin_id = 1;

            for( const auto& bin : peptide_bins )
                {
                    // get the counts for the peptides in this group
                    for( const auto& peptide : bin )
                        {
                            current_row_counts
                                .emplace_back( scores( sample_name, peptide ) );
                        }

                    std::sort( current_row_counts.begin(),
                               current_row_counts.end()
                               );
                    // separate to new function for highest density interval support
                    std::vector<double>::iterator new_begin;
                    std::vector<double>::iterator new_end;
                    if( z_opts->hdi_percent != 0.0 )
                        {
                            std::size_t num_interval = std::round( z_opts->hdi_percent * bin.size() );
                            std::size_t count = 0;
                            std::size_t ref = current_row_counts.at( count + num_interval ) - current_row_counts.at( count );
                            std::size_t index;
                            std::size_t curr_value;
                            for( index = 0; index < bin.size() - num_interval; index++ )
                                {
                                    curr_value = current_row_counts.at( index + num_interval ) - current_row_counts.at( index );
                                    if( curr_value < ref )
                                        {
                                            ref = curr_value;
                                            count = index;
                                        }
                                }
                            new_begin = current_row_counts.begin() + count;
                            new_end = current_row_counts.begin() + count + num_interval;
                        }
                    else
                        {
                            std::size_t trim_length = current_row_counts.size() * ( 0.01 * z_opts->trim_percent );
                            new_begin = current_row_counts.begin() + trim_length;
                            new_end = current_row_counts.end() - trim_length;
                        }
                    //
                    double mean = stats::arith_mean( new_begin,
                                                     new_end
                                                   );

                    double stdev = stats::stdev( new_begin,
                                                 new_end,
                                                 mean
                                                 );


                    for( const auto& peptide : bin )
                        {

                            auto& current_value = scores( sample_name, peptide );
                            double zscore = static_cast<double>( 0.0 );

                            if( static_cast<double>( 0.0 ) == mean )
                                {
                                    current_value = std::numeric_limits<double>::quiet_NaN();
                                    *nan_report_location = nan_report( peptide,
                                                                       sample_name,
                                                                       bin_id
                                                                       );
                                    ++nan_report_location;
                                }
                            else if( stdev == 0.0 && mean != 0.0 )
                                {
                                    current_value = 0.0;
                                }
                            else
                                {
                                    zscore = stats::zscore( current_value, mean, stdev );
                                    current_value = zscore;
                                }
                        }

                    ++bin_id;
                    current_row_counts.clear();
                }
    }

};

#endif // MODULE_ZSCORE_HH_INCLUDED
