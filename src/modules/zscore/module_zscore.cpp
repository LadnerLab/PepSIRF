#include "module_zscore.h"
#include "matrix.h"
#include "peptide_scoring.h"
#include "time_keep.h"
#include "stats.h"
#include <iostream>
#include <fstream>


module_zscore::module_zscore() = default;

void module_zscore::run( options *opts )
{
    options_zscore *z_opts = (options_zscore*) opts;
    time_keep::timer time;
    time.start();

    // parse the input
    peptide_score_data_sample_major input;

    peptide_scoring::parse_peptide_scores( input, z_opts->in_fname );

    // transpose the input matrix
    input.scores = input.scores.transpose();

    // calculate the zscores
    for( std::uint32_t row_idx = 0; row_idx < input.scores.nrows(); ++row_idx )
        {
            stats::zscore( input.scores.row_begin( row_idx ),
                           input.scores.row_end( row_idx ),
                           input.scores.row_begin( row_idx )
                         );
        }

    // transpose the output matrix
    input.scores = input.scores.transpose();


    // write output
    peptide_scoring::write_peptide_scores( z_opts->out_fname, input );

    time.stop();

    std::cout << "Took " << time.get_elapsed() << " seconds.\n";
}
