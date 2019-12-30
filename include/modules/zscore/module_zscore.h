#ifndef MODULE_ZSCORE_HH_INCLUDED
#define MODULE_ZSCORE_HH_INCLUDED
#include "module.h"
#include "options_zscore.h"
#include "peptide_scoring.h"
#include "peptide_bin.h"
#include "matrix.h"

class module_zscore : public module
{
 public:
    module_zscore();
    void run( options *opts );

    void
        calculate_zscores( const bin_collection& peptide_bins,
                           const double trim_percent,
                           labeled_matrix<double,std::string>& scores,
                           const std::uint32_t sample_num
                         );
};

#endif // MODULE_ZSCORE_HH_INCLUDED
