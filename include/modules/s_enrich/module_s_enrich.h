#ifndef MODULE_S_ENRICH_HH_INCLUDED
#define MODULE_S_ENRICH_HH_INCLUDED
#include "module.h"
#include "options_s_enrich.h"
#include "peptide_score.h"
#include "peptide_scoring.h"
#include <vector>
#include <fstream>


class module_s_enrich : public module
{
public:

    void run( options *opts );

    /**
     * Get the score data about a peptide, returning each of its scores in 
     * a vector of peptide_score types. 
     * @param zscore_data contains data for peptide zscores.
     * @param norm_score_data Data containing normalized scores for each peptide
     * @param raw_score_data Data containing raw scores for each peptide.
     * @pre zscore_data, norm_score_data, and (if included) raw_score_data 
     *      all contain data for the same peptides and samples.
     * @pre zscore_data and norm_score_data are not null
     * @note if raw_score_data is null, the return peptide_scores will 
     *       have their raw_score member set to 0.
     **/
    std::vector<peptide_score<std::string>>
    get_enrichment_candidates( const peptide_score_data_sample_major *zscore_data,
                               const peptide_score_data_sample_major *norm_score_data,
                               const peptide_score_data_sample_major *raw_score_data,
                               const std::string sample_id
                             );

    /**
     * Write the names of probes to an output stream, one per line..
     * @param stream The stream to write probe names to
     * @param probes A vector of peptide_scores to whose names
     *        should be written to the stream.
     **/
    void write_probe_names( std::ostream &stream,
                            const std::vector<peptide_score<std::string>>& probes
                          );


};

#endif // MODULE_S_ENRICH_HH_INCLUDED 
