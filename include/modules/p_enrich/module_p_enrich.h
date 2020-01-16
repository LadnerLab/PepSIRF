#ifndef MODULE_P_ENRICH_HH_INCLUDED
#define MODULE_P_ENRICH_HH_INCLUDED
#include "module.h"
#include "options_p_enrich.h"
#include "peptide_scoring.h"
#include <string>

class module_p_enrich : public module
{

public:

    using sample_type = std::pair<std::string,std::string>;
    /**
     * Run the bin module, passing 
     * the options specified by a user.
     **/
    void run( options *opts );

    /**
     * Parse a list of sample names from a stream. 
     * @param file The stream from which samples should be read.
     * @throws std::runtime_exception if the file is not formatted correctly.
     * @pre file.good()
     * @returns a vector of pairs of strings representing the input 
     *          sample name pairs.
     **/
    std::vector<sample_type> parse_samples( std::istream& file );
};

#endif // MODULE_P_ENRICH_HH_ENCLUDED
