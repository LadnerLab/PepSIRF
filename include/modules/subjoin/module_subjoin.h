#ifndef MODULE_SUBJOIN_HH_DEFINED
#define MODULE_SUBJOIN_HH_DEFINED
#include<fstream>

#include "module.h"
#include "options_subjoin.h"
#include "peptide_scoring.h"

/**
 * A module to join multiple score matrices into one.
 * Filters can be included to determine which peptide names will be output.
 **/
class module_subjoin : public module
{
 public:
    module_subjoin();

    /**
     * Run the subjoin module.
     **/
    void run( options *opts );

    /**
     * Parse a list of names from the input stream.
     * There should be one name per line in the stream.
     * @param dest The location to store names, 
     *        this method will store items found in file at 
     *        the end of the vector.
     * @param file A stream to read names from.
     **/
    void parse_namelist( std::vector<std::string>& dest,
                         std::istream& file
                       );


};

#endif // MODULE_SUBJOIN_HH_INCLUDED
