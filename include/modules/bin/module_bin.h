#ifndef MODULE_BIN_HH_INCLUDED
#define MODULE_BIN_HH_INCLUDED

#include "peptide_scoring.h"
#include "module.h"
#include "options_bin.h"
#include <vector>

class module_bin : public module
{
 public:
    module_bin();

    /**
     * Run the bin module, passing 
     * the options specified by a user.
     **/
    void run( options *opts );

    /**
     * Sum the counts of the rows in a labeled_matrix.
     * @param data The labeled_matrix whose row values should be summed
     * @returns a vector containing the sum of the counts of the rows in the matrix. 
     *          The n'th entry in this vector is the sum of the n'th row in the data 
     *          matrix.
     **/
    std::vector<double>
    sum_counts( const labeled_matrix<double,std::string>& data );


};

#endif // MOUDLE_BIN_HH_INCLUDED
