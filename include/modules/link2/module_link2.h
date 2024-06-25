#ifndef MODULE_LINK2_HH_INCLUDED
#define MODULE_LINK2_HH_INCLUDED
#include "module.h"
#include "options_link2.h"
#include "scored_entity.h"
#include "sequence.h"
#include "metadata_map2.h"
#include "time_keep.h"
#include <regex>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "omp_opt.h"


class module_link2 : public module
{
 public:
    const std::vector<char> FILTER = {'X', 'B', 'J', 'Z'};

    module_link2();
    /**
     * Run the link module, passing the options
     * specified by a user.
     **/
    void run( options *opts );

    /**
     * Write outputs for the linkage file generation.
     * @param fname The name of the file to write output to.
     * @param peptide_sp_vec vector containing the information to write
     *        file output to.
     * @note Writes a header to the file.
     * @note Each peptide will get a line in the file. Each line follows
     *       this format:
     *       pep_name\tid:score,id:score,id:score
     **/
    void write_outputs( std::string fname, 
            const std::unordered_map<std::string, std::unordered_map<std::string, int>> out_scores );

    std::vector<std::vector<std::size_t>> generate_patterns( std::size_t kmer_size, std::size_t span );

    void get_patterned_kmers( std::unordered_set<std::string> &target_kmers, std::string seq, 
                                 std::vector<std::size_t> patt, std::vector<char> filter
                             );
};


#endif // MODULE_LINK2_HH_INCLUDED
