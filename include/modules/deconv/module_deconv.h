#ifndef MODULE_DECONV_HH_INCLUDED 
#define MODULE_DECONV_HH_INCLUDED 
#include <string>
#include <vector>
#include <algorithm>
#include <memory>

#include "module.h"
#include "options_deconv.h"
#include "maps.h"

struct species_data
{
public:
    species_data();
    int count;
    std::vector<std::string> 
    peptides;

};

class module_deconv : public module
{
 public:

    /**
     * Default constructor.
     **/
    module_deconv();

    /**
     * Get the name of the module
     * @returns the string 'Deconv'
     **/
    std::string get_name();

    /**
     * Run the 'deconv' module.
     * @param opts pointer to options_deconv 
     *        configuration for this run.
     **/
    void run( options *opts );

    /**
     * Parses a 'linked' file. 
     * This file should consist of two tab-separated columns
     * where the first column is the name of the peptide, 
     * and the second column contains comma-separated values.
     * Each value in the comma-separated column should be 
     * the species id of a species that shares a kmer with the 
     * peptide in the first tab-delimited column.
     * @param fname Name of file to parse.
     * @returns vector of pairs, where the first entry is the name 
     *          of the peptide, and the second a vector of size_t ids 
     *          that represent the id of the species that share a 7mer 
     *          with the peptide.
     **/
    std::vector<std::pair<std::string,std::vector<std::size_t>>>
        parse_linked_file( std::string fname );

    /**
     *
     **/
    int get_score( std::size_t size );

    /**
     * Write output to a file that will be named 'out_name'
     * @param out_name Name of file that output will be written to.
     * @param out_counts vector of pairs, where the first item in each pair 
     *        is the id of the species, and the second is the count of that 
     *        species. 
     **/
    void write_outputs( std::string out_name,
                        std::vector<std::pair<std::size_t,std::size_t>>& out_counts
                      );

    /**
     * Reads the list of enriched peptides from a file.
     * @param f_name string filename to read
     * @returns sequential set, where each entry in the set 
     *          is the name of an enriched peptide.
     * @note Currently the file cannot have a header,
     *       it is assumed that the first line in the file is the 
     *       first enriched peptide.
     **/
    sequential_set<std::string>
        parse_enriched_file( std::string f_name );

 private:

    /**
     * Populate a map with pairings of <species_id, vector of peptide names> 
     * entries.
     * @param id_pep_map sequential_map to populate. Each entry in the map will 
     *        have a key id and a value vector containing the names of the peptides 
     *        this species shares a 7mer with.
     * @param pep_species_vec a vector containing pairs with the first entry 
     *        a species id, and the second a vector of string peptide names.
     **/
    void id_to_pep( sequential_map<std::size_t, std::vector<std::string>>&
                    id_pep_map,
                    std::vector<std::pair<std::string, std::vector<std::size_t>>>&
                    pep_species_vec );

    /**
     * Populate a map with pairings of <peptide name, vector of species ids>
     * entries. 
     * @param pep_id_map Map to populate with entries. This each key in this map
     *        is a peptide name, and each value is the species that this 
     *        peptide shares a 7mer with.
     * @param pep_species_vec a vector containing pairs with the first entry 
     *        a species id, and the second a vector of string peptide names.
     **/
    void pep_to_id( sequential_map<std::string, std::vector<std::size_t>>&
                    pep_id_map,
                    std::vector<std::pair<std::string, std::vector<std::size_t>>>&
                    pep_species_vec
                  );

    /**
     * Count the number of peptides that share a 7mer with 
     * a peptide.
     * @param id_counts vector in which the counts will be stored.
     * @param id_count_map input map containing id->peptides mapping.
     **/
    void count_species( std::vector<std::pair<std::size_t, std::size_t>>&
                        id_counts,
                        sequential_map<std::size_t,std::vector<std::string>>&
                        id_count_map
                      );

/**
 * Filter counts that do not have a high enough score out of the id_counts vector.
 * @param id_counts Vector containing <species_id, score> pairs.
 * @param thresh The threshold value. Any pairs in id_counts whose second 
 *        entry is strictly less than this value will be removed.
 * @note This method has the side effect of removing items from id_counts
 **/
 void filter_counts( std::vector<std::pair<std::size_t, std::size_t>>& id_counts,
                     std::size_t thresh
                   );


};

template <class K>
struct compare_pair
{
    bool operator()( const std::pair<K,K>& first,
                     const std::pair<K,K>& second
                   )
    {
        return std::get<1>( first ) > std::get<1>( second );
    }
};

#endif // MODULE_DECONV_HH_INCLUDED 
