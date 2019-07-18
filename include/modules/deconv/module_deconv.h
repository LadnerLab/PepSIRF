#ifndef MODULE_DECONV_HH_INCLUDED 
#define MODULE_DECONV_HH_INCLUDED 
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <memory>

#include "module.h"
#include "options_deconv.h"
#include "sequence.h"
#include "maps.h"

namespace score_method
{
    enum score_strategy
    {
        INTEGER_SCORING = 0,
        FRACTIONAL_SCORING,
        SUMMATION_SCORING
    };
}; //namespace score_method

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
     *          that represent the id of the species that share a kmer 
     *          with the peptide.
     **/
    std::vector<std::pair<std::string,std::vector<std::pair<std::size_t,std::size_t>>>>
        parse_linked_file( std::string fname );


    /**
     * The the score for a certain spec_count map. 
     * This means that spec_count_map contains the peptides that share a 
     * kmer with a certain species. 
     * @param spec_count_map Map that maps peptides to the species that 
     *        a peptides shares an ID with.
     * @param id the species id that is being searched. Note that this is only 
     *        used in summation scoring.
     * @param peptides A list of enriched peptides.
     * @param strat The scoring strategy to use for scoring peptides.
     * @returns The score of the species
     **/
    double get_score( sequential_map<std::string,std::vector<std::pair<std::size_t,std::size_t>>>&
                      spec_count_map,
                      std::size_t id,
                      std::vector<std::string>& peptides,
                      score_method::score_strategy strat
                    );

    /**
     * Parse a map that will provide name->tax id mappings. This map should be formatted 
     * in the same manner as that of 'lineage.dmp' from NCBI. 
     * @param fname The name fo the file to parse
     * @param name_map the destination map that will store the mappings of id->name
     **/
    void
        parse_name_map( std::string fname, std::map<std::size_t,std::string>& name_map );


    /**
     * Write output to a file that will be named 'out_name'
     * @param out_name Name of file that output will be written to.
     * @param out_counts vector of pairs, where the first item in each pair 
     *        is the id of the species, and the second is the count of that 
     *        species. 
     **/
    void write_outputs( std::string out_name,
                        std::map<std::size_t,std::string>*
                        id_name_map,
                        std::vector<std::pair<std::size_t,double>>& out_counts
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
                    std::vector<std::pair<std::string, std::vector<std::pair<std::size_t,std::size_t>>>>&
                    pep_species_vec
                  );

    /**
     * Populate a map with pairings of <peptide name, vector of species ids>
     * entries. 
     * @param pep_id_map Map to populate with entries. This each key in this map
     *        is a peptide name, and each value is the species that this 
     *        peptide shares a 7mer with.
     * @param pep_species_vec a vector containing pairs with the first entry 
     *        a species id, and the second a vector of string peptide names.
     **/
    void pep_to_id( sequential_map<std::string, std::vector<std::pair<std::size_t,std::size_t>>>&
                    pep_id_map,
                    std::vector<std::pair<std::string, std::vector<std::pair<std::size_t,std::size_t>>>>&
                    pep_species_vec
                  );

    /**
     * Count the number of peptides that share a 7mer with 
     * a peptide.
     * @param id_counts vector in which the counts will be stored.
     * @param id_count_map input map containing id->peptides mapping.
     **/
    void count_species( std::vector<std::pair<std::size_t, double>>&
                        id_counts,
                        sequential_map<std::size_t,std::vector<std::string>>&
                        id_count_map,
                        sequential_map<std::string,std::vector<std::pair<std::size_t,std::size_t>>>&
                        spec_count_map,
                        score_method::score_strategy strat
                      );

    /**
     * Choose the 'best' kmers as defined by the scoring options passed to the program.
     **/
    void choose_kmers( options_deconv *opts );


    /**
     * Get the species ID from a name that matches
     * the pattern 'OXX=([0-9]+),([0-9]*),([0-9]*),([0-9])',
     * i.e. 'OXX=' followed by some ids.
     * @param name The name from which to grab the id. 
     * @note The species id is the second group of the 
     *       above regex.
     **/
    std::size_t get_id( std::string name );

    /**
     * Create the linkage file to be used by 'choose_kmers' method.
     **/
    void create_linkage( options_deconv *opts );

    /**
     * Create a map linking kmers to the species that they 
     * appear in. Each kmer in the sequences in 'sequences' is 
     * given a count for each species that the kmer appears in.
     * @note After completion of this function 'map' will contain 
     *       mappings of the form: 'kmer' -> 'species_id' -> 'count'
     * @param map The map that will store mappings of kmer -> species -> count.
     * @param sequences The sequences to analyze.
     **/
    void create_prot_map( sequential_map<std::string,
                          sequential_map<std::size_t,std::size_t>>&
                         map,
                          std::vector<sequence>& sequences,
                          std::size_t k
                       );

    /**
     * Create a map that maps peptides to the 
     * number of times that peptide shares a kmer with
     * a certain species.
     * @param kmer_sp_map Map populated by module_deconv::create_prot_map,
     *        mapping kmers to species identifiers.
     * @param peptide_sp_vec vector to which output will be 
     *        written.
     * @param peptides Vector containing peptides to 
     *        analyze.
     * @param k The kmer size to use. Each peptide in the 
     *        peptides vector will be broken down into its 
     *        component kmers.
     **/
    void create_pep_map( sequential_map<std::string,
                         sequential_map<std::size_t,std::size_t>>&
                         kmer_sp_map,
                         std::vector<std::tuple<std::string,sequential_map<std::size_t,std::size_t>>>&
                         peptide_sp_vec,
                         std::vector<sequence>&
                         peptides,
                         std::size_t k
                       );
    /**
     * Write outputs for the linkage file generation. 
     * @param fname The name of the file to write output to.
     * @param peptide_sp_vec vector containing the information to write 
     *        file output to. Each item in the vector is a tuple with 
     *        the first entry being a string identifying the peptide.
     *        The second maps peptides to a score.
     * @note Writes a header to the file.
     * @note Each peptide will get a line in the file. Each line follows
     *       this format: 
     *       pep_name\tid:score,id:score,id:score
     **/
    void write_outputs( std::string fname,
                        std::vector<std::tuple<std::string,sequential_map<std::size_t,std::size_t>>>&
                        peptide_sp_vec
                      );



    /**
     * Filter counts that do not have a high enough score out of the id_counts vector.
     * @param id_counts Vector containing <species_id, score> pairs.
     * @param thresh The threshold value. Any pairs in id_counts whose second 
     *        entry is strictly less than this value will be removed.
     * @note This method has the side effect of removing items from id_counts
     **/
    void filter_counts( std::vector<std::pair<std::size_t, double>>&
                        id_counts,
                        double thresh
                        );


    /**
     * The the score strategy to use for scoring.
     * Determines which to used based on the values of 
     * fractional_scoring and summation_scoring
     * @param opts A pointer to 'options_deconv' object.
     **/
    score_method::score_strategy
        get_score_method( options_deconv *opts );

};

template <class K, class V>
struct compare_pair
{
    bool operator()( const std::pair<K,V>& first,
                     const std::pair<K,V>& second
                   )
    {
        return std::get<1>( first ) > std::get<1>( second );
    }
};

template <class K, class V>
struct compare_pair_non_decreasing
{
    bool operator()( const std::pair<K,V>& first,
                     const std::pair<K,V>& second
                   )
    {
        return std::get<1>( first ) < std::get<1>( second );
    }
};


#endif // MODULE_DECONV_HH_INCLUDED 
