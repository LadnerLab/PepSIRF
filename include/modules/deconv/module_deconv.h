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


namespace tie_data
{
    enum tie_type
    {
        SINGLE_WAY_TIE = 0,
        TWO_WAY_TIE,
        K_WAY_TIE
    };

};
namespace evaluation_strategy
{
    enum score_strategy
    {
        INTEGER_SCORING = 0,
        FRACTIONAL_SCORING,
        SUMMATION_SCORING
    };

    enum filter_strategy
    {
        SCORE_FILTER = 0,
        COUNT_FILTER
    };

    enum tie_eval_strategy
    {
        PERCENT_TIE_EVAL = 0,
        INTEGER_TIE_EVAL,
        SUMMATION_SCORING_TIE_EVAL
    };
}; //namespace evaluation_strategy

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
     * Additionally, each item in the comma-separated column
     * may contain a key-value pair separated by a colon, 
     * where the first item is the id of the species, and the second
     * is the number of 7mers this peptide shares with the species.
     * An example is: peptide_1 TAB 134:22,54:1
     * This means that 22 of peptide_1's kmers are found in species
     * 134, and only one is found in species 4. 
     * If the line is in the form peptide_1 TAB 134,54
     * then each species is assigned a score of 1.
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
                      evaluation_strategy::score_strategy strat
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
                        std::vector<
                        std::tuple<std::size_t,std::size_t,double>
                        >&
                        out_counts
                      );

    /**
     * Write the map detailing which peptides were assigned 
     * to which species. This map is formatted as a tab-delimited file
     * where the first entry in a column is the peptide name, and the second
     * is a comma-delimited list of species this peptide was assigned to.
     * @note The comma-delimited list will only have more than one 
     *       entry in the event of a tie.
     * @param fname The name of the file to write output to
     * @param out_map The map that specifies which peptides were assigned 
     *        to each species. 
     **/
    void
        write_species_assign_map( std::string fname,
                                  sequential_map<std::string,std::vector<std::size_t>>&
                                  out_map
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



    /**
     * Determine how 'overlapped' two species are.
     * Here for species a and b, we define overlap in one 
     * of two was. We can define overlap by the number of shared 
     * peptides by each species or by the percentage of 
     * peptides shared between species.
     * @param id_peptide_map Map relating species ids to the 
     *        peptides that share a kmer with that species.
     * @param first The id of the first species to check. 
     * @param second The id of the second species to check.
     * @param ev_strat Tie evaluation strategy to use. 
     * @pre Both first and second should be a key in id_peptide_map
     * @returns The overlap amount, expressed as either a ratio or 
     *          a count as determined by ev_strat.
     **/
    bool 
        sufficient_overlap( sequential_map<std::size_t,std::vector<std::string>>&
                            id_peptide_map,
                            sequential_map<std::string,sequential_map<std::size_t,std::size_t>>&
                            pep_species_map_wcounts,
                            std::size_t first,
                            std::size_t second,
                            evaluation_strategy::tie_eval_strategy
                            ev_strat,
                            double threshold
                          );


    /**
     * Get and report any potential species that may be tied.
     * For two species to be tied, their scores must be within 
     * score_threshold of each other. 
     * Once two species have been determined to be tied, a tie evaluation 
     * strategy is considered.
     **/
    void
        handle_ties( std::vector<std::pair<std::size_t,double>>&
                     dest_vec,
                     sequential_map<std::size_t, std::vector<std::string>>&
                     id_pep_map,
                     sequential_map<std::string,sequential_map<std::size_t,std::size_t>>&
                     pep_species_map_wcounts,
                     std::vector<std::pair<std::size_t,double>>&
                     tie_candidates,
                     evaluation_strategy::tie_eval_strategy
                     tie_evaluation_strategy,
                     tie_data::tie_type tie_type,
                     double overlap_threshold
                   );

    tie_data::tie_type
        get_tie_candidates( std::vector<std::pair<std::size_t,double>>&
                            candidates,
                            std::vector<std::pair<std::size_t,double>>&
                            scores,
                            double threshold
                          );
    tie_data::tie_type
        get_tie_type( std::size_t to_convert );


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
    void score_species( std::vector<std::pair<std::size_t, double>>&
                        id_counts,
                        sequential_map<std::size_t,std::vector<std::string>>&
                        id_count_map,
                        sequential_map<std::string,std::vector<std::pair<std::size_t,std::size_t>>>&
                        spec_count_map,
                        evaluation_strategy::score_strategy strat
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
    template<template<class,class> class T, class K, class V>
        void filter_counts( T<K,V>& id_counts,
                            V thresh
                          )
    {
        for( auto& item : id_counts )
            {
                if( item.second < thresh )
                    {
                        id_counts.erase( item.first );
                    }
            }
    }

    template<class K, class V>
        void filter_counts( std::vector<std::pair<K,V>> in_vec,
                       V thresh
                    )
    {
        auto it = std::remove_if( in_vec.begin(),
                                  in_vec.end(),
                                  [&]( const std::pair<K,V>& p ) -> bool
                                  {
                                      return p.second < thresh;
                                  }
                                );

        in_vec.erase( it, in_vec.end() );
    }

    
    /**
     * The the score strategy to use for scoring.
     * Determines which to used based on the values of 
     * fractional_scoring and summation_scoring
     * @param opts A pointer to 'options_deconv' object.
     **/
    evaluation_strategy::score_strategy
        get_evaluation_strategy( options_deconv *opts );

    evaluation_strategy::filter_strategy
        get_filter_method( options_deconv *opts );


    /**
     * Get the number of peptides each species shares a 7mer with.
     * @param id_pep_map The map to scan, contains mappings of species id 
     *        to a vector of peptides containing shared kmers. 
     * @param pep_counts Destination map that will contain the counts for each species.
     *        The map will contain one entry per species. The key is the species id
     *        and the value will be the the count.
     **/
    void
        get_species_counts_per_peptide( sequential_map<std::size_t, std::vector<std::string>>&
                                        id_pep_map,
                                        sequential_map<std::size_t,std::size_t>& pep_counts
                                      );


void
handle_kway_tie( std::vector<std::pair<std::size_t,double>>& tie_candidates,
                   sequential_map<std::size_t, std::vector<std::string>>& id_pep_map,
                   evaluation_strategy::tie_eval_strategy
                   eval_strat
                 );
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
