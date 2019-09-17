#ifndef MODULE_DECONV_HH_INCLUDED 
#define MODULE_DECONV_HH_INCLUDED 
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <memory>
#include <unordered_map>

#include "module.h"
#include "options_deconv.h"
#include "sequence.h"
#include "maps.h"
#include "util.h"

// TODO: Remove magic numbers for id index


/**
 * Used for data that is relevant to determining and 
 * communicating data about ties.
 **/
namespace tie_data
{
    /**
     * Enum declaring the types of ties that 
     * can be considered. Here we consider 
     * 3 types.
     **/
    enum tie_type
    {
        /**
         * A single-way tie happens when only one species 
         * meets the criteria that considers ties.
         **/
        SINGLE_WAY_TIE = 0,

        /**
         * A two-way tie happens when exactly two species
         * meet the criteria that considers ties.
         **/
        TWO_WAY_TIE,

        /**
         * A k-way tie happens when at least three species
         * meet the criteria that considers ties.
         **/
        K_WAY_TIE
    };

};

/**
 * Defines certain evaluation strategies for
 * handling different situations. 
 * If we need to score or filter species 
 * we have different ways we can do that. 
 **/
namespace evaluation_strategy
{
    /**
     * Strategy for scoring species. 
     **/
    enum score_strategy
    {
        /**
         * For integer scoring, the score of 
         * each species is the number of peptides 
         * that species shares a kmer with.
         **/
        INTEGER_SCORING = 0,

        /**
         * For fractional scoring,
         * each species is given a score based 
         * on the ratio of species its shared peptides
         * share a kmer with. So a peptide that shares a kmer
         * with 3 species is given a score of 1/3, and a peptide
         * who only shares a kmer with a given species is given a 
         * score of 1/1. In this way peptides that are more unique to
         * a species are given a higher score.
         **/
        FRACTIONAL_SCORING,

        /**
         * Summation scoring is used in conjuction 
         * with peptides that have count information.
         * So if species1 shares 7 kmers with a peptide,
         * then that species is given a score of 7. 
         * Note that this is done for each peptide a species 
         * shares kmers with and the species is given the score
         * of the sum of these scores.
         **/
        SUMMATION_SCORING
    };

    /**
     * Strategies for filtering peptides, 
     * i.e. if either a species' score or 
     * count falls below a threshold that 
     * species is removed from consideration.
     **/
    enum filter_strategy
    {
        /**
         * Species whose scores fall below
         * a certain threshold will be removed.
         **/
        SCORE_FILTER = 0,

        /**
         * Species whose count falls 
         * below a certain threshold will be 
         * removed.
         **/
        COUNT_FILTER
    };

    /**
     * Strategies specifying how to break ties between species. 
     * Generally one of two outcomes will result from a tie:
     * -Two species are tied and have enough overlapped peptides to be reported 
     *  together. 
     * -Two or more species are tied and the species that shares 
     * the smallest overlap with others is reported.
     **/
    enum tie_eval_strategy
    {
        /**
         * Species that share a certain percentage
         * of their kmers are reported together, or 
         * if they do not meet a given threshold they are reported
         * separately.
         **/
        PERCENT_TIE_EVAL = 0,

        /**
         * Species that share a certain number of 
         * peptides are considered tied.
         **/
        INTEGER_TIE_EVAL,

        /**
         * Species must have a number of 
         * unique shared peptides
         **/
        SUMMATION_SCORING_TIE_EVAL
    };

    enum tie_distance_strategy
    {
        INTEGER_DISTANCE = 0,
        RATIO_DISTANCE
    };

}; //namespace evaluation_strategy

/**
 * The species deconvolution module of the 
 * PepSIRF package. This module is used to 
 * reduce ambiguities when there are 
 * multiple similar species whose probes 
 * are considered enriched. 
 **/
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
     * is the number of kmers this peptide shares with the species.
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
     * @TODO: Document type changes
     **/
    std::vector<std::pair<std::string,std::vector<std::pair<std::string,std::size_t>>>>
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

    double get_score( std::unordered_map<std::string,std::vector<std::pair<std::string,std::size_t>>>&
                      spec_count_map,
                      std::string id,
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
        parse_name_map( std::string fname, std::map<std::string,std::string>& name_map );


    /**
     * Write output to a file that will be named 'out_name'
     * @param out_name Name of file that output will be written to.
     * @param out_counts vector of pairs, where the first item in each pair 
     *        is the id of the species, and the second is the count of that 
     *        species. 
     * @param original_scores Mapping of id->( original_count, original_score ) pairs,
     *        used for outputting what the original count/score of each enriched species is.
     *        This is useful for seeing how the scores and counts have changed for each 
     *        of the enriched species.
     **/
    void write_outputs( std::string out_name,
                        std::map<std::string,std::string>*
                        id_name_map,
                        std::vector<
                        std::tuple<std::string,std::size_t,double,bool>
                        >&
                        out_counts,
                        std::unordered_map<std::string,std::pair<std::size_t,double>>&
                        original_scores
                      );

    /**
     * Write the map detailing which peptides were assigned 
     * to which species. This map is formatted as a tab-delimited file
     * where the first entry in a column is the peptide name, the second
     * is a comma-delimited list of species this peptide was assigned to,
     * and the third is the species the peptide originally shared a peptide 
     * with (including that which it was assigned to).
     * @note The comma-delimited list will only have more than one 
     *       entry in the event of a tie.
     * @param fname The name of the file to write output to
     * @param peptide_assign_original Map that relates peptide names to 
     *        the species it originally shared kmers with. Also included 
     *        are the counts that indicate how many kmers that peptide 
     *        shared with the species.
     * @param out_map The map that specifies which peptides were assigned 
     *        to each species. 
     **/
    void
        write_species_assign_map( std::string fname,
                                  std::unordered_map<std::string,std::vector<std::pair<std::string,std::size_t>>>&
                                  peptide_assign_original,
                                  std::unordered_map<std::string,std::vector<std::string>>&
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
     * @param pep_species_map_wcounts Reference to a map that 
     *        maps string peptide names to a map of species that 
     *        share a kmer with that peptide. Each species is also
     *        linked to its score. In summary, we have a mapping 
     *        that looks like this:
     *        peptide->species->score
     * @param first The id of the first species to check. 
     * @param second The id of the second species to check.
     * @param ev_strat Tie evaluation strategy to use. 
     * @param threshold The threshold that is set to determine
     *        whether the overlap between species is sufficient. 
     * @pre Both first and second should be a key in id_peptide_map
     * @returns The overlap amount, expressed as either a ratio or 
     *          a count as determined by ev_strat.
     **/
    bool 
        sufficient_overlap(  std::unordered_map<std::string,std::vector<std::string>>&
                             id_peptide_map,
                             std::unordered_map<std::string,std::unordered_map<std::string,std::size_t>>&
                             pep_species_map_wcounts,
                             std::string first,
                             std::string second,
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
     * @param dest_vec The vector to write tied species to. 
     * @param id_pep_map A map that associates species ids with 
     *        the peptides that it shares a kmer with.
     * @param pep_species_map_wcounts Reference to a map that 
     *        maps string peptide names to a map of species that 
     *        share a kmer with that peptide. Each species is also
     *        linked to its score. In summary, we have a mapping 
     *        that looks like this:
     *        peptide->species->score
     * @param tie_candidates The species whose scores fall within 
     *        a certain threshold of one another. This method 
     *        determines how these species should be handled.
     * @param tie_eval_strategy Strategy to use when evaluating ties.
     *        See evaluation_strategy::tie_eval_strategy for more information.
     * @param tie_type The type of tie, either one, two or k-way ties are 
     *        supported. See tie_data::tie_type for more
     * @param overlap_threshold Threshold that is either an integer number 
     *        or percentage of 
     *        shared kmers for ties to be reported together.
     * @post dest_vec contains one or two species ids that 
     *       are tied. 
     **/
    void
        handle_ties( std::vector<std::pair<std::string,double>>&
                     dest_vec,
                     std::unordered_map<std::string, std::vector<std::string>>&
                     id_pep_map,
                     std::unordered_map<std::string,std::unordered_map<std::string,std::size_t>>&
                     pep_species_map_wcounts,
                     std::vector<std::pair<std::string,double>>&
                     tie_candidates,
                     evaluation_strategy::tie_eval_strategy
                     tie_evaluation_strategy,
                     tie_data::tie_type tie_type,
                     double overlap_threshold
                   );

    /**
     * Determine which species may be tied. 
     * Two species may be tied if their scores are within 
     * threshold of eachother. N species are considered tied 
     * if all of their distances are within threshold of eachother.
     * For example, if threshold is 
     * 0.0, then species are considered tied only if their 
     * scores are exactly the same. 
     * @param candidates The vector to which the candidates will be output.
     *        After this method candidates will contain all species whose  
     *        scores are within threshold of eachother. 
     * @param scores The vector of species and scores to search 
     *        for ties. This vector must be sorted in non-increasing order
     *        of the scores.
     * @param threshold The threshold for species to be enriched.
     *        This method terminates if if reaches a species whose score does 
     *        not reach this threshold.
     * @param ovlp_threshold The score threshold species scores must be within 
     *        for them to considered tied. 
     * @param distance A struct or class that calculates distance 
     *        between scores. For example, it may be good to 
     *        calculate ratios of scores vs. the scores themselves.
     * @note scores MUST be sorted in non-increasing order.
     * @returns The 'type' of tie that is a result. There are three 
     *          types of tie defined by tie_data::tie_type.
     **/
    template<class DistanceCalc>
        tie_data::tie_type
        get_tie_candidates( std::vector<std::pair<std::string,double>>&
                            candidates,
                            std::vector<std::pair<std::string,double>>&
                            scores,
                            double threshold,
                            double ovlp_threshold,
                            DistanceCalc distance
                          )
        {
            double curr_score = 0;
            std::size_t index = 0;

            auto score_diff = [&]( const std::pair<std::string, double>& first,
                                   const std::pair<std::string, double>& second
                                   ) -> double 
                {
                    return distance( first.second, second.second );
                };

            // D( a, a ) = 0
            candidates.push_back( scores[ index ] );

            for( index = 1;
                 index < scores.size()
                     && scores[ index ].second >= threshold;
                 ++index
                 )
                {
                    // remember scores[ 0 ] >= scores[ index ]
                    curr_score = score_diff( scores[ 0 ], scores[ index ] );

                    if( !util::is_integer( ovlp_threshold ) )
                        {
                            if( curr_score >= ovlp_threshold )
                                {
                                    candidates.push_back( scores[ index ] );
                                }
                            else
                                {
                                    break;
                                }
                        }
                    else
                        {

                            // the score of these two species is close
                            // enough to warrant a further check
                            if( curr_score <= ovlp_threshold )
                                {
                                    candidates.push_back( scores[ index ] );
                                }
                            else
                                {
                                    break;
                                }
                        }
                    std::cout << scores[ index ].first << " " << scores[ index ].second << "\n";
                }

            return get_tie_type( candidates.size() );
        }

    /**
     * Determine the way in which 
     * species. There are three types of ties:
     * one-way, two-way and k-way ties. See 
     * tie_data::tie_type for more information
     * @param to_convert The number of species that 
     *        appear to be tied. 
     **/
    tie_data::tie_type
        get_tie_type( std::size_t to_convert );


    /**
     * Populate a map with pairings of <species_id, vector of peptide names> 
     * entries.
     * @param id_pep_map std::unordered_map to populate. Each entry in the map will 
     *        have a key id and a value vector containing the names of the peptides 
     *        this species shares a kmer with.
     * @param pep_species_vec a vector containing pairs with the first entry 
     *        a species id, and the second a vector of string peptide names.
     * @TODO: Document argument changes
     **/
    void id_to_pep( std::unordered_map<std::string, std::vector<std::string>>&
                    id_pep_map,
                    std::vector<std::pair<std::string, std::vector<std::pair<std::string,std::size_t>>>>&
                    pep_species_vec
                  );


    /**
     * Populate a map with pairings of <peptide name, vector of species ids>
     * entries. 
     * @param pep_id_map Map to populate with entries. This each key in this map
     *        is a peptide name, and each value is the species that this 
     *        peptide shares a kmer with.
     * @param pep_species_vec a vector containing pairs with the first entry 
     *        a species id, and the second a vector of string peptide names.
     * @TODO  Document argument changes
     **/
    void pep_to_id( std::unordered_map<std::string, std::vector<std::pair<std::string,std::size_t>>>&
                    pep_id_map,
                    std::vector<std::pair<std::string, std::vector<std::pair<std::string,std::size_t>>>>&
                    pep_species_vec
                  );

    /**
     * Count the number of peptides that share a kmer with 
     * a peptide.
     * @param id_counts vector in which the counts will be stored.
     * @param id_count_map input map containing id->peptides mapping.
     * @param spec_count_map Map that associates peptides with the species 
     *        that share a kmer with the peptide.
     **/
    void score_species( std::vector<std::pair<std::string, double>>&
                        id_counts,
                        std::unordered_map<std::string,std::vector<std::string>>&
                        id_count_map,
                        std::unordered_map<std::string,std::vector<std::pair<std::string,std::size_t>>>&
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
     * @TODO Document id_index
     **/
    std::string get_id( std::string name, std::size_t id_index );

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
     * @param k The kmer size to use when creating the map. 
     *        A species will be linked to a peptide if a peptide shares a
     *        kmer with that species.
     * @TODO: Document id_index
     **/
    void create_prot_map( std::unordered_map<std::string,
                          std::unordered_map<std::string,std::size_t>>&
                          map,
                          std::vector<sequence>& sequences,
                          std::size_t k,
                          std::size_t id_index
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
    void create_pep_map( std::unordered_map<std::string,
                         std::unordered_map<std::string,std::size_t>>&
                         kmer_sp_map,
                         std::vector<std::tuple<std::string,std::unordered_map<std::string,std::size_t>>>&
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
                        std::vector<std::tuple<std::string,std::unordered_map<std::string,std::size_t>>>&
                        peptide_sp_vec
                      );

    /**
     * Filter counts that do not have a high enough score out of the id_counts structure.
     * @param id_counts Structure containing <species_id, score> pairs.
     * @param thresh The threshold value. Any pairs in id_counts whose second 
     *        entry is strictly less than this value will be removed.
     * @note This method has the side effect of removing items from id_counts
     **/
    template<class K, class V>
        void filter_counts( std::unordered_map<K,V>& id_counts,
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

    /**
     * Filter counts that do not have a high enough score out of the id_counts structure.
     * This is a partial specialization for filter_counts.
     * @param id_counts Structure containing <species_id, score> pairs.
     * @param thresh The threshold value. Any pairs in id_counts whose second 
     *        entry is strictly less than this value will be removed.
     * @note This method has the side effect of removing items from id_counts
     **/
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
     * The strategy to use for scoring.
     * Determines which to used based on the values of 
     * fractional_scoring and summation_scoring
     * @param opts A pointer to 'options_deconv' object.
     * @returns The evaluation strategy to use.
     **/
    evaluation_strategy::score_strategy
        get_evaluation_strategy( options_deconv *opts );

    /**
     * Get the filter method to use based upon the 
     * options the program was started with.
     * @param opts The options object that 
     *        has been initialized.
     * @returns the Filter method to use.
     **/
    evaluation_strategy::filter_strategy
        get_filter_method( options_deconv *opts );


    /**
     * Determine how ties should be evaluated 
     * based upon the options provided by the user.
     * @param opts Options that have been set by the 
     *        user
     * @returns The strategy to use when evaluating ties.
     **/
    evaluation_strategy::tie_eval_strategy
        get_tie_eval_strategy( options_deconv *opts );


    /**
     * Get the number of peptides each species shares a kmer with.
     * @param id_pep_map The map to scan, contains mappings of species id 
     *        to a vector of peptides containing shared kmers. 
     * @param pep_counts Destination map that will contain the counts for each species.
     *        The map will contain one entry per species. The key is the species id
     *        and the value will be the the count.
     **/
    void
        get_species_counts_per_peptide( std::unordered_map<std::string, std::vector<std::string>>&
                                         id_pep_map,
                                         std::unordered_map<std::string,std::size_t>& pep_counts
                                       );



    /**
     * Handle a k-way tie, where k >= 3. This method will
     * remove all but one item from tie_candidates. 
     * The remaining item is that which has the least integer overlap
     * with the other item in tie_candidates.
     * @param tie_candidates A vector containing 
     *        pairs of species_id:count items. 
     * @param id_pep_map Map that associates species with the 
     *        peptides they share a kmer with.
     **/
    void
        handle_kway_tie( std::vector<std::pair<std::string,double>>& tie_candidates,
                         std::unordered_map<std::string, std::vector<std::string>>& id_pep_map
                       );


    /**
     * Combine the count and the score for species.
     * Creates a mapping of id->( species_count, species_score ) pairs.
     * @param map A map to which output will be written.
     * @param c_map A mapping of species_id->species_count pairs.
     * @param scores A vector of pairs, each of the form (species_id, species_score)
     * @pre The set of keys in c_map must equal the set of 
     *      first items in the pairs of scores.
     **/
    template<template<class, class, class...> class Mtype>
        void combine_count_and_score( Mtype<std::string,std::pair<std::size_t,double>>& map,
                                      Mtype<std::string,std::size_t>& c_map,
                                      std::vector<std::pair<std::string,double>>& scores
                                    )
        {

            for( const auto& score : scores )
                {
                    const auto count = c_map.find( score.first );

                    map.emplace( score.first,
                                 std::make_pair
                                 ( count->second, score.second )
                               );

                }

        }


    /**
     * Either find an item in a map, or return a default-constructed
     * value from that map. This useful if you want to print the value 
     * if it is found in the map, or an empty string otherwise.
     * @param map A map of type I that maps K->V
     * @param search The item to search for in map.
     * @returns Either the value of search in the map, or 
     *          a default-constructed value if not found. 
     **/
    template<template<class, class, class...> class I, class K, class V>
        V get_map_value( const I<K,V>& map,
                         const K& search
                         )
    {
        if( map.find( search ) != map.end() )
            {
                // return the found value
                return map.find( search )->second;
            }
        // default-construct a value
        return V();
    }

    /**
     * Same function as module_deconv::get_map_value, 
     * but this one takes a Pointer to the map instead of 
     * map itself.
     * @note This is equivalent to calling:
     *       get_map_value( *map, search )
     **/
    template<template<class, class, class...> class I, typename K, typename V>
        V get_map_value( const I<K,V>* map,
                         const K& search
                         )
    {
        // just dereference it and pass along to
        // the other template function
        return get_map_value( *map, search );
    }


    /**
     * The base case for to_stream_if.
     * @param stream The stream to write output to (if cond is true)
     * @param arg The argument to output to the stream.
     * @note If recursion has gone past the first call then 
     *       cond is true by the nature of to_stream_if. 
     *       However, we still check because otherwise we have 
     *       an unused parameter.
     **/
    template<typename T>
        void to_stream_if( std::ostream& stream,
                           bool cond,
                           T arg
                         )
        {
            if( cond )
                {
                    stream << arg;
                }
        }

    /**
     * Recursive template function that writes a set of 
     * arguments to a stream if cond is true.
     * @param T the type of the current arg.
     * @param Args the typenames of the remaining
     *        arguments
     * @param stream The stream to write to if cond is true
     * @param cond The condition to check for truth value.
     *        if cond is true, then arg and args will be output
     *        to the stream. Otherwise, nothing will be output.
     * @param arg The current arg to output
     * @param args The remaining args to output,
     *        one will be output at each recursive call of 
     *        this function.
     **/
    template<typename T, typename... Args>
        void to_stream_if( std::ostream& stream,
                           bool cond,
                           T arg,
                           Args... args
                         )
    {
        if( cond )
            {
                stream << arg;
                to_stream_if( stream, cond, args... );
            }
    }

    /**
     * Used to determine whether the tie threshold should be 
     * treated as a ratio. 
     * @param threshold The threshold to treat as a ratio.
     * @note this is equivalent to calling "!util::is_integer( threshold )"
     * @returns boolean true if threshold should be treated as an ratio,
     *          false otherwise
     **/
    bool use_ratio_score_tie_thresh( double threshold );

    /**
     * Used to determine whether the tie overlap threshold should 
     * be treated as a ratio.
     * @param threshold the overlap threshold to consider
     * @note this is equivalent to calling "!util::is_integer( threshold )"
     * @returns boolean true if the value should be treated as a ratio,
     *          false otherwise
     **/
    bool use_ratio_overlap_threshold( double threshold );
};

/**
 * Compare pairs in non-decreasing order,
 * i.e. for x[ i ], x[ j ] when i < j, then 
 * x[ i ].second >= x[ j ].second
 * 
 * @param first the first pair to check
 * @param second the second pair to check
 * @note operator '>' must be defined for type V
 * @returns true if first.second > second.first
 **/
template <class K, class V>
struct compare_pair_non_increasing
{
    bool operator()( const std::pair<K,V>& first,
                     const std::pair<K,V>& second
                   )
    {
        return first.second > second.second;
    }
};

/**
 * Compare pairs in non-increasing order,
 * i.e. for x[ i ], x[ j ] when i < j, then 
 * x[ i ] <= x[ j ]
 * @param first the first pair to check
 * @param second the second pair to check
 * @note operator '<' must be defined for type V
 * @returns true if first.second < second.first
 **/
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

/**
 * Returns euclidean distance between 
 * 1-dimensional points a and b, i.e. a - b.
 * @param a Item of type V
 * @param b Item of type V
 * @returns a - b
 **/
template<class V>
struct difference
{
    V operator()( const V a, const V b )
    {
        return a - b;
    }
};

/**
 * Returns a / b if a < b,
 * b / a otherwise.
 **/
template<class V>
struct ratio
{
    V operator()( const V a, const V b )
    {
        // we return something <= 1.0, so
        // want to make sure the larger is in the
        // denominator
        return a > b ? b / a : a / b;
    }
};

#endif // MODULE_DECONV_HH_INCLUDED 
