#ifndef MODULE_DECONV_HH_INCLUDED 
#define MODULE_DECONV_HH_INCLUDED 
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <memory>
#include <unordered_map>
#include <map>
#include <unordered_set>

#include "overlap_data.h"
#include "module.h"
#include "options_deconv.h"
#include "sequence.h"
#include "maps.h"
#include "util.h"
#include "setops.h"
#include "evaluation_strategy.h"
#include "tie_data.h"
#include "scored_peptide.h"
#include "scored_entity.h"
#include "species_data.h"

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
     *          of the peptide, and the second a vector of string ids 
     *          that represent the id of the species that share a kmer 
     *          with the peptide.
     **/
    std::vector<std::pair<std::string,std::vector<std::pair<std::string,double>>>>
        parse_linked_file( std::string fname );


    /**
     * The the score for a certain spec_count map. 
     * This means that spec_count_map contains the peptides that share a 
     * kmer with a certain species. 
     * @param spec_count_map Map that maps peptides to the species that 
     *        a peptides shares an kmer with.
     * @param id the species id that is being searched. Note that this is only 
     *        used in summation scoring.
     * @param peptides A list of enriched peptides.
     * @param strat The scoring strategy to use for scoring peptides.
     * @returns The score of the species
     **/
    double get_score( std::unordered_map<std::string,std::vector<std::pair<std::string,double>>>&
                      spec_count_map,
                      std::string id,
                      std::vector<std::string>& peptides,
                      evaluation_strategy::score_strategy strat
                    );

    /**
     * Get the score that a single peptide contributes to a species.
     * @param peptide The peptide whose score is found.
     * @param spec_count_map Map that associates peptides to the species 
     *        the peptide shares a kmer with
     * @param id the species id that is being searched. 
     * @param score_strat The scoring strategy to use.
     * @note id is only used for summation scoring
     * @returns The score the id contributes to the species, as defined by 
     *          score_strat
     **/
    double score_peptide_for_species( const peptide& peptide,
                                      std::unordered_map
                                      <std::string,
                                      std::vector<std::pair<std::string,double>
                                      >>&
                                      spec_count_map,
                                      std::string id,
                                      evaluation_strategy::score_strategy score_strat
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
                        std::pair<species_data, bool>
                        >&
                        out_counts,
                        std::unordered_map<std::string,std::pair<double,double>>&
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
                                  std::unordered_map<std::string,std::vector<std::pair<std::string,double>>>&
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
     * of two ways. We can define overlap by the number of shared 
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
     * @returns overlap_data<double>, class storing the overlap of first
     *          and second, and the overlap of second and first.
     **/
    template<template<typename...> class MapType>
        overlap_data<double>
        calculate_overlap( MapType<std::string,std::vector<std::string>>&
                           id_peptide_map,
                           MapType<std::string,MapType<std::string,double>>&
                           pep_species_map_wcounts,
                           std::string first,
                           std::string second,
                           evaluation_strategy::tie_eval_strategy
                           ev_strat
                          )
        {
            auto& first_peptides  = id_peptide_map.find( first )->second;
            auto& second_peptides = id_peptide_map.find( second )->second;
            // reserve max( map[ first ].size, map[ second ].size ) in
            std::size_t intersection_size = 0;

            sequential_set<std::string> intersection;

            // we do either union or intersection, alias
            // the intersection set for clear naming 
            sequential_set<std::string>& set_union = intersection;

            intersection.reserve( std::max( first_peptides.size(),
                                            second_peptides.size()
                                            )
                                  );


            if(  ev_strat
                 != evaluation_strategy
                 ::tie_eval_strategy
                 ::SUMMATION_SCORING_TIE_EVAL
                 )
                {
                    setops::set_intersection( intersection,
                                              first_peptides,
                                              second_peptides
                                              );
                    intersection_size = intersection.size();
                }

            if( ev_strat
                == evaluation_strategy
                ::tie_eval_strategy
                ::INTEGER_TIE_EVAL
                )
                {
                    return overlap_data<double>( (double) intersection_size,
                                                 (double) intersection_size
                                               );
                }

            else if( ev_strat
                     == evaluation_strategy
                     ::tie_eval_strategy
                     ::SUMMATION_SCORING_TIE_EVAL
                     )
                {
                    setops::set_union( set_union,
                                       first_peptides,
                                       second_peptides
                                       );

                    // where 'a' and 'b' are for species first and
                    // species second
                    double a_score, b_score, a_num, a_denom,
                        b_num, b_denom;
                    a_score = b_score = a_num = a_denom =
                        b_num = b_denom = 0;

                    for( const auto& peptide : set_union )
                        {
                    
                            auto pep_species_map = pep_species_map_wcounts
                                .find( peptide )->second;

                            if( pep_species_map.find( first ) != pep_species_map.end() )
                                {
                                    a_score = pep_species_map
                                        .find( first )->second;
                                }

                            if( pep_species_map.find( second ) != pep_species_map.end() )
                                {

                                    b_score = pep_species_map
                                        .find( second )->second;
                                }

                            if( b_score > 0 )
                                {
                                    a_num += a_score;
                                }
                            if( a_score > 0 )
                                {
                                    b_num += b_score;
                                }
                            a_denom += a_score;
                            b_denom += b_score;

                            // reset the scores
                            a_score = b_score = 0;
                        }
                    
                    return overlap_data<double>( util::divide( a_num, a_denom ),
                                                 util::divide( b_num, b_denom )
                                               );
                }

            // percent tie evaluation strategy 
            double a_denom = id_peptide_map.find( first )->second.size();
            double b_denom = id_peptide_map.find( second )->second.size();
            double i_size = static_cast<double>( intersection_size );

            return overlap_data<double>( util::divide( i_size, a_denom ),
                                         util::divide( i_size, b_denom )
                                       );
        }


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
                     std::unordered_map<std::string,std::unordered_map<std::string,double>>&
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
     * @param thresholds Map of the thresholds for each species to be enriched.
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
                            std::unordered_map<std::string, std::size_t> thresholds,
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
                     && scores[ index ].second >= thresholds[scores[ index ].first];
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
     **/
    void id_to_pep( std::unordered_map<std::string, std::vector<std::string>>&
                    id_pep_map,
                    std::vector<std::pair<std::string, std::vector<std::pair<std::string,double>>>>&
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
     **/
    void pep_to_id( std::unordered_map<std::string, std::vector<std::pair<std::string,double>>>&
                    pep_id_map,
                    std::vector<std::pair<std::string, std::vector<std::pair<std::string,double>>>>&
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
                        std::unordered_map<std::string,std::vector<std::pair<std::string,double>>>&
                        spec_count_map,
                        evaluation_strategy::score_strategy strat
                      );

    /**
     * Calculate the score that each peptide contributes to a species.
     * @param dest The map to write outputs to. For each species id,
     *        a vector of scored_peptides will be created. Each entry in 
     *        this vector will contain a scored peptide whose is score is the score 
     *        the peptide contributes to the species.
     * @param id_count_map Map associating species to the peptides they share a kmer with
     * @param spec_count_map Map that associates peptides with the species 
     *        that share a kmer with the peptide.
     * @param strat The scoring strategy to use when calculating scores for peptides.
     **/
    void score_species_peptides(
                   std::unordered_map<std::string,
                   std::vector<scored_peptide<double>>
                   >& dest,
                   std::unordered_map<std::string,std::vector<std::string>>&
                   id_count_map,
                   std::unordered_map<std::string,std::vector<std::pair<std::string,double>>>&
                   spec_count_map,
                   evaluation_strategy::score_strategy strat
                                );

    /**
     * Choose the 'best' kmers as defined by the scoring options passed to the program.
     **/
    void choose_kmers( options_deconv *opts );

    /**
     * Write species ids, names (if included), counts, and scores 
     * of species to a stream. 
     * @param id_name_map Pointer to map associating string species ids
     *        with the name of the species. If a species name is not found 
     *        in this map, the species id will be written in place of 
     *        the name. If this is a null pointer, only species ids will 
     *        be written to the stream.
     * @param score_map An iterable (map, vector of pairs, etc) containing 
     *        first values that are the string species id, and the second being 
     *        a pair consisting with the first item being the count 
     *        for a species, and the second being the species' score.
     * @note This method flushes the buffer after it has finished 
     *       writing to it. 
     **/
    template<
        template<typename...> class MapType,
        template<typename...> class IterType
        >
        void write_scores( std::ostream& stream,
                           MapType<std::string,std::string> const *id_name_map,
                           IterType<std::string,std::pair<double,double>>& scores
                         )
        {
            stream << "Species ID\tSpecies Name\tCount\tScore\n";

            for( auto& species : scores )
                {
                    stream << species.first << "\t";
                    if( id_name_map != nullptr )
                        {
                            stream << get_map_value( (*id_name_map),
                                                     species.first,
                                                     species.first
                                                   );
                        }
                    else
                        {
                            stream << species.first;
                        }

                    stream << "\t"
                           << species.second.first
                           << "\t"
                           << species.second.second
                           << "\n";
                }
            std::flush( stream );
        }

    /**
     * Filter counts that do not have a high enough score out of the id_counts structure.
     * @param id_counts Structure containing <species_id, score> pairs.
     * @param thresholds Map of thresholds for each species. Any pairs in id_counts whose second 
     *        entry is strictly less than this value will be removed.
     * @note This method has the side effect of removing items from id_counts
     **/
    template<class K, class V, class T>
        void filter_counts( std::unordered_map<K,V>& id_counts,
                            std::unordered_map<K,T> thresholds
                          )
    {
        for( auto& item : id_counts )
            {
                if( item.second < thresholds[item.first] )
                    {
                        id_counts.erase( item.first );
                    }
            }
    }

    /**
     * Filter counts that do not have a high enough score out of the id_counts structure.
     * This is a partial specialization for filter_counts.
     * @param id_counts Structure containing <species_id, score> pairs.
     * @param thresholds Map of thresholds for each species. Any pairs in id_counts whose second 
     *        entry is strictly less than this value will be removed.
     * @note This method has the side effect of removing items from id_counts
     **/
    template<class K, class V, class T>
        void filter_counts( std::vector<std::pair<K,V>>& in_vec,
                       std::unordered_map<K,T> thresholds
                    )
    {
        auto it = std::remove_if( in_vec.begin(),
                                  in_vec.end(),
                                  [&]( const std::pair<K,V>& p ) -> bool
                                  {
                                      return p.second < thresholds[p.first];
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
                                         std::unordered_map<std::string,double>& pep_counts
                                       );



    /**
     * Handle a k-way tie, where k >= 3. This method will
     * remove all but one item from tie_candidates. 
     * The remaining item is that which has the least integer overlap
     * with the other item in tie_candidates.
     * @param id_pep_map A map that associates species ids with 
     *        the peptides that it shares a kmer with.
     * @param pep_species_map_wcounts Reference to a map that 
     *        maps string peptide names to a map of species that 
     *        share a kmer with that peptide. Each species is also
     *        linked to its score. In summary, we have a mapping 
     *        that looks like this:
     *        peptide->species->score
     * @param tie_candidates A vector containing 
     *        pairs of species_id:count items. 
     * @param tie_eval_strategy Strategy to use when evaluating ties.
     *        See evaluation_strategy::tie_eval_strategy for more information.
     * @param ovlp_threshold Overlap threshold for species to be considered tied.
     *        The overlaps are calculated by module_deconv::calculate_overlap
     **/
    void
        handle_kway_tie( std::vector<std::pair<std::string,double>>& tie_outputs,
                         std::unordered_map<std::string, std::vector<std::string>>& id_pep_map,
                         std::unordered_map<std::string,std::unordered_map<std::string,double>>&
                         pep_species_map_wcounts,
                         std::vector<std::pair<std::string,double>>& tie_candidates,
                         evaluation_strategy::tie_eval_strategy eval_strat,
                         const double ovlp_threshold
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
        void combine_count_and_score( Mtype<std::string,std::pair<double,double>>& map,
                                      Mtype<std::string,double>& c_map,
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
     * Get the highest scoring peptide for each species found in 
     * species_peptide_scores. Each species is mapped to its 
     * highest-scoring peptide in dest.
     * @param dest The location to store each species with its 
     *        highest-scoring peptide
     * @param species_peptide_scores A map associating species with 
     *        the peptides it shares a kmer with. Each scored_peptide
     *        has a peptide and a score.
     **/
    void get_highest_score_per_species( std::unordered_map<std::string,
                                        scored_peptide<double>>& dest,
                                        const std::unordered_map
                                        <std::string, std::vector<scored_peptide<double>>
                                        >&
                                        species_peptide_scores
                                      );



    /**
     * Either find an item in a map, or return a default-constructed
     * value from that map. This useful if you want to print the value 
     * if it is found in the map, or a default string otherwise.
     * @param map A map of type I that maps K->V
     * @param search The item to search for in map.
     * @param default_return The value to return if search
     *        is not found in map. This defaults to a default-constructed
     *        V.
     * @returns Either the value of search in the map, or 
     *          the item provided in default_return value if not found. 
     **/
    template<template<class...> class I, class K, class V>
        V get_map_value( const I<K,V>& map,
                         const K& search,
                         V default_return = V()
                         )
    {
        if( map.find( search ) != map.end() )
            {
                // return the found value
                return map.find( search )->second;
            }
        // default-construct a value
        return default_return;
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
                         const K& search,
                         V default_return = V()
                         )
    {
        // just dereference it and pass along to
        // the other template function
        return get_map_value( *map, search, default_return );
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

    /**
     * Created a map of threshold for each taxaID from linkage map
     * @param thresh_map to populated with taxID as key and threshold as value
     * @param filename filepath to tab-delimited file
     * @returns map of taxID and threshold
     **/
    void thresh_file_to_map( std::unordered_map<std::string, std::size_t>& thresh_map, std::string filename );
};

#endif // MODULE_DECONV_HH_INCLUDED 
