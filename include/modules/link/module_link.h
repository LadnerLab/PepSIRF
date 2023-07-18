#ifndef MODULE_LINK_HH_INCLUDED
#define MODULE_LINK_HH_INCLUDED
#include "module.h"
#include "options_link.h"
#include "scored_entity.h"
#include "sequence.h"
#include "metadata_map.h"
#include <regex>
#include <vector>
#include <unordered_set>
#include <unordered_map>


class module_link : public module
{
 public:
    module_link();
    /**
     * Run the link module, passing the options
     * specified by a user.
     **/
    void run( options *opts );

    /**
     * Create a map linking kmers to the species that they
     * appear in. Each kmer in the sequences in 'sequences' is
     * given a count for each species that the kmer appears in.
     * @note After completion of this function 'map' will contain
     *       mappings of the form: 'kmer' -> 'scored_entity, where
     *       the scored_entity represents the score of a species.
     * @param scores_map The map that will store mappings of kmer -> scored_entity.
     * @param sequences The sequences to analyze.
     * @param id Template object, either an unsigned integer or memory address to an unordered map
     * @param k The kmer size to use when creating the map.
     *        A species will be linked to a peptide if a peptide shares a
     *        kmer with that species.
     **/
    template< typename retrieve_id >
    void create_prot_map( std::unordered_map<std::string,
                          std::unordered_set<scored_entity<std::string,double>>>&
                          scores_map,
                          std::vector<sequence>& sequences,
                          std::size_t k,
                          retrieve_id id
                        );

    /**
     * Create a map that maps peptides to the
     * number of times that peptide shares a kmer with
     * a certain species.
     * @param kmer_sp_map Map populated by module_link::create_prot_map,
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
                         std::unordered_set<scored_entity<std::string,double>>>&
                         kmer_sp_map,
                         std::vector<std::tuple<std::string,
                         std::unordered_set<scored_entity<std::string,double>>>>&
                         peptide_sp_vec,
                         std::vector<sequence>&
                         peptides,
                         std::size_t k
                       );
    /**
     * Create a map that scores peptides based
     * on the scores of its component kmers. Here, the score of a kmer is
     * the 1 / ( the number of times the kmer appears in the 'peptides' vector )
     * @param kmer_sp_map Map populated by module_link::create_prot_map,
     *        mapping kmers to species identifiers.
     * @param peptide_sp_vec vector to which output will be
     *        written.
     * @param peptides Vector containing peptides to
     *        analyze.
     * @param k The kmer size to use. Each peptide in the
     *        peptides vector will be broken down into its
     *        component kmers.
     **/
    void create_pep_map_with_kmer_penalty( std::unordered_map<std::string,
                                           std::unordered_set<scored_entity<std::string,double>>>&
                                           kmer_sp_map,
                                           std::vector<std::tuple<std::string,
                                           std::unordered_set<scored_entity<std::string,double>>>>&
                                           peptide_sp_vec,
                                           std::vector<sequence>&
                                           peptides,
                                           std::size_t k
                                         );

    /**
     * Get the species ID from a name that matches
     * the pattern 'OXX=([0-9]+),([0-9]*),([0-9]*),([0-9])',
     * i.e. 'OXX=' followed by some ids.
     * @param name The name from which to grab the id.
     * @param id_index The index (0-based) of id to choose.
     * @note The species id is the second group of the
     *       above regex.
     **/
    std::string get_id( std::string name, std::size_t id_index );

    /**
     * Either an unsigned integer (provided index_id ) or an address (provided metadata file reference)
     * is provided. When given an integer, the get_id function is called, otherwise a metadata map is constructed and function
     * to fill map is called.
     * @param name The name from which to grab the id.
     * @param id_index value used to retrieve species id.
     * @note Overloaded method. Allows for use of unsighned integer or reference to unordered map.
    **/
    std::string verify_id_type( std::string name, std::size_t id_index);

    /**
     * Either an unsigned integer (provided index_id ) or an address (provided metadata file reference)
     * is provided. When given an integer, the get_id function is called, otherwise a metadata map is constructed and function
     * to fill map is called.
     * @param name The name from which to grab the id.
     * @param map unordered map address for reference used to retrieve species id.
     * @note Overloaded method. Allows for use of unsighned integer or reference to unordered map.
    **/
    std::string verify_id_type( std::string name, std::unordered_map<std::string, std::string> *map );

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
                        std::vector<std::tuple<std::string,
                        std::unordered_set<scored_entity<std::string,double>>>>&
                        peptide_sp_vec
                      );
};


#endif // MODULE_LINK_HH_INCLUDED
