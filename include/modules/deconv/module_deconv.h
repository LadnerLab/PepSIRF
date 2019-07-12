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
    module_deconv();
    std::string get_name();
    void run( options *opts );

    std::vector<std::pair<std::string,std::vector<std::size_t>>>
        parse_linked_file( std::string fname );

    int get_score( std::size_t size );

 private:

    void id_to_pep( sequential_map<std::size_t, std::vector<std::string>>&
                    id_pep_map,
                    std::vector<std::pair<std::string, std::vector<std::size_t>>>&
                    pep_species_vec );

    void pep_to_id( sequential_map<std::string, std::vector<std::size_t>>&
                    pep_id_map,
                    std::vector<std::pair<std::string, std::vector<std::size_t>>>&
                    pep_species_vec
                  );
void count_species( std::vector<std::pair<std::size_t, std::size_t>>&
                    id_counts,
                    sequential_map<std::size_t,std::vector<std::string>>&
                    id_count_map
                    );

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
