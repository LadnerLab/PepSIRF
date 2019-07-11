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
    void count_species( sequential_map<std::size_t, std::size_t>&
                        id_count_map,
                        std::vector<std::pair<std::string, std::vector<std::size_t>>>&
                        vector
                      );
};

template <class K>
struct compare_vec
{
    bool operator()( std::vector<K>& first, std::vector<K>& second )
    {
        return first.size() < second.size();
    }
};

#endif // MODULE_DECONV_HH_INCLUDED 
