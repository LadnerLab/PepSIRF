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
