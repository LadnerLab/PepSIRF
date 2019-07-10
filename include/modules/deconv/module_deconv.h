#ifndef MODULE_DECONV_HH_INCLUDED 
#define MODULE_DECONV_HH_INCLUDED 
#include <string>
#include <vector>
#include <algorithm>

#include "module.h"
#include "options_deconv.h"
#include "maps.h"

class module_deconv : public module
{
 public:
    module_deconv();
    std::string get_name();
    void run( options *opts );

    parallel_map<std::size_t, int>
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
