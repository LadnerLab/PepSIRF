#ifndef FS_TOOLS_HH_INCLUDED
#define FS_TOOLS_HH_INCLUDED
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>

namespace fs_tools
{
    using namespace boost::filesystem;

    std::string to_dir_name( const std::string& str );

    std::string& to_dir_name( std::string& dest,
                              const std::string& str
                            );
};


#endif // FS_TOOLS_HH_INCLUDED
