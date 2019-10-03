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

    template<typename T>
        std::string& create_fname( std::string& dest,
                                   const std::string& path_base,
                                   T next
                                 )
        {
            std::string path_cpy = path_base;

            if( dest.length() == 0 )
                { 
                    if( !has_trailing_slash( path_cpy ) )
                        {
                            to_dir_name( path_cpy, path_base );
                        }
                    dest += path_cpy;
                }

            dest += boost::lexical_cast<std::string>( next );

            return dest;
        }

    template<typename T, typename... Args>
        std::string& create_fname( std::string& dest,
                                   const std::string& path_base,
                                   T next,
                                   Args... remaining 
                                 )
        {
            std::string path_cpy = path_base;
            if( dest.length() == 0 )
                { 
                    if( !has_trailing_slash( path_cpy ) )
                        {
                            to_dir_name( path_cpy, path_base );
                        }
                    dest += path_cpy;
                }

            dest += boost::lexical_cast<std::string>( next );

            return create_fname( dest, path_cpy, remaining... );
        }

};


#endif // FS_TOOLS_HH_INCLUDED
