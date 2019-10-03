#include "fs_tools.h"

std::string& fs_tools::to_dir_name( std::string& dest,
                                    const std::string& str
                                  )
{
    dest = str;
    if( str.back() != '/' )
        {
            dest += "/";
        }
    return dest;
}

std::string fs_tools::to_dir_name( const std::string& str )
{
    std::string ret;

    fs_tools::to_dir_name( ret, str );
    return ret;
}

bool fs_tools::has_trailing_slash( const std::string& query )
{
    return query.back() == '/';
}
