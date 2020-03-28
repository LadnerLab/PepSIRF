#include "file_io.h"

bool pepsirf_io::is_gzipped( std::istream& to_check )
{
    to_check.clear();
    to_check.seekg( 0 );

    unsigned char bytes[ 3 ];

    to_check.get( (char*) bytes, 3 );

    bool is_gzipped = ( bytes[ 0 ] == 0x1F )
        && ( bytes[ 1 ] == 0x8B );

    to_check.clear();
    to_check.seekg( 0 );

    return is_gzipped;
}
