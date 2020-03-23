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

#ifdef ZLIB_ENABLED

std::ostream& pepsirf_io::gzip_writer::get_raw_stream()
{
    return *stream_ptr;
}


namespace std
{
    pepsirf_io::gzip_reader& getline( pepsirf_io::gzip_reader& zip, string& s )
    {
        char str[ 10000 ];
        zip.getline( str, 10000 );

        s = str;
        return zip;
    }

};


#endif
