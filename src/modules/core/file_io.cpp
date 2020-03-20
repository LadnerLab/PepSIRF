#include "file_io.h"

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
