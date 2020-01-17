#ifndef FILE_IO_HH_INCLUDED
#define FILE_IO_HH_INCLUDED
#include <fstream>

namespace pepsirf_io
{

template
<typename Iterator,
 typename JoinSequence
>
void write_file( std::ostream& output,
                 Iterator begin,
                 Iterator end,
                 const JoinSequence join
               )
{
    for( auto current = begin; current != end; ++current )
        {
            output << *current << join;
        }
}

}; // namespace pepsirf_io

#endif // FILE_IO_HH_INCLUDED
