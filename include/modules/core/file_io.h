#ifndef FILE_IO_HH_INCLUDED
#define FILE_IO_HH_INCLUDED
#include <fstream>

namespace pepsirf_io
{

    /**
     * Simple method to write joined items in the range [begin, end) 
     * to an output stream.
     * @tparam Iterator the iterator to be iterated over.
     * @tparam JoinSequence The sequence that will be written to 
     *         the stream between each item in the range [begin, end)
     * @param output The ostream being written to
     * @param begin The first item in the range
     * @param end the last item in the range
     * @param join The JoinSequence on which items in the range should be 
     *        joined.
     **/
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
