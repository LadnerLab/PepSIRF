#ifndef FILE_IO_HH_INCLUDED
#define FILE_IO_HH_INCLUDED
#include <fstream>
#include <boost/algorithm/string.hpp>

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

    /**
     * Read a file from an input stream, outputting 
     * parsed values to an Iterator. Can perform type 
     * initialization.
     * @tparam SplitPredicate a predicate that returns true if 
     *         the current character is the split character that 
     *         delimits entries within one line of the file.
     * @tparam IteratorInitializer a function that takes iterators
     *         that form the range [ begin, end ) and initialize a value 
     *         that can be written to the OutputIterator
     * @tparam OutputIterator an iterator to write parsed values to.
     * @param input The stream to take data from
     * @param split a split predicate, as described above.
     * @param init_pred the IteratorInitializer, as described above.
     * @param end the output iterator
     **/
    template
        <
        typename SplitPredicate,
        typename IteratorInitializer,
        typename OutputIterator
       >
    void read_file( std::istream& input,
                    const SplitPredicate split,
                    IteratorInitializer init_pred,
                    OutputIterator end
                  )
    {
        std::vector<std::string> split_line;
        std::string current_line;

        while( std::getline( input, current_line ).good() )
            {
                boost::algorithm::split( split_line, current_line, split );
                *end = init_pred( split_line.begin(), split_line.end() );
                ++end;
            }
    }


}; // namespace pepsirf_io

#endif // FILE_IO_HH_INCLUDED
