#ifndef FILE_IO_HH_INCLUDED
#define FILE_IO_HH_INCLUDED
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <iostream>
#include <string>

#ifdef ZLIB_ENABLED

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <memory>

#endif

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
     * @note Empty lines in the input are consumed and ignored.
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

        while( std::getline( input, current_line ) )
            {
                boost::algorithm::trim( current_line );
                if( !current_line.empty() )
                    {
                        boost::algorithm::split( split_line, current_line, split );
                        *end = init_pred( split_line.begin(), split_line.end() );
                        ++end;
                    }
            }
    }

    /**
     * Check to see whether an istreamjk is gzipped.
     * @note False positives are positive, as this method simply performs a
     *       check for the 'magic number' (0x1F8B) at the beginning of a gzipped file.
     *       Any file beginning with this magic number will be found to be gzipped
     *       by this function. However, the likelihood that a non-gzipped file,
     *       especially one used by PepSIRF, is rather small.
     * @param to_check The stream that we want to check for gzip-edness.
     * @post istream is ready to begin at the first item in the file.
     * @returns true if to_check is gzipped, false otherwise.
     **/
    bool is_gzipped( std::istream& to_check );

// double check that we can use ZLIB, otherwise
// this is not useful
#ifdef ZLIB_ENABLED

    using gzip_reader_filter =
        boost::iostreams::filtering_streambuf<boost::iostreams::input>;
    using gzip_writer_filter =
        boost::iostreams::filtering_streambuf<boost::iostreams::output>;

    /**
     * A class that is used to read gzipped data from an istream.
     * Provides methods to retrieve data from the internal stream.
     **/
    class gzip_reader
        : public std::istream
    {

        /**
         * Filter that allows us to read gzipped data.
         **/
        gzip_reader_filter read_filter;

    public:

        /**
         * Argument constructor.
         * @param input The input stream to read gzipped data from.
         **/
        gzip_reader( std::istream& input )
            {
                using namespace boost::iostreams;
                read_filter.push( gzip_decompressor() );
                read_filter.push( input );

                this->rdbuf( &read_filter );
            }

        /**
         * Ensure the stream is closed.
         **/
        ~gzip_reader()
            {
                boost::iostreams::close( read_filter );
            }
    };

    /**
     * A class that is used to write gzipped data to an ostream.
     **/
    class gzip_writer : public std::ostream
    {
        /**
         * Filter that is used to compress data written to files.
         **/
        gzip_writer_filter write_filter;

    public:

        /**
         * Construct the writer with an output stream.
         * @param output the stream to which gzipped data will be written to.
         **/
        gzip_writer( std::ostream& output )
            {
                using namespace boost::iostreams;
                write_filter.push( gzip_compressor() );
                write_filter.push( output );
                this->rdbuf( &write_filter );
            }

        /**
         * Destructor to ensure the stream is closed.
         **/
        ~gzip_writer()
            {
                boost::iostreams::close( write_filter );
            }
    };

#endif

}; // namespace pepsirf_io

#endif // FILE_IO_HH_INCLUDED
