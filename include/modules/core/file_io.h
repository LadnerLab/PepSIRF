#ifndef FILE_IO_HH_INCLUDED
#define FILE_IO_HH_INCLUDED
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <string>

#ifdef ZLIB_ENABLED

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <memory>
#include <iostream>

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

// double check that we can use ZLIB, otherwise
// this is not useful
#ifdef ZLIB_ENABLED

    using gzip_reader_filter =
        boost::iostreams::filtering_streambuf<boost::iostreams::input>;

    /**
     * A class that is used to read gzipped data from an istream.
     * Provides methods to retrieve data from the internal stream.
     **/
    class gzip_reader
        : public std::istream
    {

        /**
         * Pointer to the istream object that contains the filter that 
         * will be used for reading. 
         **/
        std::unique_ptr<std::istream>
            stream_ptr;

        /**
         * A filter to read gzipped data.
         **/
        gzip_reader_filter filter;

    public:

        /**
         * Construct the gzip reader, given an input stream.
         * @note input does not have to be open at the time of construction, 
         *       but MUST be opened before any reading is attempted.
         * @param input The istream to read from. This will be used internally 
         *        as the source of gzipped content for the gzip_reader.
         **/
        gzip_reader( std::istream& input )
            {
                using namespace boost::iostreams;
                filter.push( gzip_decompressor() );
                filter.push( input );
                stream_ptr = std::unique_ptr<std::istream>( new std::istream( &filter ) );
            }

        /**
         * Get a line from the gzipped stream.
         * @note follows the semantics of istream::getline
         * @note Requires the parameters that are necessary for a call to 
         *       istream::getline.
         **/
        template<typename... Args>
        gzip_reader& getline( Args&&... par )
            {
            this
                ->stream_ptr
                ->getline( std::forward<Args>( par )... );
            return *this;

            }

        /**
         * Get a character from the stream.
         * @note follows the semantics of istream::getline.
         **/
        template<typename... Args>
        gzip_reader& get( Args&&... par )
        {
            return this
                ->stream_ptr
                ->get( std::forward<Args>( par )... );
        }

        /**
         * Extract formatted data from the gzip_reader.
         * @note follows semantics of std::operator>>
         **/
        template<typename Args>
            gzip_reader& operator>>( Args&& par )
            {
                std::operator>>( *(this->stream_ptr),
                                 std::forward<Args>( par )
                               );
                return *this;
            }

    };


#endif 

}; // namespace pepsirf_io

#ifdef ZLIB_ENABLED

namespace std
{
    /**
     * Override the getline method when a gzip_reader is used. 
     * This allows us to ensure the gzip_decompressor filter is used when 
     * reading the data. Otherwise, behavior would be undefined.
     * @note follows semantics of std::getline for regular std::istream.
     *       The only difference is that this method first decompresses the data
     *       that is referred to by the stream
     * @param zip The gzip reader to get a line from.
     * @param s The location to store the unzipped line from the file.
     * @returns zip, the original gzip_reader that was passed in.
     **/
    pepsirf_io::gzip_reader& getline( pepsirf_io::gzip_reader& zip, string& s );
};

#endif

#endif // FILE_IO_HH_INCLUDED
