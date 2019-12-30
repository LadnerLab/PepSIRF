#ifndef PEPTIDE_BIN_HH_INCLUDED
#define PEPTIDE_BIN_HH_INCLUDED
#include <iostream>
#include <vector>

class peptide_bin
{
 public:


};

class bin_collection
{
 public:

    /**
     * Actual store of our bins
     **/
    std::vector<peptide_bin> bins;

    using iterator = typename std::vector<peptide_bin>::iterator;
    using const_iterator = typename std::vector<peptide_bin>::const_iterator;

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    /**
     * Construct a bin collection from items in the range
     * [begin, end).
     * @tparam Iterator the type of iterator to read from
     * @param begin The first item in the range
     * @param end the (exclusive) last item in the range
     **/
    template<typename Iterator>
        bin_collection( const Iterator begin,
                        const Iterator end
                      )
        : bins{ begin, end } {}
        

};

namespace peptide_bin_io
{
    bin_collection parse_bins( std::istream& bin_source );
    void write_bins( std::ostream& dest, const bin_collection& bins );
};

#endif // PEPTIDE_BIN_HH_INCLUDED
