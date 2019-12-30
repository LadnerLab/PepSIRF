#ifndef PEPTIDE_BIN_HH_INCLUDED
#define PEPTIDE_BIN_HH_INCLUDED
#include <iostream>
#include <vector>
#include <unordered_set>

class peptide_bin
{
 public:

    std::unordered_set<std::string> peptide_names;

    using iterator = typename std::unordered_set<std::string>::iterator;
    using const_iterator = typename std::unordered_set<std::string>::const_iterator;
    using value_type = std::string;

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    void add_peptide( const value_type to_add );

    /**
     * Add each peptide in the range [begin,end) to the bin.
     **/
    template<typename Iterator>
        void add_peptides( const Iterator begin,
                           const Iterator end
                         )
        {
            for( auto x = begin; x != end; ++x )
                {
                    add_peptide( *x );
                }
        }

    /**
     * Construct a peptide bin from items in the range
     * [begin, end).
     * @tparam Iterator the type of iterator to read from
     * @param begin The first item in the range
     * @param end the (exclusive) last item in the range
     **/
    template<typename Iterator>
        peptide_bin( const Iterator begin,
                     const Iterator end
                   )
        : peptide_names{ begin, end } {}

    peptide_bin();
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
    using value_type = peptide_bin;

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    /**
     * Add a single bin to the bin_collection.
     **/
    void add_bin( const peptide_bin bin );

    /**
     * Add each bin in the range [begin,end) to the
     * collection.
     **/
    template<typename Iterator>
        void add_bins( const Iterator begin,
                       const Iterator end
                     )
        {
            for( auto x = begin; x != end; ++x )
                {
                    add_bin( *x );
                }
        }

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

    bin_collection();
        

};

namespace peptide_bin_io
{
    /**
     *
     **/
    bin_collection parse_bins( std::istream& bin_source );

    /**
     *
     **/
    void write_bins( std::ostream& dest, const bin_collection& bins );
};

#endif // PEPTIDE_BIN_HH_INCLUDED
