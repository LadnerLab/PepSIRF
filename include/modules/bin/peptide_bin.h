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
    using size_type = std::unordered_set<std::string>::size_type;

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    void add_peptide( const value_type to_add );

    bool contains( const value_type& look_for );

    /**
     * Determine whether this bin equals other.
     * @param other the bin to compare this with.
     * @returns true if this and other contain the same
     *          peptides, false otherwise.
     **/
    bool operator==( const peptide_bin& other ) const;

    /**
     * Get the size of the bin.
     * @returns the number of peptides in the bin
     **/
    size_type size() const;

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
    using size_type = std::vector<peptide_bin>::size_type;

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
     * Get the size of the bin_collection.
     * @returns the number of bins in the collection.
     **/
    size_type size() const;

    /**
     * Get the smallest bin in the bin_collection
     * @returns the smallest sized bin
     **/
    peptide_bin smallest();

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

    /**
     * Determine whether bin collections are equal.
     * @param other The bin to compare this with.
     * @returns true if every bin in this collection is
     *          equal to every bin in other.
     * @note order of the bins matters.
     **/
    bool operator==( const bin_collection& other ) const;


};

namespace peptide_bin_io
{
    /**
     * Parse bins from the stream bin_source.
     * Each bin is a tab-separated list of peptide names.
     * Bins are separated by newline characters.
     * @param bin_source The stream containing bins
     * @returns A bin_collection containing each bin from bin_source
     **/
    bin_collection parse_bins( std::istream& bin_source );

    /**
     *
     **/
    void write_bins( std::ostream& dest, const bin_collection& bins );
};

#endif // PEPTIDE_BIN_HH_INCLUDED
