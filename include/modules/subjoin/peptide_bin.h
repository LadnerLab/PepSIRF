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
    std::vector<peptide_bin> bins;
    using iterator = typename std::vector<peptide_bin>::iterator;

    iterator begin() { return bins.begin(); }
    iterator end() { return bins.end(); }

};

namespace peptide_bin_io
{
    bin_collection parse_bins( std::istream& bin_source );
    void write_bins( std::ostream& dest, const bin_collection& bins );
};

#endif // PEPTIDE_BIN_HH_INCLUDED
