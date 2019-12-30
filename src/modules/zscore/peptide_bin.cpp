#include "peptide_bin.h"

// peptide_bin
peptide_bin::iterator peptide_bin::begin() { return peptide_names.begin(); }
peptide_bin::iterator peptide_bin::end() { return peptide_names.end(); }

peptide_bin::const_iterator peptide_bin::begin() const { return peptide_names.begin(); }
peptide_bin::const_iterator peptide_bin::end() const { return peptide_names.end(); }

peptide_bin::peptide_bin() = default;

bool peptide_bin::operator==( const peptide_bin& other ) const
{
    return peptide_names == other.peptide_names;
}

void peptide_bin::add_peptide( value_type to_add )
{
    peptide_names.emplace( to_add );
}

// bin_collection 
bin_collection::iterator bin_collection::begin() { return bins.begin(); }
bin_collection::iterator bin_collection::end() { return bins.end(); }

bin_collection::const_iterator bin_collection::begin() const { return bins.begin(); }
bin_collection::const_iterator bin_collection::end() const { return bins.end(); }

bin_collection::bin_collection() = default;

bool bin_collection::operator==( const bin_collection& other ) const
{
    return bins == other.bins;
}

void bin_collection::add_bin( const peptide_bin bin )
{
    bins.emplace_back( bin );
}

// namespace peptide_bin_io
bin_collection peptide_bin_io::parse_bins( std::istream& bin_source )
{


}
