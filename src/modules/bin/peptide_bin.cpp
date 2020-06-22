#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/trim.hpp>
#include "peptide_bin.h"

// peptide_bin
peptide_bin::iterator peptide_bin::begin() { return peptide_names.begin(); }
peptide_bin::iterator peptide_bin::end() { return peptide_names.end(); }

peptide_bin::const_iterator peptide_bin::begin() const { return peptide_names.begin(); }
peptide_bin::const_iterator peptide_bin::end() const { return peptide_names.end(); }

peptide_bin::peptide_bin() = default;

peptide_bin::size_type peptide_bin::size() const
{
    return peptide_names.size();
}

bool peptide_bin::operator==( const peptide_bin& other ) const
{
    return peptide_names == other.peptide_names;
}

void peptide_bin::add_peptide( value_type to_add )
{
    peptide_names.emplace( to_add );
}

bool peptide_bin::contains( const value_type& look_for )
{
    return peptide_names.find( look_for ) != peptide_names.end();
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

bin_collection::size_type bin_collection::size() const
{
    return bins.size();
}

peptide_bin bin_collection::smallest()
{
    peptide_bin smallest = *bins.begin();
    for( const auto& bin : bins )
        {
            if( smallest.size() > bin.size() )
                {
                    smallest = bin;
                }
        }
    return smallest;
}

// namespace peptide_bin_io
bin_collection peptide_bin_io::parse_bins( std::istream& bin_source )
{
    std::string line;
    std::vector<std::string> peptides_in_bin;
    bin_collection bins;

    while( std::getline( bin_source, line ) )
        {

            boost::algorithm::trim( line );

            boost::split( peptides_in_bin,
                          line,
                          boost::is_any_of( "\t" )
                        );
            if( !( line.empty() || peptides_in_bin.empty() ) )
                {

                    peptide_bin new_bin( peptides_in_bin.begin(),
                                         peptides_in_bin.end()
                                         );
                    bins.add_bin( new_bin );
                }

        }

    return bins;
}

void peptide_bin_io::write_bins( std::ostream &dest, const bin_collection &bins )
{
    for( auto bin = bins.begin();
         bin != bins.end();
         ++bin
       )
        {
            dest << boost::algorithm::join( *bin, "\t" );
            dest << "\n";
        }
}
