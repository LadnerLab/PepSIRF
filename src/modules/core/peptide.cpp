#include "peptide.h"

const std::string& peptide::get_name() const
{
    return name;
}

const std::string& peptide::get_sequence() const
{
    return sequence;
}

void peptide::set_name( const std::string& new_val )
{
    name = new_val;
}

void peptide::set_sequence( const std::string& new_val )
{
    sequence = new_val;
}


bool peptide::operator==( const peptide& other ) const
{
    return !get_sequence().compare( other.get_sequence() );
}

bool peptide::operator!=( const peptide& other ) const
{
    return !( *this == other );
}

