#include "peptide.h"

const std::string& peptide::get_name() { return name; }
const std::string& peptide::get_sequence() { return sequence; }

void peptide::set_name( const std::string& new_val )
{
    name = new_val;
}

void peptide::set_sequence( const std::string& new_val )
{
    sequence = new_val;
}
