#include "sequence.h"

sequence::sequence()
{
    name = std::string();
    seq  = std::string();
}
sequence::sequence( const std::string& in_name, const std::string& in_seq )
{
    name = in_name;
    seq  = in_seq;
}

sequence::sequence( const sequence& s )
{
    name = s.name;
    seq = s.seq;
}

sequence::~sequence() = default;

void sequence::set_name( std::string& new_name )
{
    name = new_name;
}

void sequence::set_seq( std::string& new_seq )
{
    seq = new_seq;
}

int sequence::count( const char val )
{
    unsigned int index = 0;
    int count = 0;

    for( index = 0; index < seq.length(); ++index )
        {
            if( seq[ index ] == val )
                {
                    ++count;
                }
        }
    return count;
}

std::size_t sequence::length()
{
    return seq.length();
}

bool sequence::operator==( const sequence& s ) const
{
    return !s.seq.compare( seq );
}

bool sequence::operator<( const sequence& s ) const
{
    return seq < s.seq;
}
