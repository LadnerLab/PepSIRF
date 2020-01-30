#ifndef NOOP_ITERATOR_HH_INCLUDED
#define NOOP_ITERATOR_HH_INCLUDED
#include <iterator>

/**
 * No operation (noop) iterator.
 * Assignment, incrementation, etc. all do nothing.
 * Useful if you don't want to use the result assigned to 
 * an output iterator.
 **/
struct noop_iterator
    : std::iterator<std::output_iterator_tag,
                    noop_iterator
                   >
{

    template<typename T>
    void operator=( const T& ) {}
    noop_iterator &operator++() {return *this; }
    noop_iterator &operator++( int ) {return *this; }
    noop_iterator &operator*(){ return *this; };

};


#endif 
