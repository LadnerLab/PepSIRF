#ifndef DISTANCE_MATRIX_HH_INCLUDED
#define DISTANCE_MATRIX_HH_INCLUDED
#include <vector>

/**
 * Struct used for tracking location access
 * of the distance matrix.
 **/
struct loc
{
    /**
     * The row id to select
     **/
    int row_id;

    /**
     * The column id to select.
     **/
    int column_id;
};

template<typename T>
class distance_matrix
{
    
 public:
    distance_matrix() = default;

    typedef T value_type;

    distance_matrix( std::size_t n ){ reserve( n ); }

    std::vector<T>& operator[]( const std::size_t idx )
    {
        return data[ idx ];
    }

    T& operator[]( loc const& cLoc )
    {
        return data[ cLoc.row_id ][ cLoc.column_id ];
    }

    typename std::vector<std::vector<value_type>>::iterator
        const begin()
        {
            return data.begin();
        }

    typename std::vector<std::vector<value_type>>::iterator
        const end()
        {
            return data.end();
        }

    std::size_t size()
        {
            return data.size();
        }

    void reserve( std::size_t n )
    {
        std::size_t index = 0;
        for( index = 0; index < n; ++index )
            {
                data.emplace_back( std::vector<value_type>() );
                data[ index ].reserve( n );
            }
    }


 private:
    std::vector<std::vector<T>> data;
};

#endif // DISTANCE_MATRIX_HH_INCLUDED
