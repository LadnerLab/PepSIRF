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

    std::vector<T>& operator[]( const std::size_t idx )
    {
        return data[ idx ];
    }

    T& operator[]( loc const& cLoc )
    {
        return data[ cLoc.row_id ][ cLoc.column_id ];
    }

    std::size_t size()
        {
            return data.size();
        }


 private:
    std::vector<std::vector<T>> data;
};

#endif // DISTANCE_MATRIX_HH_INCLUDED
