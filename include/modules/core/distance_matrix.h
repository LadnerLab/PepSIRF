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

    typedef T value_type;

    distance_matrix()
        {
            data = new std::vector<std::vector<value_type>*>();
        }


    distance_matrix( std::size_t n )
        {
            data = new std::vector<std::vector<value_type>*>();
            reserve( n );
        }

    std::vector<T>& operator[]( const std::size_t idx )
    {
        assert( idx < N );
        return *data->at( idx );
    }

    T& operator[]( loc const& cLoc )
    {
        return *(data[ cLoc.row_id ])[ cLoc.column_id ];
    }

    typename std::vector<std::vector<value_type>>::iterator
        const begin()
        {
            return data->begin();
        }

    typename std::vector<std::vector<value_type>>::iterator
        const end()
        {
            return data->end();
        }

    std::size_t size()
        {
            return data->size();
        }

    void reserve( std::size_t n )
    {
        N = n;
        std::size_t index = 0;
        for( index = 0; index < n; ++index )
            {
                data->push_back( new std::vector<value_type>() );
                data->at( index )->reserve( n );
            }
    }

    ~distance_matrix()
        {
            for( std::size_t index = 0; index < N; ++index )
                {
                    delete (*data)[ index ];
                }
            delete data;
        }


 private:
    std::vector<std::vector<T>*> *data;
    std::size_t N;
};

#endif // DISTANCE_MATRIX_HH_INCLUDED
