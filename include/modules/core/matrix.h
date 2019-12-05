#ifndef MATRIX_HH_INCLUDED
#define MATRIX_HH_INCLUDED
#include <cstdint>
#include <new>
#include <stdexcept>

#ifdef MATRIX_CHECK_BOUNDS

#define ACCESS_MATRIX( x, y ) \
({                        \
 if( ( ( (x) * N ) + (y) ) >= ( N * M ) ) \
     { \
       throw std::out_of_range( "Matrix access too large" ); \
     } \
    arr[ access_to_1d( (x), (y) ) ]; \
})

#else
#define ACCESS_MATRIX( x, y ) \
({                        \
    arr[ access_to_1d( (x), (y) ) ];  \
})

#endif

/**
 * A simple matrix to implement common matrix 
 * operations.
 **/
template <typename ValType>
class matrix
{
    /**
     * Construct an N x M matrix
     * @param in_N the x-dimension of the matrix to be created.
     *        (the number of rows)
     * @param in_M the y-dimension of the matriy to be created.
     *        (the number of columns)
     * @throws std::bad_alloc if allocating data for the matrix fails
     **/
    matrix( const std::uint32_t in_N, const std::uint32_t in_M )
        : N( in_N ), M( in_M )
        {
            arr = (ValType*) malloc( sizeof( ValType ) * N * M );

            if( !arr )
                {
                    throw std::bad_alloc();
                }
        }

    /**
     * Return the allocated memory
     **/
    ~matrix()
        {
            delete arr;
        }


 private:
    std::uint32_t N;
    std::uint32_t M;

    ValType *arr;

    std::uint32_t access_to_1d( const std::uint32_t x, const std::uint32_t y )
        {
            return ( x * N ) + y;
        }
    
};

#endif // MATRIX_HH_INCLUDED
