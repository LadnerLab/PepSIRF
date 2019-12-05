#ifndef MATRIX_HH_INCLUDED
#define MATRIX_HH_INCLUDED
#include <cstdint>
#include <new>
#include <stdexcept>

#ifdef MATRIX_CHECK_BOUNDS

#define ACCESS_MATRIX( x, y ) \
({                        \
    if( !in_range( (x), (y) ))                  \
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
     * Access the constant (x,y) element of the array, checking that we will be 
     * performing a check to ensure the access is in bounds.
     * @param x The index of the row to acccess
     * @param y the index of the column to access.
     * @returns constant reference to the (x,y) item in the matrix.
     * @throws std::out_of_range if (x,y) is greater than the shape of the 
     *         matrix.
     **/
    const ValType& at( const std::uint32_t x, const std::uint32_t y ) const
    {
        if( !in_range( x, y ) )
            {
                throw std::out_of_range( "Matrix access too large" );
            }
        return arr[ access_to_1d( x, y ) ] ;
    }


    /**
     * Access the mutable (x,y) element of the array, checking that we will be 
     * performing a check to ensure the access is in bounds.
     * @param x The index of the row to acccess
     * @param y the index of the column to access.
     * @returns mutable reference to the (x,y) item in the matrix.
     * @throws std::out_of_range if (x,y) is greater than the shape of the 
     *         matrix.
     **/
    ValType& at( const std::uint32_t x, const std::uint32_t y )
    {
        if( !in_range( x, y ) )
            {
                throw std::out_of_range( "Matrix access too large" );
            }
        return arr[ access_to_1d( x, y ) ] ;
    }

    /**
     * Access the mutable (x,y) element of the array.
     * @param x The index of the row to acccess
     * @param y the index of the column to access.
     * @returns mutable reference to the (x,y) item in the matrix.
     * @note if the compile-time constant 'MATRIX_CHECK_BOUNDS' is defined, 
     *       this method throws std::out_of_range on invalid access.
     **/
    ValType& operator()( const std::uint32_t x, const std::uint32_t y )
        {
            ACCESS_MATRIX( x, y );
        }

    /**
     * Access the constant (x,y) element of the array.
     * @param x The index of the row to acccess
     * @param y the index of the column to access.
     * @returns constant reference to the (x,y) item in the matrix.
     * @note if the compile-time constant 'MATRIX_CHECK_BOUNDS' is defined, 
     *       this method throws std::out_of_range on invalid access.
     **/
    const ValType& operator()( const std::uint32_t x, const std::uint32_t y ) const
        {
            ACCESS_MATRIX( x, y );
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

    /**
     * Turn an (x,y) coordinate into a 1-dimensional coordinate.
     **/
    std::uint32_t access_to_1d( const std::uint32_t x, const std::uint32_t y )
        {
            return ( x * N ) + y;
        }

    /**
     * Check if accessing the (x,y)th element of this matrix will 
     * be valid.
     **/    
    bool in_range( const std::uint32_t x, const std::uint32_t y )
    {
        return access_to_1d( x, y ) < ( N * M );
    }
    
};

#endif // MATRIX_HH_INCLUDED
