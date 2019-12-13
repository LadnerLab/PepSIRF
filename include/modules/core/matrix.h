#ifndef MATRIX_HH_INCLUDED
#define MATRIX_HH_INCLUDED
#include <cstdint>
#include <new>
#include <stdexcept>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <iostream>

#ifdef MATRIX_CHECK_BOUNDS

#define ACCESS_MATRIX( x, y ) \
({                        \
    if( !this->in_range( (x), (y) )) \
     { \
       throw std::out_of_range( "Matrix access too large" ); \
     } \
    return this->arr[ this->access_to_1d( (x), (y) ) ]; \
})

#else
#define ACCESS_MATRIX( x, y ) \
({                        \
    return this->arr[ this->access_to_1d( (x), (y) ) ]; \
})

#endif

/**
 * A simple matrix to implement common matrix 
 * operations.
 **/
template <typename ValType>
class matrix
{
 public:
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
     * Set all of the members of the matrix to 
     * a certain value.
     * @param val The value to set all members of the matrix to.
     **/
    void set_all( const ValType& val)
    {
        std::uint32_t idx = 0;

        for( idx = 0; idx < ( access_to_1d( N, M ) ); ++idx )
            {
                arr[ idx ] = val;
            }
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
     * Get the shape of the matrix, where 
     * the shape of an N x M matrix is (N,M).
     * @returns a pair containing the number of rows 
     *          and number of columns in the matrix.
     **/
    const std::pair<std::uint32_t,std::uint32_t>
        get_shape() const
        {
            return std::make_pair( N, M );
        };

    /**
     * Return the allocated memory
     **/
    ~matrix()
        {
            delete arr;
        }


 protected:
    std::uint32_t N;
    std::uint32_t M;

    ValType *arr;

    /**
     * Turn an (x,y) coordinate into a 1-dimensional coordinate.
     **/
    std::uint32_t access_to_1d( const std::uint32_t x, const std::uint32_t y )
        {
            return ( x * M ) + y;
        }

    /**
     * Check if accessing the (x,y)th element of this matrix will 
     * be valid.
     **/    
    bool in_range( const std::uint32_t x, const std::uint32_t y )
    {
        return access_to_1d( x, y ) < ( N * M );
    }

    /**
     * Returns a pointer to the first item in this matrix's
     * data. For internal use only.
     * @note Does not perform any bounds or safety checks.
     * @param row_idx The 0-based row number to return.
     * @returns const pointer to the first item in the row_idx'th row
     **/
    const ValType *get_row_ptr( std::uint32_t row_idx ) const
    {
        return &(this->arr[ row_idx * this->M ]);
    }
    
};


/**
 * A matrix whose rows and columns are labeled. 
 * A label can be any type, but is likely a string. 
 * @tparam ValType The type of the values stored in the matrix.
 * @tparam LabelType The type of the labels.
 **/
template <typename ValType,
          typename LabelType
        >
class labeled_matrix : public matrix<ValType>
{
 public:
    
    /**
     * Constructor that initializes the matrix itself without any labels.
     * Initializes in_N \times in_M values of type ValType
     * @param in_N The number of rows the matrix has.
     * @param in_M The number of columns the matrix has.
     **/
    labeled_matrix( const std::uint32_t in_N, const std::uint32_t in_M )
        : matrix<ValType>{ in_N, in_M } {}

    /**
     * Initialize a labeled matrix with row and column labels.
     * @tparam RowLabelContainerType The container that holds the values for 
     *         the row labels. This container must be ordered, containing 
     *         begin() and end() functions such that for LabelContainerType L,
     *         L.begin() + N is the label for the N'th row.
     * @tparam ColLabelContainerType The container that holds the values for 
     *         the column labels. This container must be ordered, containing 
     *         begin() and end() functions such that for LabelContainerType L,
     *         L.begin() + N is the label for the N'th column.
     * @param in_N The number of rows in the matrix.
     * @param in_M The number of columns in the matrix.
     * @param row_labels The labels for each row. If any row has a label,
     *        then each row should have a label, so row_labels.size() = in_N.
     * @param col_labels The labels for each column. If any column has a label,
     *        then each column should have a label, so col_labels.size() = in_M.
     **/
    template<typename RowLabelContainerType,
             typename ColLabelContainerType
           >
        labeled_matrix( const std::uint32_t in_N, const std::uint32_t in_M,
                         const RowLabelContainerType& row_labels,
                         const ColLabelContainerType& col_labels
                       )
        : matrix<ValType>{ in_N, in_M }
    {
        init_labels( row_labels,
                     this->row_labels
                   );

        init_labels( col_labels,
                     this->col_labels
                   );

    }

    using matrix<ValType>::operator();
    using matrix<ValType>::at;

    /**
     * Access a mutable member of the matrix using its row/column labels instead 
     * of integer indices. 
     * @note Performs validity check on the labels and bounds checks on the access
     * @param row_lab The label of the row to access
     * @param col_lab The label of the column to access
     * @throws std::out_of_range if row_lab and col_lab are not valid row/column labels.
     * @throws std::out_of_range if (row_lab, col_lab) resolves to indices 
     *         that are not in range of the matrix.
     * @returns mutable reference to the matrix at (row_lab, col_lab) 
     **/
    ValType &at( const LabelType& row_lab, const LabelType& col_lab )
        {
            const auto row_val = row_labels.find( row_lab );
            const auto col_val = col_labels.find( col_lab );

            if( row_val == row_labels.end() )
                {
                    throw std::out_of_range( "Bad row access label! Item not found." );
                }
            else if( col_val == col_labels.end() )
                {
                    throw std::out_of_range( "Bad column access label! Item not found." );
                }

            return this->at( row_val->second, col_val->second );
        }

    /**
     * Access a constant member of the matrix using its row/column labels instead 
     * of integer indices. 
     * @note Performs validity check on the labels and bounds checks on the access
     * @param row_lab The label of the row to access
     * @param col_lab The label of the column to access
     * @throws std::out_of_range if row_lab and col_lab are not valid row/column labels.
     * @throws std::out_of_range if (row_lab, col_lab) resolves to indices 
     *         that are not in range of the matrix.
     * @returns constant reference to the matrix at (row_lab, col_lab) 
     **/
    const ValType &at( const LabelType& row_lab, const LabelType& col_lab ) const
        {
            const auto& row_val = row_labels[ row_lab ];
            const auto& col_val = col_labels[ col_lab ];

            if( row_val == row_labels.end() )
                {
                    throw std::out_of_range( "Bad row access label! Item not found." );
                }
            else if( col_val == col_labels.end() )
                {
                    throw std::out_of_range( "Bad column access label! Item not found." );
                }

            return this->at( row_val->second, col_val->second );
        }

    /**
     * Access a mutable member of the matrix using its row/column labels instead of 
     * integer indices.
     * @note Does NOT perform bounds check on the access, behavior is undefined if 
     *       the access is invalid.
     * @returns Mutable reference to the (row_lab, col_lab) member of the matrix.
     **/
    ValType &operator()( const LabelType& row_lab, const LabelType& col_lab  )
        {
            const std::uint32_t row_idx = row_labels[ row_lab ];
            const std::uint32_t col_idx = col_labels[ col_lab ];

            ACCESS_MATRIX( row_idx, col_idx );
        }

    /**
     * Access a constant member of the matrix using its row/column labels instead of 
     * integer indices.
     * @note Does NOT perform bounds check on the access, behavior is undefined if 
     *       the access is invalid.
     * @returns Constant reference to the (row_lab, col_lab) member of the matrix.
     **/
    const ValType &operator()( const LabelType& row_lab, const LabelType& col_lab  ) const
        {
            std::uint32_t row_idx = row_labels[ row_lab ];
            std::uint32_t col_idx = col_labels[ col_lab ];

            ACCESS_MATRIX( row_idx, col_idx );
        }

    /**
     * Initialize the locally-stored labels for the matrix.
     * Stores each label with its corresponding row/column number in label_dest
     * @tparam LabelContainerType the type of the container holding labels.
     * @param label_dest The place to store the label/index mappings.
     * @post label_dest.size() == labels.size()
     * @post For label L and index J, label_dest[ L ] = J
     **/
    template<typename LabelContainerType>
    void init_labels( const LabelContainerType &labels,
                      std::unordered_map<LabelType, std::uint32_t>& label_dest
                    )
        {
            std::uint32_t idx = 0;
            label_dest.reserve( std::distance( labels.begin(), labels.end() ) );

            for( const auto& iter_loc : labels )
                {
                    label_dest.emplace( iter_loc, idx++ );
                }
        }


 private:
    /**
     * The labels for each row of the matrix.
     **/
    std::unordered_map<LabelType, std::uint32_t>
        row_labels;
    /**
     * The labels for each column of the matrix.
     **/
    std::unordered_map<LabelType, std::uint32_t>
        col_labels;

};

#endif // MATRIX_HH_INCLUDED
