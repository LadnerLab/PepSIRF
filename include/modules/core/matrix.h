#ifndef MATRIX_HH_INCLUDED
#define MATRIX_HH_INCLUDED
#include <cstdint>
#include <new>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include <iostream>
#include <boost/algorithm/string/join.hpp>

#include "setops.h"

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
    // TODO: Constant/mutable iterators

    /**
     * An iterator providing functionality for iterating over rows/columns
     * of a matrix.
     **/
    class mutable_iterator
    {

    public:
        /**
         * Initialize the mutable_iterator with a matrix,
         * sets mutable_iterator pointer to zero.
         **/
    mutable_iterator(  ::matrix<ValType> *matr )
            : matr( matr ),
              current_idx( 0 )
              {}

        /**
         * Initialize with a matrix and starting (row, column) position.
         * @param matr The matrix to iterate over.
         * @param start_row The row to begin iteration on
         * @param start_col the column to begin iteration on
         **/
    mutable_iterator(  ::matrix<ValType> *matr,
              std::size_t start_row,
              std::size_t start_col
                       )
            : matr( matr ),
              current_idx( matr->access_to_1d( start_row, start_col ) )
              {}


        /**
         * Determines whether the mutable_iterator has reached
         * at least the end of the matrix.
         * @returns boolean True if the current pointer is greater than
         *          the size of the matrix, False otherwise.
         **/
        bool end()
        {
            return current_idx >=
                 this->matr.access_to_1d( this->matr.N, this->matr.M );
        }

        /**
         * ant access to the current value pointed to by the
         * mutable_iterator.
         * @returns ant reference to the current iterated value
         **/
         ValType& operator*()
        {
            return this->matr->operator()( this->current_idx );
        }

        /**
         * In-place advance the mutable_iterator forward.
         * @param advance The number of positions to
         *        advance.
         **/
        void operator+=( int advance )
        {
            current_idx += advance;
        }

        /**
         * Increment in-place, return reference to
         * this mutable_iterator.
         **/
        mutable_iterator &operator++()
            {
                *this += 1; return *this;
            }

        /**
         * Advance forward an entire row from the
         * current position.
         * @note This does not advance to the beginning of
         *       the row, but advances to the current position in
         *       the next row.
         **/
        mutable_iterator& next_row()
            {
                this += this->matr.M;
                return *this;
            }

        /**
         * Advance forward an entire column from the
         * current position.
         * @note This does not advance to the beginning of
         *       the column, but advances to the current position in
         *       the next column.
         **/
        mutable_iterator &next_col()
            {
                this += this->matr.N;
                return *this;
            }

        /**
         * Determines whether this mutable_iterator is NOT equal to
         * another.
         * @note We say two mutable_iterators are equal if they point to the
         *       same position in the same matrix. See matrix::mutable_iterator::operator==
         * @returns True if !( this == comp ), false otherwise
         **/
        bool operator!=(  const mutable_iterator& comp ) const
        {
            return !( *this == comp );
        }

        /**
         * Determine whether two mutable_iterators are equal.
         * @note We define two mutable_iterators as equivalent if
         *       they point to the same position in the same
         *       matrix.
         * @returns True if this == comp, false otherwise.
         **/
        bool operator==( const mutable_iterator& comp ) const
        {
            return this->matr == comp.matr
                     && this->current_idx == comp.current_idx;
        }

        /**
         * Assign this mutable_iterator to another.
         * @param other The mutable_iterator to copy.
         * @note After calling this method, this == other will be true,
         *       and it is thus important to note that this method does not
         *       copy the matrix pointed to by other, it copies only its
         *       pointer.
         **/
        void operator=( const typename matrix<ValType>::mutable_iterator& other )
            {
                this->matr = other.matr;
                this->current_idx = other.current_idx;
            }

        using difference_type = std::ptrdiff_t;
        using iterator_category	= std::forward_iterator_tag;
        using value_type = ValType;
        using pointer = ValType*;
        using reference = ValType&;

    private:
        /**
         * The matrix being iterated over.
         **/
         ::matrix<ValType> *matr;

        /**
         * The current location of the mutable_iterator.
         **/
        std::size_t current_idx;


    };


    class iterator
    {

    public:
        /**
         * Initialize the iterator with a matrix,
         * sets iterator pointer to zero.
         **/
    iterator( const ::matrix<ValType> *matr )
            : matr( matr ),
              current_idx( 0 )
              {}

        /**
         * Initialize with a matrix and starting (row, column) position.
         * @param matr The matrix to iterate over.
         * @param start_row The row to begin iteration on
         * @param start_col the column to begin iteration on
         **/
    iterator( const ::matrix<ValType> *matr,
              std::size_t start_row,
              std::size_t start_col
                       )
            : matr( matr ),
              current_idx( matr->access_to_1d( start_row, start_col ) )
              {}


        /**
         * Determines whether the iterator has reached
         * at least the end of the matrix.
         * @returns boolean True if the current pointer is greater than
         *          the size of the matrix, False otherwise.
         **/
        bool end()
        {
            return current_idx >=
                 this->matr.access_to_1d( this->matr.N, this->matr.M );
        }

        /**
         * Constant access to the current value pointed to by the
         * iterator.
         * @returns constant reference to the current iterated value
         **/
        const ValType& operator*() const
        {
            return this->matr->operator()( this->current_idx );
        }

        /**
         * In-place advance the iterator forward.
         * @param advance The number of positions to
         *        advance.
         **/
        void operator+=( int advance )
        {
            current_idx += advance;
        }

        /**
         * Increment in-place, return reference to
         * this iterator.
         **/
        iterator &operator++()
            {
                *this += 1; return *this;
            }

        /**
         * Advance forward an entire row from the
         * current position.
         * @note This does not advance to the beginning of
         *       the row, but advances to the current position in
         *       the next row.
         **/
        iterator& next_row()
            {
                this += this->matr.M;
                return *this;
            }

        /**
         * Advance forward an entire column from the
         * current position.
         * @note This does not advance to the beginning of
         *       the column, but advances to the current position in
         *       the next column.
         **/
        iterator &next_col()
            {
                this += this->matr.N;
                return *this;
            }

        /**
         * Determines whether this iterator is NOT equal to
         * another.
         * @note We say two iterators are equal if they point to the
         *       same position in the same matrix. See matrix::iterator::operator==
         * @returns True if !( this == comp ), false otherwise
         **/
        bool operator!=( const iterator& comp ) const
        {
            return !( *this == comp );
        }

        /**
         * Determine whether two iterators are equal.
         * @note We define two iterators as equivalent if
         *       they point to the same position in the same
         *       matrix.
         * @returns True if this == comp, false otherwise.
         **/
        bool operator==( const iterator& comp ) const
        {
            return this->matr == comp.matr
                     && this->current_idx == comp.current_idx;
        }

        /**
         * Assign this iterator to another.
         * @param other The iterator to copy.
         * @note After calling this method, this == other will be true,
         *       and it is thus important to note that this method does not
         *       copy the matrix pointed to by other, it copies only its
         *       pointer.
         **/
        void operator=( const typename matrix<ValType>::iterator& other )
            {
                this->matr = other.matr;
                this->current_idx = other.current_idx;
            }

    private:
        /**
         * The matrix being iterated over.
         **/
        const ::matrix<ValType> *matr;

        /**
         * The current location of the iterator.
         **/
        std::size_t current_idx;

    };

    /**
     * Returns an iterator to
     * this matrix whose pointer is at the
     * first element in the matrix.
     **/
    iterator get_iterator() const
    {
        return iterator( this );
    }

    iterator get_end_iterator() const
    {
        return std::end( this );
    }

    /**
     * Get an iterator that is pointing to the end of
     * the n'th row.
     * @note No safety checks are made to ensure
     *       the chosen row is actually in the matrix.
     * @param row The row to point to the end of.
     **/
    const iterator row_end( std::size_t row ) const
    {
        return iterator( this, row + 1, 0 );
    }

    /**
     * Get an iterator to the beginning of the n'th row.
     **/
    const iterator row_begin( std::size_t row ) const
    {
        return iterator( this, row, 0 );
    }


    const mutable_iterator row_end( std::size_t row )
    {
        return mutable_iterator( this, row + 1, 0 );
    }

    /**
     * Get an iterator to the beginning of the n'th row.
     **/
    const mutable_iterator row_begin( std::size_t row )
    {
        return mutable_iterator( this, row, 0 );
    }


    /**
     * Default constructor.
     **/
    matrix()
        : N{ 0 }, M{ 0 }, arr{ nullptr } {};

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
            arr = new ValType[ N * M ]();

            if( !arr )
                {
                    throw std::bad_alloc();
                }
        }

    matrix( const matrix<ValType>& other )
        : N( other.N ), M( other.M ), arr( new ValType[ N * M ] )
    {
        for( std::uint32_t idx = 0; idx < N * M; ++idx )
            {
                arr[ idx ] = other.arr[ idx ];
            }
    }

    /**
     * Move constructor, takes ownership of to_move's data.
     **/
    matrix( matrix<ValType>&& to_move )
        : matrix()
        {
            swap( *this, to_move );
        }

    /**
     * Swap this matrix with another.
     **/
    friend void swap( matrix<ValType>& first, matrix<ValType>& second )
    {
        using std::swap;

        swap( first.arr, second.arr );
        swap( first.M, second.M );
        swap( first.N, second.N );
    }

    /**
     * Assign this matrix to another.
     **/
    matrix<ValType>& operator=( matrix<ValType> other )
        {
            swap( *this,other );
            return *this;
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

        for( idx = 0; idx <  N * M ; ++idx )
            {
                arr[ idx ] = val;
            }
    }

        /**
         * Determine whether this matrix and other are equal.
         * For two matrices, A and B, we say A == B if A and B have
         * the same shape, and all of the pairwise entries of A are equal
         * to all of the pairwise entries of B.
         * @param other The matrix to compare this with.
         * @returns true if this and other are equal, false otherwise.
         **/
        bool operator==( const matrix<ValType>& other ) const
        {
            bool same = ( this->M == other.M ) && ( this->N == other.N );

            for( std::uint32_t row = 0; row < this->N && same; ++row )
                {
                    for( std::uint32_t col = 0; col < this->M && same; ++col )
                        {
                            same = this->operator()( row, col ) == other( row, col );
                        }
                }
            return same;
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
     * Access a constant member of the matrix using 1-dimensional
     * access.
     * @note no safety checks are made to ensure the access is valid.
     * @returns constant reference to the access_idx'th item in
     *          the matrix.
     **/
    const ValType& operator()( const std::uint32_t access_idx ) const
        {
            return this->arr[ access_idx ];
        }

    /**
     * Access a mutable member of the matrix using 1-dimensional
     * access.
     * @note no safety checks are made to ensure the access is valid.
     * @returns mutable reference to the access_idx'th item in
     *          the matrix.
     **/
    ValType& operator()( const std::uint32_t access_idx )
        {
            return this->arr[ access_idx ];
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
     * Get the number of rows in the matrix.
     **/
    std::size_t nrows() const
        {
            return N;
        }

    /**
     * Get the number of columns in the matrix.
     **/
    std::size_t ncols() const
        {
            return M;
        }

    /**
     * Get the total size (total number of possible elements)
     * of the matrix.
     *
     **/
    std::size_t total_size() const
        {
            return N * M;
        }


    /**
     * Compare the row_idx'th row of this matrix
     * with the row_idx'th row of comp.
     * For the rows we say they are equal if they have
     * the same number of columns AND each item in the row
     * is equal to the item in comp's row at that index.
     * @param comp The matrix to grab a row from to compare
     * @param row_idx The row to compare
     * @returns boolean true if the matrices have the same
     *          number of columns, and each item in the row
     *          of each are the same.
     **/
    bool compare_row( const matrix<ValType>& comp,
                      std::uint32_t row_idx
                    )
    {
        const ValType *row_start_self = this->get_row_ptr( row_idx );
        const ValType *row_start_comp = comp.get_row_ptr( row_idx );

        bool diff = this->M == comp.M;

        for( std::uint32_t col_idx = 0;
             col_idx < this->M && !diff;
             ++col_idx
           )
            {
                diff = *row_start_self == *row_start_comp;
                ++row_start_self;
                ++row_start_comp;
            }

        return diff;
    }

    /**
     * Return the allocated memory
     **/
    ~matrix()
        {
            delete[] arr;
        }


 protected:
    /**
     * The number of rows
     * the matrix has
     **/
    std::uint32_t N;

    /**
     * The number of columns the matrix has.
     **/
    std::uint32_t M;

    /**
     * The values stored in the matrix.
     **/
    ValType *arr;

    /**
     * Turn an (x,y) coordinate into a 1-dimensional coordinate.
     **/
    std::uint32_t access_to_1d( const std::uint32_t x, const std::uint32_t y ) const
        {
            return ( x * M ) + y;
        }

    /**
     * Turn a 1-dimensional coordinate in to an (x,y) coordinate.
     **/
    std::pair<std::uint32_t,std::uint32_t> access_to_2d( const std::uint32_t idx ) const
        {
            std::pair<std::uint32_t,std::uint32_t> transform_acc_vec;
            transform_acc_vec.first = idx / this->M;

            transform_acc_vec.second = idx / this->N;
            return transform_acc_vec;
        }

    /**
     * Check if accessing the (x,y)th element of this matrix will
     * be valid.
     **/
    bool in_range( const std::uint32_t x, const std::uint32_t y ) const
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

    labeled_matrix()
        : matrix<ValType>{ 0, 0 } {};
    /**
     * Constructor that initializes the matrix itself without any labels.
     * Initializes in_N \times in_M values of type ValType
     * @param in_N The number of rows the matrix has.
     * @param in_M The number of columns the matrix has.
     **/
    labeled_matrix( const std::uint32_t in_N, const std::uint32_t in_M )
        : matrix<ValType>{ in_N, in_M }, row_labels{}, col_labels{} {}

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

    /**
     * Construct a labeled_matrix given row and column labels.
     **/
 labeled_matrix( const std::uint32_t in_N, const std::uint32_t in_M,
                 const std::unordered_map<LabelType,std::uint32_t> row_labels,
                 const std::unordered_map<LabelType,std::uint32_t> col_labels
                 )
     : matrix<ValType>{ in_N, in_M }
    {
        this->row_labels = row_labels;
        this->col_labels = col_labels;
    }


    /**
     * Swap first and second.
     **/
    void friend swap( labeled_matrix<ValType,LabelType>& first,
                      labeled_matrix<ValType,LabelType>& second
                    )
    {
        using std::swap;
        swap( static_cast<matrix<ValType>&>( first ),
              static_cast<matrix<ValType>&>( second )
            );
        swap( first.row_labels,
              second.row_labels
            );
        swap( first.col_labels,
              second.col_labels
            );
    }

    /**
     * Move constructor.
     **/
    labeled_matrix( labeled_matrix<ValType,LabelType>&& to_move )
        : labeled_matrix()
        {
            swap( *this, to_move );
        }

    labeled_matrix( const labeled_matrix<ValType,LabelType>& copy )
        : matrix<ValType>{ static_cast<const matrix<ValType>>( copy ) },
          row_labels{ copy.row_labels },
          col_labels{ copy.col_labels }
    {}

    /**
     * Assignment operator that uses move-and-swap idiom.
     **/
    labeled_matrix<ValType,LabelType>& operator=( labeled_matrix<ValType,LabelType>
                                                  to_assign
                                                )
        {
            swap( *this, to_assign );
            return *this;
        }

    using matrix<ValType>::operator();
    using matrix<ValType>::at;

    /**
     * Set the row labels for this matrix.
     * @tparam ContainerType the container labels are currently stored in.
     *         this should be an ordered container.
     * @pre the values in labels are in the same order as the
     *      rows of this matrix.
     * @post This matrix's row labels are set to the labels in containers
     **/
    template<typename ContainerType>
        void set_row_labels( const ContainerType& labels )
        {
            init_labels( labels,
                         this->row_labels
                       );
        }
    /**
     * Update the row labels for this matrix.
     * @param labels A vector containing the replacement labels.
     * @post The vector labels must be in the same order as the original labels
     **/
    void update_row_labels( const std::vector<std::string> labels )
        {
            update_labels( labels,
                           this->row_labels );
        }

    /**
     * Update the col labels for this matrix.
     * @param labels A vector containing the replacement labels.
     * @post The vector labels must be in the same order as the original labels.
     *       If these are not in the same order, then the labels may be
     *       incorrectly assigned.
     **/
    void update_col_labels( const std::vector<std::string> labels )
        {
            update_labels( labels,
                           this->col_labels );
        }

    /**
     * Get the row labels of the matrix, store them in
     * dest.
     * @tparam ContainerType The type of container to store labels in.
     *         This type should have an 'insert( hint, val )' method.
     * @param dest The location to store the labels in.
     * @note This method does not force any ordering of the labels,
     *       do not expect the n'th item in dest to be the label of the
     *       n'th row in the matrix.
     **/
    template<typename ContainerType>
        ContainerType& get_row_labels( ContainerType& dest ) const
        {
            get_labels( dest,
                        this->row_labels
                      );
            return dest;
        }

    /**
     * Get the column labels of the matrix, store them in
     * dest.
     * @tparam ContainerType The type of container to store labels in.
     *         This type should have an 'insert( hint, val )' method.
     * @param dest The location to store the labels in.
     * @note This method does not force any ordering of the labels,
     *       do not expect the n'th item in dest to be the label of the
     *       n'th column in the matrix.
     **/
    template<typename ContainerType>
        ContainerType& get_col_labels( ContainerType& dest ) const
        {
            get_labels( dest,
                        this->col_labels
                        );

            return dest;
        }

    /**
     * Get the labels from a map, store them in the container dest.
     **/
    template<typename ContainerType>
        void get_labels( ContainerType& dest,
                         const std::unordered_map<LabelType,std::uint32_t>& src
                         ) const
        {
            for( const auto lab : src )
                {
                    dest.insert( dest.end(), lab.first );
                }
        }

    const std::unordered_map<LabelType,std::uint32_t>&
        get_row_labels()
        {
            return row_labels;
        }

    const std::unordered_map<LabelType,std::uint32_t>&
        get_col_labels()
        {
            return col_labels;
        }

    /**
     * Full outer join this matrix with other.
     * Creates a matrix with all of the rows and columns of
     * this and other. Any values that are not in the join will
     * be zero-initialized, as defined by ValType.
     * @param other The matrix to join with this.
     * @returns a labeled matrix with the rows and columns of this
     *          and other.
     **/
    labeled_matrix<ValType, LabelType>
        full_outer_join( const labeled_matrix<ValType, LabelType>& other )
    {
        // find a_columns union b_cols -> number of columns in output
        std::unordered_set<LabelType> col_lab_union;
        setops::set_union( col_lab_union,
                           this->col_labels,
                           other.col_labels,
                           setops
                           ::get_key<
                           LabelType,std::uint32_t
                           >()
                         );

        // find a_rows union b_cols -> number of rows in output
        std::unordered_set<LabelType> row_lab_union;
        setops::set_union( row_lab_union,
                           this->row_labels,
                           other.row_labels,
                           setops
                           ::get_key<
                           LabelType,std::uint32_t
                           >()
                         );

        std::uint32_t join_N = row_lab_union.size(),
                      join_M = col_lab_union.size();

        // create the output with the correct number of rows, columns
        labeled_matrix<ValType,LabelType> joined_matr( join_N, join_M,
                                                       row_lab_union,
                                                       col_lab_union
                                                     );

        for( const auto& row_lab : this->row_labels )
            {
                for( const auto& col_lab : this->col_labels )
                    {
                                joined_matr( row_lab.first, col_lab.first ) =
                                    this->operator()( row_lab.first, col_lab.first );
                    }
            }

        for( const auto& row_lab : other.row_labels )
            {
                for( const auto& col_lab : other.col_labels )
                    {
                                joined_matr( row_lab.first, col_lab.first ) =
                                    other( row_lab.first, col_lab.first );
                    }
            }
        return joined_matr;
    }

    typename matrix<ValType>::iterator find( const LabelType& row_lab,
                                             const LabelType& col_lab
                                           )
    {
        const auto& row_label = this->row_labels.find( row_lab );
        const auto& col_label = this->col_labels.find( col_lab );
        if( row_label == this->row_labels.end()
            || col_label == this->col_labels.end()
          )
            {
                return this->end();
            }

        return matrix<ValType>::iterator( this->arr, row_label.second, col_label.second );
    }

    /**
     * Set the column labels for this matrix.
     * @tparam ContainerType the container labels are currently stored in.
     *         this should be an ordered container.
     * @pre the values in labels are in the same order as the
     *      columns of this matrix.
     * @post This matrix's column labels are set to the labels in containers
     **/
    template<typename ContainerType>
        void set_col_labels( const ContainerType& labels )
        {
            init_labels( labels,
                         this->col_labels
                       );
        }

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
            std::uint32_t row_idx = row_labels.find( row_lab )->second;
            std::uint32_t col_idx = col_labels.find( col_lab )->second;

            ACCESS_MATRIX( row_idx, col_idx );
        }

    /**
     * Access an item in the matrix, given a label for the row and the index
     * of the column.
     * @param row_lab The label of the row to access
     * @param col_idx the index of the column to access
     **/
    ValType &operator()( const LabelType& row_lab, std::uint32_t col_idx )
        {
            std::uint32_t row_idx = row_labels.find( row_lab )->second;

            ACCESS_MATRIX( row_idx, col_idx );
        }

    /**
     * Access an item in the matrix, given an index for the row
     * and a label for the column.
     * @param row_idx The index of the row to grab
     * @param col_lab the label of the column to access
     **/
    ValType &operator()( std::uint32_t row_idx, const LabelType& col_lab )
        {
            std::uint32_t col_idx = col_labels.find( col_lab )->second;
            ACCESS_MATRIX( row_idx, col_idx );
        }

    void set_row_label( const LabelType& orig_label,
                        const LabelType& new_label
                      )
    {
        set_label( this->row_labels,
                   orig_label,
                   new_label
                 );
    }

    void set_col_label( const LabelType& orig_label,
                        const LabelType& new_label
                      )
    {
        set_label( this->col_labels,
                   orig_label,
                   new_label
                 );
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
    /**
     * Clears the locally-stored labels for the matrix and emplaces each label from labels.
     * Each label is stored in label_dest with the value corresponding to the labels index.
     * @param labels An ordered list of labels to be stored in matrix.
     * @param label_dest The place to store the label/index mappings.
     * @post labels is set as vector because an ordered list is necessary to associate index to label.
     **/
    void update_labels( const std::vector<std::string> &labels,
                        std::unordered_map<LabelType, std::uint32_t>& label_dest )
        {
            std::uint32_t idx = 0;
            label_dest.clear();
            label_dest.reserve( std::distance( labels.begin(), labels.end() ) );

            for( const auto& iter_loc : labels )
                {
                    label_dest.emplace( iter_loc, idx++ );
                }
        }

    /**
     * Filter the rows in this matrix, keeping only the rows whose
     * labels are found in the provided container.
     * @tparam ContainerType the type of container the labels are
     *         stored in. Must be iterable by a range-based for loop.
     * @param row_labels
     * @note row_labels the order of the labels in row_labels is
     *       preserved.
     * @returns An K x M matrix, where K is the size of row_labels. The returned
     *          matrix will contain the rows in this matrix whose labels are found in
     *          row labels. The columns and their labels will remain unchanged.
     **/
    template<typename ContainerType>
    labeled_matrix<ValType,LabelType>
        filter_rows( const ContainerType& row_labels )
        {
            std::vector<LabelType> column_labels;
            column_labels.resize( this->col_labels.size() );

            // assign the column labels to the a standalone data structure,
            // taking care to remember std::unordered_map is unordered
            for( const auto& col_data : this->col_labels )
                {
                    column_labels[ col_data.second ] = col_data.first;
                }

            std::uint32_t num_rows = row_labels.size();
            labeled_matrix<ValType,LabelType> return_matrix( num_rows, this->M,
                                                             row_labels,
                                                             column_labels
                                                           );

            std::uint32_t row_idx_new = 0;

            for( const auto& row_lab : row_labels )
                {
                    std::uint32_t row_idx_original = this->row_labels.at( row_lab );
                    const ValType *row_ptr = this->get_row_ptr( row_idx_original );

                    for( std::uint32_t row_val = 0; row_val < this->M; ++row_val )
                        {
                            return_matrix( row_idx_new, row_val ) = *row_ptr;
                            ++row_ptr;
                        }
                    ++row_idx_new;
                }
            return return_matrix;
        }

    /**
     * Filter the cols in this matrix, keeping only the rows whose
     * labels are found in the provided container.
     * @tparam ContainerType the type of container the labels are
     *         stored in. Must be iterable by a range-based for loop.
     * @param col_labels The order of the labels in col_labels is
     *        preserved.
     * @returns An N x K matrix, where K is the size of col_labels. The returned
     *          matrix will contain the cols in this matrix whose labels are found in
     *          col labels. The rows and their labels will remain unchanged.
     **/
    template<typename ContainerType>
    labeled_matrix<ValType,LabelType>
        filter_cols( const ContainerType& col_labels )
        {
            std::vector<LabelType> row_labels;
            row_labels.resize( this->row_labels.size() );

            // assign the row labels to the a standalone data structure,
            // taking care to remember std::unordered_map is unordered
            for( const auto& row_data : this->row_labels )
                {
                    row_labels[ row_data.second ] = row_data.first;
                }
            // resizing cols to only include filter results
            std::uint32_t num_cols = col_labels.size();
            labeled_matrix<ValType,LabelType> return_matrix( this->N, num_cols,
                                                             row_labels,
                                                             col_labels
                                                           );

            std::uint32_t col_idx_new = 0;

            for( const auto& col_lab : col_labels )
                {
                    std::uint32_t col_idx_original = this->col_labels.at( col_lab );
                    const ValType *col_ptr = &(this->arr[col_idx_original]);

                    for( std::uint32_t col_val = 0; col_val < this->N; ++col_val )
                        {
                            return_matrix( col_val, col_idx_new ) = *col_ptr;
                            ++col_ptr;
                        }
                    ++col_idx_new;
                }
            return return_matrix;
        }


    /**
     * Convert this matrix to a string.
     **/
    std::string to_string() const
        {
            std::stringstream output_str_buf;
            output_str_buf << "Sequence name\t";
            std::vector<LabelType> row_labs;
            std::vector<LabelType> col_labs;

            row_labs.reserve( this->row_labels.size() );
            col_labs.reserve( this->col_labels.size() );

            this->get_row_labels( row_labs );
            this->get_col_labels( col_labs );

            output_str_buf << boost::algorithm::join( col_labs, "\t" ) << "\n";

            for( const auto& row_lab : row_labs  )
                {
                    output_str_buf << row_lab << "\t";

                    std::uint32_t col_idx = 0;
                    for( const auto& col_lab : col_labs )
                        {
                            output_str_buf << this->operator()( row_lab, col_lab );

                            if( col_idx < this->ncols() - 1 )
                                {
                                    output_str_buf << "\t";
                                }
                            ++col_idx;
                        }
                    output_str_buf << "\n";
                }

            return output_str_buf.str();
        }

    bool operator==( const labeled_matrix<ValType,LabelType>& other ) const
    {
        return matrix<ValType>::operator==( other )
               && this->row_labels == other.row_labels
                  && this->col_labels == other.col_labels;
    }

    /**
     * Transpose this matrix, returning the result.
     * @returns this^T, the transposition of this matrix.
     **/
    labeled_matrix<ValType,LabelType>
    transpose()
    {
        labeled_matrix<ValType,LabelType>
            ret_val{ this->M,
                     this->N,
                     this->col_labels,
                     this->row_labels
                   };

        for( std::uint32_t row = 0; row < this->N; ++row )
            {
                for( std::uint32_t col = 0; col < this->M; ++col )
                    {
                        ret_val( col, row ) = this->operator()( row, col );
                    }
            }

        return ret_val;
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

    /**
     * Set a label in one of this class' map members.
     * @param map The map whose label to set
     * @param original_label The original label, must be a key in map.
     * @param new_label The new label to set
    **/
    void set_label( std::unordered_map<LabelType,std::uint32_t>&
                    map,
                    const LabelType& original_label,
                    const LabelType& new_label
                  )
    {
        auto original = map.find( original_label );
        auto value = original->second;
        map.erase( original_label );
        map.emplace( new_label, value );
    }

};

template<typename Stream, typename Lab, typename Val>
Stream& operator<<( Stream& out,
                    const labeled_matrix<Lab,Val>& matr
                  )
{
    out << matr.to_string();
    return out;
}

#endif // MATRIX_HH_INCLUDED
