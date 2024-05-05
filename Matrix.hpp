#ifndef MATRIX_C76222F2_8354_4C29_AB07_1C2D3389FAC4
#define MATRIX_C76222F2_8354_4C29_AB07_1C2D3389FAC4

#include <map>
#include <array>
#include <vector>
#include <cstddef>
#include <type_traits>
#include <concepts>
#include <fstream>
#include <sstream>
#include <string>
#include <complex>

namespace algebra
{

    enum class Ordering
    {
        RowMajor,
        ColumnMajor
    };

    enum class Norm
    {
        One,
        Infinity,
        Frobenius
    };

    // Concepts needed

    /**
     * This concept uses std::is_convertible_v<T, double> to check if T can be converted to double, which includes int, long,
     * long long, and floating-point types. . The requires clause with std::enable_if and std::is_same ensures that T can only
     * be std::complex<U> for some scalar type U. This types all support the underlying operations we'll need: +, +=, * scalar,
     * std::abs() etc;
     */

    template <typename U>
    concept ScalarOrComplex = std::is_convertible_v<U, double> ||
                              std::is_same_v<std::decay_t<U>, typename std::complex<typename std::decay_t<U>::value_type>>;

    template <ScalarOrComplex T, Ordering StorageOrder> // forward declaration needed for friend method
    class Matrix;

    template <ScalarOrComplex T, Ordering StorageOrder>
    std::vector<T> operator*(const Matrix<T, StorageOrder> &, const std::vector<T> &);

    template <ScalarOrComplex T, Ordering StorageOrder>
    Matrix<T, StorageOrder> operator*(const Matrix<T, StorageOrder> &, const Matrix<T, StorageOrder> &);

    template <ScalarOrComplex T, Ordering StorageOrder>
    class Matrix
    {

    private:
        // Custom comparator for std::map that translates ColumnMajor / RowMajor by modifiyng the behaviour of '<'
        struct Comparator
        {
            bool operator()(const std::array<std::size_t, 2> &lhs, const std::array<std::size_t, 2> &rhs) const
            {
                if constexpr (StorageOrder == Ordering::ColumnMajor)
                {
                    // Comparison logic for ColumnMajor
                    if (lhs[1] < rhs[1])
                        return true;
                    if (lhs[1] > rhs[1])
                        return false;
                    return lhs[0] < rhs[0];
                }
                else
                {
                    // Default comparison logic for RowMajor
                    return lhs < rhs;
                }
            }
        };

        // Dimensions
        std::size_t Rows;
        std::size_t Cols;

        // Dynamic
        std::map<std::array<std::size_t, 2>, T, Comparator> Data;

        // Compressed
        /**
         * @note since both row and column storage order are conteplated i avoid using the term row and column,
         * intead i use inner and outer. Inner refers to the ordering thus contains the starting index for the
         * elements of each row (col) while outer contains the corresponding column (row) index. The ordering in
         * itself will be decided by the StorageOrder parameter and implemented through the Comparator
         */
        std::vector<T> Values;
        std::vector<std::size_t> Inner_ptr;
        std::vector<std::size_t> Outer_idxs;

        // State
        bool Compressed; // Check what's the default value for bool, in case modify the Default constructor

    public:
        Matrix() = default;
        Matrix(Matrix &) = default;
        Matrix(Matrix &&) = default;

        // Constructors
        Matrix(std::size_t, std::size_t);

        // Default assignements
        Matrix<T, StorageOrder> &operator=(const Matrix<T, StorageOrder> &) = default;
        Matrix<T, StorageOrder> &operator=(Matrix<T, StorageOrder> &&) = default;

        // Getters
        std::size_t N_Rows() const;
        std::size_t N_Cols() const;

        // Resize
        void resize(const std::size_t, const std::size_t);

        // Read from matrix market
        void read_from_matrix_market(const std::string);

        // Print
        void print() const;

        // Call Operators
        const T &operator()(std::size_t row, std::size_t col) const;
        T &operator()(std::size_t row, std::size_t col);

        // State Check
        bool is_compressed() const { return Compressed; };

        // Method That Compresses The Matrix
        void compress();

        // Method That Uncompresses The Matrix
        void uncompress();

        // Resize The matrix
        void resize();

        // Norm
        template <Norm normType>
        double norm() const;

        // Vec multiplication
        friend std::vector<T> operator* <T, StorageOrder>(const Matrix<T, StorageOrder> &, const std::vector<T> &);

        // Matrix by Matrix
        friend Matrix<T, StorageOrder> operator* <T, StorageOrder>(const Matrix<T, StorageOrder> &, const Matrix<T, StorageOrder> &);
    };
}

#include "Matrix_implementation.hpp"

#endif /* MATRIX_C76222F2_8354_4C29_AB07_1C2D3389FAC4 */
