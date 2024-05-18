#ifndef MATRIX_IMPLEMENTATION_CCC1361C_AF40_4256_829C_5ED86F87C4F0
#define MATRIX_IMPLEMENTATION_CCC1361C_AF40_4256_829C_5ED86F87C4F0
// clang-format off
#include "Matrix.hpp"
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <cassert>
#include <algorithm>
#include <cmath>

namespace algebra
{
    //---------------------------------------------------------------------------
    template <ScalarOrComplex T, Ordering StorageOrder>
    Matrix<T, StorageOrder>::Matrix(std::size_t rows, std::size_t cols) : Rows(rows), Cols(cols), Compressed(false){};
//@note take the habit of using {} insread of (): Rows{rows}, Cols{cols}, Compressed{false}. To get used to it. Here is irrelevant but in general is a good practice

    //---------------------------------------------------------------------------

    template <ScalarOrComplex T, Ordering StorageOrder>
    std::size_t Matrix<T, StorageOrder>::N_Rows() const
    {
        return Rows;
    }
    template <ScalarOrComplex T, Ordering StorageOrder>
    std::size_t Matrix<T, StorageOrder>::N_Cols() const
    {
        return Cols;
    }

    //---------------------------------------------------------------------------

    template <ScalarOrComplex T, Ordering StorageOrder>
    void Matrix<T, StorageOrder>::resize(const std::size_t new_rows, const std::size_t new_cols)
    {

        if (Compressed)
            uncompress();

        if (new_rows >= Rows && new_cols >= Cols)
        {
        
            Rows = new_rows;
            Cols = new_cols;
            return;
        }
        else
        {
            std::cerr << "Attention: resizing with smaller dimension could've have eliminated some elements" << std::endl;
            // Clear elements that are out of the new bounds
            for (auto it = Data.begin(); it != Data.end();)
            {
                const auto &[row, col] = it->first;
                if (row >= new_rows || col >= new_cols)
                    it = Data.erase(it); // Remove element and move to the next
                else
                    ++it; // Move to the next element
            }

            // Update matrix dimensions
            Rows = new_rows;
            Cols = new_cols;
        }
        //@note If the matrix was in a compressed state, now is uncompressed. You should recompress again!

        return;
    }

    //---------------------------------------------------------------------------

    template <ScalarOrComplex T, Ordering StorageOrder>
    void Matrix<T, StorageOrder>::read_from_matrix_market(const std::string filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Unable to open file");
        }

        std::string line;
        std::size_t rows, cols, non_zero_entries;
        bool header_processed = false;
        std::map<std::array<std::size_t, 2>, T, Comparator> entries;

        while (std::getline(file, line))
        {
            if (line.empty() || line[0] == '%')
                continue; // Skip comments and empty lines

            std::istringstream iss(line); // it accounts for scientific notation when dealing with floating point numbers
            if (!header_processed)
            {
                iss >> rows >> cols >> non_zero_entries;
                header_processed = true;
            }
            else
            {
                std::size_t row, col;
                T value;
                iss >> row >> col >> value;
                entries[{row - 1, col - 1}] = value; // Adjust for 0-based indexing
            }
        }

        if (entries.size() != non_zero_entries)
        {
            throw std::runtime_error("The number of non-zero entries does not match the expected count.");
        }

        file.close();

        // Clear existing data and update matrix dimensions
        Data.clear();
        Rows = rows;
        Cols = cols;
        Compressed = false;

        // Update the data member with new entries
        Data = std::move(entries); //@note nice using move here. I only do not understand why not filling Data directly

        return;
    }

    //---------------------------------------------------------------------------

    template <ScalarOrComplex T, Ordering StorageOrder>
    void Matrix<T, StorageOrder>::print() const
    {
        if (Compressed)
        {
            std::cout << "Matrix is compressed. Decompress before printing." << std::endl;
            //@note Not necessary to decompress before printing, you can print the compressed matrix as well
            return;
        }

        for (std::size_t i = 0; i < Rows; ++i)
        {
            for (std::size_t j = 0; j < Cols; ++j)
            {
                auto it = Data.find({i, j});
                if (it != Data.end())
                {
                    std::cout << it->second << "    ";
                }
                else
                {
                    std::cout << "0    "; // Print zero for non-existent elements in sparse representation
                }
            }
            std::cout << "\n"
                      << std::endl;
        }
    }

    //---------------------------------------------------------------------------

    template <ScalarOrComplex T, Ordering StorageOrder>
    const T &Matrix<T, StorageOrder>::operator()(std::size_t row, std::size_t col) const
    {
        //@note There is no need to return a const T&. You can return directly a T, whihc is the most common treturn type for const methods
        if (row >= Rows || col >= Cols)
            throw std::out_of_range("Matrix indices out of range");

        std::size_t inner, outer;

        if constexpr (StorageOrder == Ordering::RowMajor) //@note nice
        {
            inner = row;
            outer = col;
        }
        else
        {
            inner = col;
            outer = row;
        }

        if (Compressed)
        {
            // Compressed state
            //@note If the inner indexes are ordered you can use the binary search algorithm to find the element
            for (std::size_t i = Inner_ptr[inner]; i < Inner_ptr[inner + 1]; ++i)
            {
                if (Outer_idxs[i] == outer)
                {
                    return Values[i]; // Found element
                }
            }
            return T(); // Element not found (value is 0) @note better to return T{0}, to be sure! It works also for complex numbers
        }
        else
        {
            // Uncompressed state
            auto it = Data.find({row, col});

            if (it != Data.end())
                return it->second; // Element found
            else
                return T(); // Element not found (value is 0)
        }
    }

    //---------------------------------------------------------------------------

    template <ScalarOrComplex T, Ordering StorageOrder>
    T &Matrix<T, StorageOrder>::operator()(std::size_t row, std::size_t col)
    {
        if (row >= Rows || col >= Cols)
            throw std::out_of_range("Matrix indices out of range");

        if (Compressed)
            throw std::logic_error("Cannot modify elements inner compressed state");
        else
        {
            auto it = Data.find({row, col});

            if (it != Data.end())
                return it->second; // Return reference to existing element
            else
                // Element not found, add new element
                return Data[{row, col}]; // Creates a new element with default value (0 for numeric types)
                //@note if you want to be sure that the element is zero, you can use return (Data[{row, col}] = T{0}); instead
        }
    }

    //---------------------------------------------------------------------------

    template <ScalarOrComplex T, Ordering StorageOrder>
    void Matrix<T, StorageOrder>::compress()
    {
        // Already compressed or no data to compress
        if (Compressed)
            return;

        if (Data.empty())
        {
            Compressed = true;
            return;
        }

        // the field corresponding to the compressed case will for sure be empty given the implementation of all
        // the other methods so there is no need to clear them before proceeding

        std::size_t current_inner = 0;
        Inner_ptr.push_back(0); // Start of first inner

        auto it = Data.begin();

        if constexpr (StorageOrder == Ordering::RowMajor)
        {
            while (it != Data.end())
            {
                const auto &index = it->first;
                std::size_t inner = index[0];
                std::size_t outer = index[1];
                T value = it->second;

                // Update Inner_idxs if necessary
                /**
                 * @note you can also think of the inner_idx storage as it stores the cumulative count
                 * of non-zero elements up to each row(col)
                 */
                if (inner > current_inner)
                {
                    Inner_ptr.push_back(Values.size()); // Start of next inner
                    current_inner = inner;              // Update current_inner to the new inner index
                }

                // Store value and column index
                Values.push_back(value);
                Outer_idxs.push_back(outer);

                ++it;
            }
        }
        else
        {
            while (it != Data.end())
            {
                const auto &index = it->first;
                std::size_t inner = index[1];
                std::size_t outer = index[0];
                T value = it->second;

                if (inner > current_inner)
                {
                    Inner_ptr.push_back(Values.size());
                    current_inner = inner;
                }

                // Store value and column index
                Values.push_back(value);
                Outer_idxs.push_back(outer);

                ++it;
            }
        }

        // Add one more entry to Inner_idxs to mark the end of the last inner
        Inner_ptr.push_back(Values.size());

        Data.clear();

        Compressed = true;
    }

    //---------------------------------------------------------------------------

    template <ScalarOrComplex T, Ordering StorageOrder>
    void Matrix<T, StorageOrder>::uncompress()
    {
        // Already decompressed
        if (!Compressed)
            return;

        if (Values.empty() && Inner_ptr.empty() && Outer_idxs.empty())
        {
            Compressed = false;
            return;
        }

        for (std::size_t i = 0; i < Inner_ptr.size() - 1; ++i)
        {
            std::size_t start = Inner_ptr[i];
            std::size_t end = Inner_ptr[i + 1];
            for (std::size_t j = start; j < end; ++j)
            {
                std::size_t inner = i;
                std::size_t outer = Outer_idxs[j];
                if constexpr (StorageOrder == Ordering::RowMajor)
                {
                    Data.emplace(std::array<std::size_t, 2>{inner, outer}, Values[j]);//@note good useing emplace here
                }
                else if constexpr (StorageOrder == Ordering::ColumnMajor)
                {
                    Data.emplace(std::array<std::size_t, 2>{outer, inner}, Values[j]);
                }
            }
        }

        Values.clear();
        Inner_ptr.clear();
        Outer_idxs.clear();
        Compressed = false;
    }

    //---------------------------------------------------------------------------

    template <ScalarOrComplex T, Ordering StorageOrder>
    template <Norm normType>
    double Matrix<T, StorageOrder>::norm() const
    {
//@note good work, but normally higher order operations assume that the matrix is compressed, for simplicity
        if (!Compressed)
        {
            if constexpr (normType == Norm::One)
            {
                std::vector<double> v;
                v.reserve(Cols);
                std::fill_n(std::back_inserter(v), Cols, 0.);
                double max_value = 0.0;

                for (const auto &entry : Data)
                {
                    const auto [i, j] = entry.first;
                    v[j] += std::abs(entry.second);
                    max_value = std::max(max_value, v[j]);
                }

                return max_value;
            }
            else if constexpr (normType == Norm::Infinity)
            {
                std::vector<double> v;
                v.reserve(Rows);
                std::fill_n(std::back_inserter(v), Rows, 0.);

                double max_value = 0.0;

                for (const auto &entry : Data)
                {
                    const auto [i, j] = entry.first;
                    v[i] += std::abs(entry.second);
                    max_value = std::max(max_value, v[i]);
                }

                return max_value;
            }
            else if constexpr (normType == Norm::Frobenius)
            {//@note good use of standard algorithms
                double sum = std::transform_reduce(Data.begin(), Data.end(), 0.0, std::plus<>(),
                                                   [](const auto &entry)
                                                   {
                                                       return std::abs(entry.second) * std::abs(entry.second);
                                                   });

                return std::sqrt(sum);
            }
        }
        else // ---- so if is compressed -----
        {
            if constexpr (normType == Norm::One)
            {
                std::vector<double> v;
                v.reserve(Rows);
                std::fill_n(std::back_inserter(v), Rows, 0.);

                double max_value = 0.0;

                if constexpr (StorageOrder == Ordering::RowMajor)
                {
                    std::vector<double> col_sums(Cols, 0.0);
                    for (std::size_t i = 0; i < Inner_ptr.size() - 1; ++i)
                    {
                        for (std::size_t j = Inner_ptr[i]; j < Inner_ptr[i + 1]; ++j)
                        {
                            col_sums[Outer_idxs[j]] += std::abs(Values[j]);
                        }
                    }
                    max_value = *std::max_element(col_sums.begin(), col_sums.end());
                }
                else if constexpr (StorageOrder == Ordering::ColumnMajor)
                {
                    std::vector<double> row_sums(Rows, 0.0);
                    for (std::size_t i = 0; i < Inner_ptr.size() - 1; ++i)
                    {
                        for (std::size_t j = Inner_ptr[i]; j < Inner_ptr[i + 1]; ++j)
                        {
                            row_sums[i] += std::abs(Values[j]);
                        }
                    }
                    max_value = *std::max_element(row_sums.begin(), row_sums.end());
                }

                return max_value;
            }
            else if constexpr (normType == Norm::Infinity)
            {

                std::vector<double> v;
                v.reserve(Rows);
                std::fill_n(std::back_inserter(v), Rows, 0.);
                double max_value = 0.0;

                if constexpr (StorageOrder == Ordering::RowMajor)
                {
                    for (std::size_t i = 0; i < Inner_ptr.size() - 1; ++i)
                    {
                        double row_sum = 0.0;
                        for (std::size_t j = Inner_ptr[i]; j < Inner_ptr[i + 1]; ++j)
                        {
                            row_sum += std::abs(Values[j]);
                        }
                        max_value = std::max(max_value, row_sum);
                    }
                }
                else if constexpr (StorageOrder == Ordering::ColumnMajor)
                {
                    std::vector<double> col_sums(Cols, 0.0);
                    for (std::size_t i = 0; i < Inner_ptr.size() - 1; ++i)
                    {
                        for (std::size_t j = Inner_ptr[i]; j < Inner_ptr[i + 1]; ++j)
                        {
                            col_sums[Outer_idxs[j]] += std::abs(Values[j]);
                        }
                    }
                    max_value = *std::max_element(col_sums.begin(), col_sums.end());
                }

                return max_value;
            }
            else if constexpr (normType == Norm::Frobenius)
            {

                double sum = std::transform_reduce(Values.begin(), Values.end(), 0.0, std::plus<>(),
                                                   [](const auto &entry)
                                                   {
                                                       return std::abs(entry) * std::abs(entry);
                                                   });

                return std::sqrt(sum);
            }
        }
    }

    //---------------------------------------------------------------------------

    /**
     * @note In the algorithm, when iterating over each row of the matrix, the row start and end indices in the
     * values array are determined by the row pointer array. If a row contains all zeros, its start and end indices
     * in the values array will be the same. Therefore, the loop that multiplies non-zero values by corresponding
     * elements in the vector won't execute for that row because the start and end indices will be equal, effectively
     * skipping over the row.
     * 
     * @note Teacher: The previous note is correct. One of the advantages of a sparse format is that you skip the zero elements
     *
     * @note i decided not to use operator() to save on the overhead of calling a function everytime
     * 
     * @note Teacher: and this is a good choice. In general, you should avoid using operator() for performance reasons
     */

    template <ScalarOrComplex T, Ordering StorageOrder>
    std::vector<T> operator*(const Matrix<T, StorageOrder> &M, const std::vector<T> &v)
    {
        if (M.Cols != v.size())
            throw std::invalid_argument("Vector size must match the number of matrix columns.");

        std::vector<T> result(M.Rows, T()); // i use as starting value the default value (which is zero for double, longdouble, int etc)

        if (M.Compressed)
        {

            for (std::size_t i = 0; i < M.Inner_ptr.size() - 1; ++i) // For each row (col) of the matrix
            {
                std::size_t inner_start = M.Inner_ptr[i]; // Find the values in said row (col)
                std::size_t inner_end = M.Inner_ptr[i + 1];

                if constexpr (StorageOrder == Ordering::RowMajor)
                {

                    for (std::size_t j = inner_start; j < inner_end; ++j) // Iterate over such elements perform matrix vector comutation (b_i = \sum_{k=1}^Col A_ik * v_k)

                        result[i] += M.Values[j] * v[M.Outer_idxs[j]];
                }
                else if constexpr (StorageOrder == Ordering::ColumnMajor)
                {

                    for (std::size_t j = inner_start; j < inner_end; ++j) // Iterate through each non-zero element in the column

                        result[M.Outer_idxs[j]] += M.Values[j] * v[i]; // Accumulate the result by multiplying the non-zero element with the corresponding element of the vector
                }
            }
        }
        else
        {
            for (const auto &entry : M.Data)
            {
                const auto [i, j] = entry.first;
                result[i] += entry.second * v[j];
            }
        }

        return result;
    }

    //---------------------------------------------------------------------------

    // I ask the same ordering, could be extended to cover all possible combinations
    template <ScalarOrComplex T, Ordering StorageOrder>
    Matrix<T, StorageOrder> operator*(const Matrix<T, StorageOrder> &matrix_A, const Matrix<T, StorageOrder> &matrix_B)
    {
        if (matrix_A.Cols != matrix_B.Rows)
            throw std::invalid_argument("Dimension mismatch.");

        // This is a a-priori constraint, could be extended to cover all possible combinations
        if ((matrix_A.Compressed && !matrix_B.Compressed) || (matrix_B.Compressed && !matrix_A.Compressed))
            throw std::logic_error("Make sure they are both in the same state before proceeding with Matrix Multiplication.");

        Matrix<T, StorageOrder> matrix_C(matrix_A.Rows, matrix_B.Cols); // it's a empty matrix formally in a uncompressed state

        if (matrix_A.Compressed && matrix_B.Compressed)
        {
            /** @note
             * c_ij = \sum_k a_ik* b_kj
             *
             * again as in the matrix vector product we play with the bact that A*B = (B^T * A^T)^T using
             * c_ij = \sum_k a_ik* b_kj for RowMajor
             * c_ij = o_ji = \sum_k b_jk * a_ki where O = B^T*A^T for ColMajor
             * since for column major accessing the rows is not as efficient so we play with the fact that in trasposing the columns become
             * the rows and viceversa
             * 
             * @note Teacher: good! You are exploiting the properties of the matrix multiplication
             */
            matrix_C.compress();

            if constexpr (StorageOrder == Ordering::RowMajor)
            {
                matrix_C.Inner_ptr.push_back(0); // Initialize the start of the first row

                // for each row of A
                for (std::size_t i = 0; i < matrix_A.Rows; ++i)
                {
                    std::vector<T> temp_row(matrix_B.Cols, T());

                    // Loop through non-zero elements of the i-th row of A
                    for (std::size_t n = matrix_A.Inner_ptr[i]; n < matrix_A.Inner_ptr[i + 1]; ++n)
                    {
                        std::size_t k = matrix_A.Outer_idxs[n];
                        T val_A = matrix_A.Values[n];

                        // Loop through non zero elements of the k-th row of B
                        for (std::size_t m = matrix_B.Inner_ptr[k]; m < matrix_B.Inner_ptr[k + 1]; ++m)
                        {
                            std::size_t j = matrix_B.Outer_idxs[m];
                            T val_B = matrix_B.Values[m];

                            temp_row[j] += val_A * val_B;
                        }
                    }

                    for (std::size_t c = 0; c < matrix_B.Cols; ++c)
                    {
                        if (temp_row[c] != T())
                        {
                            matrix_C.Values.push_back(temp_row[c]);
                            matrix_C.Outer_idxs.push_back(c);
                        }
                    }

                    matrix_C.Inner_ptr.push_back(matrix_C.Values.size()); // since each inner_ptr can be seen as the cumulative element of each row (col in case of ColMajor)
                }
            }
            else if constexpr (StorageOrder == Ordering::ColumnMajor)
            {

                matrix_C.Inner_ptr.push_back(0); // Initialize the start of the first col

                // for each col of B (row of B^T)
                for (std::size_t i = 0; i < matrix_B.Cols; ++i)
                {
                    std::vector<T> temp_col(matrix_A.Rows, T());

                    // Loop through non-zero elements of the i-th column of B
                    for (std::size_t n = matrix_B.Inner_ptr[i]; n < matrix_B.Inner_ptr[i + 1]; ++n)
                    {
                        std::size_t k = matrix_B.Outer_idxs[n];
                        T val_B = matrix_B.Values[n];

                        // Loop through non-zero elements of the k-th column of A
                        for (std::size_t m = matrix_A.Inner_ptr[k]; m < matrix_A.Inner_ptr[k + 1]; ++m)
                        {
                            std::size_t j = matrix_A.Outer_idxs[m];
                            T val_A = matrix_A.Values[m];

                            temp_col[j] += val_A * val_B;
                        }
                    }

                    for (std::size_t c = 0; c < matrix_A.Rows; ++c)
                    {
                        if (temp_col[c] != T())
                        {
                            matrix_C.Values.push_back(temp_col[c]);
                            matrix_C.Outer_idxs.push_back(c);
                        }
                    }

                    matrix_C.Inner_ptr.push_back(matrix_C.Values.size());
                }
            }
        }
        else
        {
            for (const auto &entryA : matrix_A.Data)
            {
                const auto &indexA = entryA.first;
                const T &valueA = entryA.second;

                for (const auto &entryB : matrix_B.Data)
                //@note you can use more modern cunstructs like for (auto const & [indexB, valueB] : matrix_B.Data) 
                // and then auto const & [i,j]=indexB;
                // to make the code more readable;
                {
                    const auto &indexB = entryB.first;
                    const T &valueB = entryB.second;

                    if (indexA[1] == indexB[0])
                    {
                        matrix_C.Data[{indexA[0], indexB[1]}] += valueA * valueB;
                    }
                }
            }
        }

        return std::move(matrix_C); 
        //@note you gain nothing from the move. Returned values are rvlaues, so they will be moved when possible
        //Actually, with std::move you block copy elision, so it may be even slower. Just return matrix_C; It is the best way.
    }

    //---------------------------------------------------------------------------
}

#endif /* MATRIX_IMPLEMENTATION_CCC1361C_AF40_4256_829C_5ED86F87C4F0 */
