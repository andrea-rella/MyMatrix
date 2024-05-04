# Sparse Matrix Computation (Small) Library

## Overview

This small Matrix Library Project is a C++ template library designed to deal with 2D sparse matrices and perform various matrix operations efficiently. It can handle any type that is convertible to double and also `std::complex<T>`. Both row-major (CSR) and column-major storage orders (CSC) are supported. The library is suitable for applications in numerical analysis, linear algebra, and other fields requiring matrix manipulation.

## Brief Theoretical Background

In extensive computational tasks, leveraging sparsity is crucial when dealing with large matrices. The memory required for fully storing a matrix increases quadratically with its dimension ($N^2$), which can quickly exhaust available memory resources. However, when matrices are sparse (the number of non-zero elements grows as $N$, being $N$ the linear size of the matrix), optimizing operations such as matrix-vector multiplication becomes feasible by omitting multiplications with zero elements, thereby enhancing efficiency.

For sparse matricies we can distinguish between two families of storage techniques:

- **Dynamic (or uncompressed) storage techniques**. They allow to add (and sometimes also eliminate) new non-zero elements easily.
- **Compressed storage techniques**: They are the most efficient in terms of memory and of the computational efficiency of basic operations like matrix-vector product. But they do not allow to change the pattern of sparsity.

For more details consult [Yousef Saad. Iterative methods for sparse linear systems. SIAM, 2003](https://www-users.cse.umn.edu/~saad/IterMethBook_2ndEd.pdf) or simply have a look at [Sparse Matrices](https://en.wikipedia.org/wiki/Sparse_matrix)

## Features

- **Generic Matrix Class**: Templated to work with any data type that supports basic arithmetic operations.
- **Storage Flexibility**: Supports both dynamic (through a `std::map`) and compressed storage (CSC or CSR) to optimize different operations.
- **Compression**: Provides functionality to compress and uncompress the matrix, optimizing storage for sparse matrices.
- **Matrix Operations**: Includes operations like matrix-vector multiplication and matrix-matrix multiplication.
- **Norm Calculations**: Supports various norms (One, Infinity, Frobenius) for matrices.
- **File I/O**: Capable of reading matrices from [Matrix Market](https://math.nist.gov/MatrixMarket/) format (.mtx).

## Requirements

- C++20 compiler (GCC, Clang, MSVC)
- Standard Template Library (STL)

## Installation

No installation is needed. Include the `Matrix.hpp` and `Matrix_implementation.hpp` header in your C++ project to use the library.

## Usage

Here is a simple example of how to use the Matrix library:

```cpp
#include "Matrix.hpp"

int main() {
    algebra::Matrix<double, algebra::Ordering::RowMajor> M(10, 10); // Create a 10x10 double matrix with RowMajor storage
    M(0, 0) = 5.0; // Set element at position (0,0)
    std::vector<double> vec(10, 1.0); // Create a vector of ones
    auto result = M * vec; // Multiply matrix by vector
    mat.print(); // Print the matrix

    algebra::Matrix<double, algebra::Ordering::RowMajor> mat;
    algebra::Matrix<double, algebra::Ordering::RowMajor> result_mat;
    mat.read_from_matrix_market("pores_1.mtx"); // read from the lnsp_131.mtx https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lns/lnsp_131.html

    auto norm1 = mat.norm<algebra::Norm::One>();

    mat.compress();

    result_mat = mat * mat;
    auto result_mat
    return 0;
}
```

## Short Documentation for the Main Methods

- `compress()`: Compresses the matrix to save space, particularly useful if the matrix is sparse.

- `uncompress()`: Converts the matrix back to its uncompressed state, allowing for normal operations that may not be supported in compressed form.

- `norm()`: Calculates the specified norm (One, Infinity, Frobenius) of the matrix.

- `operator*(const Matrix<T, StorageOrder> &, const std::vector<T> &)`: Multiplies the matrix by a vector and returns the resulting vector.

- `operator*(const Matrix<T, StorageOrder> &, the Matrix<T, StorageOrder> &)`: Multiplies the matrix by another matrix and returns the resulting matrix. Works only with matrices of the same state (Compressed or Uncompressed) and with the same ordering.

## Contributing

Contributions to the library are welcome. Please feel free to fork the repository, make changes, and submit pull requests.
