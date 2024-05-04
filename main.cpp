#include "Matrix.hpp"
#include <iostream>
#include <chrono>

void printVector(const std::vector<double> &);

int main()
{
    algebra::Matrix<double, algebra::Ordering::RowMajor> mat;
    algebra::Matrix<double, algebra::Ordering::RowMajor> result_mat;

    mat.read_from_matrix_market("pores_1.mtx");

    std::vector<double> v(mat.N_Cols(), 1.);
    std::vector<double> result;

    // ------------------------------------------------------------------------------
    // -----------------------         TEST FOR 'A*v'         -----------------------
    // ------------------------------------------------------------------------------

    auto start = std::chrono::high_resolution_clock::now();
    result = mat * v;
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Elapsed time for uncompressed Matrix-vector product: " << elapsed.count() << " microseconds\n \n";
    // printVector(result);

    std::cout << "\n \n \n " << std::endl;

    mat.compress();

    start = std::chrono::high_resolution_clock::now();
    result = mat * v;
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Elapsed time for compressed Matrix-vector product: " << elapsed.count() << " microseconds\n \n";
    // printVector(result);

    // ------------------------------------------------------------------------------
    // -----------------------         TEST FOR 'Norm'        -----------------------
    // ------------------------------------------------------------------------------

    auto norm1 = mat.norm<algebra::Norm::One>();
    auto norminf = mat.norm<algebra::Norm::Infinity>();
    auto normFR = mat.norm<algebra::Norm::Frobenius>();

    std::cout << "Norm One: " << norm1 << std::endl;
    std::cout << "Norm Infinty: " << norminf << std::endl;
    std::cout << "Norm Frobeinus: " << normFR << std::endl;

    std::cout << "\n \n"
              << std::endl;

    // ------------------------------------------------------------------------------
    // -----------------------         TEST FOR 'A*B'         -----------------------
    // ------------------------------------------------------------------------------

    mat.uncompress();

    start = std::chrono::high_resolution_clock::now();
    result_mat = mat * mat;
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Elapsed time for uncompressed Matrix-vector product: " << elapsed.count() << " microseconds\n \n";
    // result_mat.print();

    mat.compress();

    start = std::chrono::high_resolution_clock::now();
    result_mat = mat * mat;
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Elapsed time for compressed Matrix-vector product: " << elapsed.count() << " microseconds\n \n";
    // result_mat.uncompress();
    // result_mat.print();

    return 0;
};

void printVector(const std::vector<double> &vec)
{
    std::cout << "(";
    for (double elem : vec)
    {
        std::cout << elem << ", ";
    }
    std::cout << ")";
    std::cout << std::endl;
}