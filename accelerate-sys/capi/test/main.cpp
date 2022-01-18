#include <iostream>
#include <vector>
#include "../src/api.h"

// This test currently fails.
// This is a bug in the Accelerate framework, a report was submitted to Apple.
void singular_test() {
    /// Define the sparsity pattern of matrix `A`.
    // [0 0 .]
    // [. 0 .]
    // [0 . 0]
    std::vector<int> rowIndices{1, 2, 0, 1};
    std::vector<int> columnIndices{0, 1, 2, 2};
    std::vector<long> columnStarts{0, 1, 2, 4};

    /// Create the single-precision coefficient matrix `A`.
    std::vector<float> aValues{1.0, 0.0, 1.0, 1.0};
//    auto A = SparseConvertFromCoordinate_Float(3, 3,
//                                                6, 1,
//                                                SparseAttributes_t(),
//                                                rowIndices.data(), columnIndices.data(),
//                                                aValues.data());

    SparseMatrixStructure structure = {
            .rowCount = 3,
            .columnCount = 3,
            .columnStarts = columnStarts.data(),
            .rowIndices = rowIndices.data(),
            .attributes = {
                    .kind = SparseOrdinary,
                    //.triangle = SparseLowerTriangle,
            },
            .blockSize = 1,
    };

    SparseMatrix_Float A = {
            .structure = structure,
            .data = aValues.data(),
    };

    /// Compute the symbolic factorization from the structure of the coefficient matrix.
    auto symbolicFactorization = SparseFactorSymbolic(SparseFactorizationQR, structure);

    /// Factorize A using the symbolic factorization.
    auto factorization = SparseFactorNumeric_Float(symbolicFactorization, A);

    std::vector<float> bValues{30, 35, 100};
    {
        auto xb = DenseVector_Float { 3, bValues.data() };
        SparseSolveInPlace_Float(factorization, xb);

        std::cout << "x = ";
        for (int i = 0; i < 3; ++i) {
            std::cout << bValues[i] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "status = " << factorization.status << std::endl;
    assert(factorization.status == SparseStatus_t::SparseMatrixIsSingular);
}

void direct_test() {
    /// Define the sparsity pattern of matrices `A0` and `A1`.
    std::vector<int> rowIndices{ 0, 1, 1, 2};
    std::vector<int> columnIndices{ 2, 0, 2, 1};

    /// Create the single-precision coefficient matrix _A0_.
    std::vector<float> a0Values{10, 20, 5, 50};
    auto A0 = SparseConvertFromCoordinate_Float(3, 3,
                                                4, 1,
                                                SparseAttributes_t(),
                                                rowIndices.data(), columnIndices.data(),
                                                a0Values.data());

    /// Create the double-precision coefficient matrix _A1_.
    std::vector<double> a1Values{5, 10, 2.5, 25};
    auto A1 = SparseConvertFromCoordinate_Double(3, 3,
                                                 4, 1,
                                                 SparseAttributes_t(),
                                                 rowIndices.data(), columnIndices.data(),
                                                 a1Values.data());

    /// Compute the symbolic factorization from the structure of either coefficient matrix.
    auto structure = A0.structure;
    auto symbolicFactorization = SparseFactorSymbolic(SparseFactorizationQR, structure);

    /// Factorize _A0_ using the symbolic factorization.
    auto factorization0 = SparseFactorNumeric_Float(symbolicFactorization, A0);

    /// Solve _A0 · x = b0_ in place.
    std::vector<float> b0Values{30, 35, 100};
    {
        auto xb = DenseVector_Float { 3, b0Values.data() };
        SparseSolveInPlace_Float(factorization0, xb);
    }

    // Expected solution
    std::vector<float> expSol0{1.0, 2.0, 3.0};

    std::cout << "b0Values = [ ";
    for (int i = 0; i < 3; ++i) {
        std::cout << b0Values[i] << " ";
        assert(b0Values[i] == expSol0[i]);
    }
    std::cout << "] as expected" << std::endl;

    /// Factorize _A1_ using the symbolic factorization.
    auto factorization1 = SparseFactorNumeric_Double(symbolicFactorization, A1);

    /// Solve _A1 · x = b1_ in place.
    std::vector<double> b1Values{60, 70, 200};
    {
        auto xb = DenseVector_Double { 3, b1Values.data() };
        SparseSolveInPlace_Double(factorization1, xb);
    }

    // Expected solution
    std::vector<float> expSol1{4.0, 8.0, 12.0};

    std::cout << "b1Values = [ ";
    for (int i = 0; i < 3; ++i) {
        std::cout << b1Values[i] << " ";
        assert(b1Values[i] == expSol1[i]);
    }
    std::cout << "] as expected" << std::endl;

    SparseCleanupSparseMatrix_Float(A0);
    SparseCleanupSparseMatrix_Double(A1);
    SparseCleanupOpaqueNumeric_Float(factorization0);
    SparseCleanupOpaqueNumeric_Double(factorization1);

    std::cout << "Done With DIRECT test" << std::endl;
}

void iterative_test() {
    /// Define the sparsity pattern of matrix `A`.
    std::vector<int> rowIndices{ 0, 1, 3, 0, 1, 2, 3, 1, 2};
    std::vector<long> columnStarts{ 0, 3, 7, 9};

    /// Create the single-precision coefficient matrix A.
    std::vector<float> values{2.0, -0.2, 2.5, 1.0, 3.2, -0.1, 1.1, 1.4, 0.5};
    SparseMatrix_Float A = {
            .structure = {
                    .rowCount = 4,
                    .columnCount = 3,
                    .columnStarts = columnStarts.data(),
                    .rowIndices = rowIndices.data(),
                    .attributes = {
                            .kind = SparseOrdinary
                    },
                    .blockSize = 1
            },
            .data = values.data()
    };

    // Create RHS b vector.
    std::vector<float> bValues{1.200, 1.013, 0.205, -0.172};
    auto b = DenseVector_Float { 4, bValues.data() };

    // Create x vector initialized to zero.
    std::vector<float> xValues{0.0, 0.0, 0.0};
    auto x = DenseVector_Float { 3, xValues.data() };

    // Expected solution
    std::vector<float> expectedSolution{0.1, 0.2, 0.3};

    /// Solve Ax = b.
    auto status = SparseSolveIterativePrecond_Float( SparseLSMR(), A, b, x, SparsePreconditionerDiagScaling);
    if(status!=SparseIterativeConverged) {
        std::cerr << "Failed to converge. Returned with error " << status << std::endl;
    } else {
        float errorTol = 0.001;
        for(int i=0; i<x.count; i++) {
            if (std::abs(expectedSolution[i] - x.data[i]) >= errorTol) {
                std::cerr << "x[" << i << "]: expected: " << expectedSolution[i] << "; actual: " << x.data[i] << std::endl;
            }

            assert(std::abs(expectedSolution[i] - x.data[i]) < errorTol);
        }
    }

    std::cout << "Done with Iterative test" << std::endl;
}

int main() {
    direct_test();
    iterative_test();
    singular_test();
}
