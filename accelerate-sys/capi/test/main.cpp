#include <iostream>
#include <vector>
#include "../src/api.h"

int main() {
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

    std::cout << "Done" << std::endl;
}
