#include <iostream>
#include <vector>
#include <Accelerate/Accelerate.h>

// Running a numerical factorization on a singular matrix causes a segmentation fault or address
// boundary error.
//
// The expected behaviour is that the factorization returns with an appropriate failure status.
int main() {
    std::cout << "Constructing an ordinary sparse matrix..." << std::endl;
    /// Define the sparsity pattern of matrix `A`.
    // [0 0 .]
    // [. 0 .]
    // [0 . 0]
    std::vector<int> rowIndices{1, 2, 0, 1};
    std::vector<int> columnIndices{0, 1, 2, 2};
    std::vector<long> columnStarts{0, 1, 2, 4};

    /// Create the single-precision coefficient matrix `A`.
    std::vector<float> aValues{1.0, 0.0, 1.0, 1.0};

    SparseMatrixStructure structure = {
            .rowCount = 3,
            .columnCount = 3,
            .columnStarts = columnStarts.data(),
            .rowIndices = rowIndices.data(),
            .attributes = {
                    .kind = SparseOrdinary,
            },
            .blockSize = 1,
    };

    SparseMatrix_Float A = {
            .structure = structure,
            .data = aValues.data(),
    };

    std::cout << "Symbolically factorizing using QR factorization..." << std::endl;
    /// Compute the symbolic factorization from the structure of the coefficient matrix.
    auto symbolicFactorization = SparseFactor(SparseFactorizationQR, structure);

    std::cout << "Status: " << symbolicFactorization.status << std::endl;

    std::cout << "Numerically factorizing..." << std::endl;
    /// Factorize A using the symbolic factorization.
    auto factorization = SparseFactor(symbolicFactorization, A);

    std::cout << "Status: " << factorization.status << std::endl;
}
