#include "api.h"


SparseOpaqueFactorization_Double SparseFactor_Double(SparseFactorization_t type, SparseMatrix_Double Matrix) {
    return SparseFactor(type, Matrix);
}

SparseOpaqueFactorization_Float SparseFactor_Float(SparseFactorization_t type, SparseMatrix_Float Matrix) {
    return SparseFactor(type, Matrix);
}

SparseOpaqueFactorization_Double SparseFactorOpt_Double(SparseFactorization_t type, SparseMatrix_Double Matrix, SparseSymbolicFactorOptions sfoptions, SparseNumericFactorOptions nfoptions) {
    return SparseFactor(type, Matrix, sfoptions, nfoptions);
}

SparseOpaqueFactorization_Float SparseFactorOpt_Float(SparseFactorization_t type, SparseMatrix_Float Matrix, SparseSymbolicFactorOptions sfoptions, SparseNumericFactorOptions nfoptions) {
    return SparseFactor(type, Matrix, sfoptions, nfoptions);
} 

/**** Sparse Factor Functions using pre-calculated symbolic factors ***********/

SparseOpaqueFactorization_Double SparseFactorNumeric_Double(SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Double Matrix) {
    return SparseFactor(SymbolicFactor, Matrix);
} 

SparseOpaqueFactorization_Float SparseFactorNumeric_Float(SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Float Matrix) {
    return SparseFactor(SymbolicFactor, Matrix);
}

SparseOpaqueFactorization_Double SparseFactorNumericOpt_Double(SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Double Matrix, SparseNumericFactorOptions nfoptions) {
    return SparseFactor(SymbolicFactor, Matrix, nfoptions);
}

SparseOpaqueFactorization_Float SparseFactorNumericOpt_Float(SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Float Matrix, SparseNumericFactorOptions nfoptions) {
    return SparseFactor(SymbolicFactor, Matrix, nfoptions);
}

/**** Sparse Factor Functions with user-defined workspace *********************/

SparseOpaqueFactorization_Double SparseFactorNumericOptWS_Double(
  SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Double Matrix,
  SparseNumericFactorOptions nfoptions, void *_Nullable factorStorage,
  void *_Nullable workspace) {
    return SparseFactor(SymbolicFactor, Matrix, nfoptions, factorStorage, workspace);
}


SparseOpaqueFactorization_Float SparseFactorNumericOptWS_Float(
  SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Float Matrix,
  SparseNumericFactorOptions nfoptions, void *_Nullable factorStorage,
  void *_Nullable workspace) {
    return SparseFactor(SymbolicFactor, Matrix, nfoptions, factorStorage, workspace);
}

SparseMatrix_Double SparseConvertFromCoordinate_Double(int rowCount, int columnCount,
  long blockCount, uint8_t blockSize, SparseAttributes_t attributes,
  const int *row, const int *column, const double *data) {

  return SparseConvertFromCoordinate(rowCount, columnCount, blockCount, blockSize, attributes,row, column,data);
}

SparseMatrix_Float SparseConvertFromCoordinate_Float(int rowCount, int columnCount, long blockCount,
  uint8_t blockSize, SparseAttributes_t attributes, const int *row,
  const int *column, const float *data) {
  return SparseConvertFromCoordinate(rowCount, columnCount, blockCount, blockSize, attributes,row, column,data);
}

SparseMatrix_Double SparseConvertFromCoordinateWS_Double(int rowCount, int columnCount, long blockCount,
  uint8_t blockSize, SparseAttributes_t attributes, const int *row,
  const int *column, const double *data, void *storage, void *workspace) {

  return SparseConvertFromCoordinate(rowCount, columnCount, blockCount, blockSize, attributes,row, column,data, storage, workspace);
}

SparseMatrix_Float SparseConvertFromCoordinateWS_Float(int rowCount, int columnCount, long blockCount,
  uint8_t blockSize, SparseAttributes_t attributes, const int *row,
  const int *column, const float *data, void *storage, void *workspace) {

  return SparseConvertFromCoordinate(rowCount, columnCount, blockCount, blockSize, attributes,row, column,data, storage, workspace);
}

/******************************************************************************
 *  @group Sparse Direct Solve Functions (DenseMatrix)
 ******************************************************************************/
void SparseSolveMatrixInPlace_Double(SparseOpaqueFactorization_Double Factored, DenseMatrix_Double XB) {
    SparseSolve(Factored, XB);
}

void SparseSolveMatrixInPlace_Float(SparseOpaqueFactorization_Float Factored, DenseMatrix_Float XB) {
    SparseSolve(Factored, XB);
}

void SparseSolveMatrix_Double(SparseOpaqueFactorization_Double Factored, DenseMatrix_Double B, DenseMatrix_Double X) {
    SparseSolve(Factored, B, X);
}

void SparseSolveMatrix_Float(SparseOpaqueFactorization_Float Factored, DenseMatrix_Float B, DenseMatrix_Float X) {
    SparseSolve(Factored, B, X);
}

/**** Solving Systems with User Defined Workspace *****************************/

void SparseSolveMatrixInPlaceWS_Double(SparseOpaqueFactorization_Double Factored,
  DenseMatrix_Double XB, void *workspace) {
    SparseSolve(Factored, XB, workspace);
}

void SparseSolveMatrixInPlaceWS_Float(SparseOpaqueFactorization_Float Factored, DenseMatrix_Float XB,
  void *workspace) {
    SparseSolve(Factored, XB, workspace);
}

void SparseSolveMatrixWS_Double(SparseOpaqueFactorization_Double Factored,
  DenseMatrix_Double X, DenseMatrix_Double B, void *workspace) {
    SparseSolve(Factored, B, X, workspace);
}

void SparseSolveMatrixWS_Float(SparseOpaqueFactorization_Float Factored, DenseMatrix_Float X,
  DenseMatrix_Float B, void *workspace) {
    SparseSolve(Factored, B, X, workspace);
}

/******************************************************************************
 *  @group Sparse Direct Solve Functions (DenseVector)
 ******************************************************************************/

void SparseSolveInPlace_Double(SparseOpaqueFactorization_Double Factored,
  DenseVector_Double xb) {
    SparseSolve(Factored, xb);
}

void SparseSolveInPlace_Float(SparseOpaqueFactorization_Float Factored, DenseVector_Float xb) {
    SparseSolve(Factored, xb);
}

void SparseSolve_Double(SparseOpaqueFactorization_Double Factored,
  DenseVector_Double b, DenseVector_Double x) {
    SparseSolve(Factored, b, x);
}

void SparseSolve_Float(SparseOpaqueFactorization_Float Factored, DenseVector_Float b,
  DenseVector_Float x) {
    SparseSolve(Factored, b, x);
}

/**** Solving Systems with User Defined Workspace *****************************/

void SparseSolveInPlaceWS_Double(SparseOpaqueFactorization_Double Factored,
  DenseVector_Double xb, void *workspace) {
    SparseSolve(Factored, xb, workspace);
}

void SparseSolveInPlaceWS_Float(SparseOpaqueFactorization_Float Factored, DenseVector_Float xb,
                 void *workspace) {
    SparseSolve(Factored, xb, workspace);
}
void SparseSolveWS_Double(SparseOpaqueFactorization_Double Factored,
  DenseVector_Double x, DenseVector_Double b, void *workspace) {
    SparseSolve(Factored, x, b, workspace);
}

void SparseSolveWS_Float(SparseOpaqueFactorization_Float Factored, DenseVector_Float x,
  DenseVector_Float b, void *workspace) {
    SparseSolve(Factored, x, b, workspace);
}


SparseOpaqueSymbolicFactorization SparseFactorSymbolic(SparseFactorization_t type,
  SparseMatrixStructure Matrix) {

    return SparseFactor(type, Matrix);
}

SparseOpaqueSymbolicFactorization SparseFactorSymbolicOpt(SparseFactorization_t type,
  SparseMatrixStructure Matrix, SparseSymbolicFactorOptions sfoptions) {
    return SparseFactor(type, Matrix, sfoptions);
}

/**** Cleaning up resources ***************************************************/

void SparseCleanupOpaqueSymbolic(SparseOpaqueSymbolicFactorization toFree) {
    SparseCleanup(toFree);
}
void SparseCleanupOpaqueNumeric_Double(SparseOpaqueFactorization_Double toFree) {
    SparseCleanup(toFree);
}
void SparseCleanupOpaqueNumeric_Float(SparseOpaqueFactorization_Float toFree) {
    SparseCleanup(toFree);
}
void SparseCleanupOpaqueSubfactor_Double(SparseOpaqueSubfactor_Double toFree) {
    SparseCleanup(toFree);
}
void SparseCleanupOpaqueSubfactor_Float(SparseOpaqueSubfactor_Float toFree) {
    SparseCleanup(toFree);
}
void SparseCleanupSparseMatrix_Double(SparseMatrix_Double toFree) {
    SparseCleanup(toFree);
}
void SparseCleanupSparseMatrix_Float(SparseMatrix_Float toFree) {
    SparseCleanup(toFree);
}
void SparseCleanupOpaquePreconditioner_Double(SparseOpaquePreconditioner_Double Preconditioner) {
    SparseCleanup(Preconditioner);
}
void SparseCleanupOpaquePreconditioner_Float(SparseOpaquePreconditioner_Float Preconditioner) {
    SparseCleanup(Preconditioner);
}
