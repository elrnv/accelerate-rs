/*
 * The documentation in this file is taken from Apple's Sparse/Solve.h. It is reproduced here for
 * convenience.
 */

#ifndef __ACCELERATE_SYS_C_API_H__
#define __ACCELERATE_SYS_C_API_H__

#include <Accelerate/Accelerate.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! @abstract Returns the specified factorization of a sparse matrix of double
 *            values.
 *
 *  @param type The type of factorization to perform.
 *
 *  @param Matrix The matrix to factorize.
 *
 *  @returns Factorization of Matrix.                                         */
SparseOpaqueFactorization_Double SparseFactor_Double(SparseFactorization_t type, SparseMatrix_Double Matrix);

/*! @abstract Returns the specified factorization of a sparse matrix of float
 *            values.
 *
 *  @param type The type of factorization to perform.
 *
 *  @param Matrix The matrix to factorize.                                    */
SparseOpaqueFactorization_Float SparseFactor_Float(SparseFactorization_t type, SparseMatrix_Float Matrix);

/*! @abstract Returns the specified factorization of a sparse matrix of double
 *            values, using the specified options.
 *
 *  @param type The type of factorization to perform.
 *
 *  @param Matrix The matrix to factorize.
 *
 *  @param sfoptions Symbolic factor options, for example the ordering algorithm
 *         to use.
 *
 *  @param nfoptions Numeric factor options, for example pivoting parameters.
 *
 *  @returns Factorization of Matrix.                                         */
SparseOpaqueFactorization_Double SparseFactorOpt_Double(SparseFactorization_t type, SparseMatrix_Double Matrix, SparseSymbolicFactorOptions sfoptions, SparseNumericFactorOptions nfoptions);

/*! @abstract Returns the specified factorization of a sparse matrix of float
 *            values, using the specified options.
 *
 *  @param type The type of factorization to perform.
 *
 *  @param Matrix The matrix to factorize.
 *
 *  @param sfoptions Symbolic factor options, for example the ordering algorithm
 *         to use.
 *
 *  @param nfoptions Numeric factor options, for example pivoting parameters.
 *
 *  @returns Factorization of Matrix.                                         */
SparseOpaqueFactorization_Float SparseFactorOpt_Float(SparseFactorization_t type, SparseMatrix_Float Matrix, SparseSymbolicFactorOptions sfoptions, SparseNumericFactorOptions nfoptions);

/**** Sparse Factor Functions using pre-calculated symbolic factors ***********/

/*! @abstract Returns the factorization of a sparse matrix of double values
 *            corresponding to the supplied symbolic factorization.
 *
 *  @param SymbolicFactor A symbolic factorization, as returned by a call of the
 *         form SymbolicFactor = SparseFactor(Matrix.structure).
 *
 *  @param Matrix The matrix to factorize.
 *
 *  @returns Factorization of Matrix.                                         */
SparseOpaqueFactorization_Double SparseFactorNumeric_Double(SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Double Matrix);

/*! @abstract Returns the factorization of a sparse matrix of float values
 *            corresponding to the supplied symbolic factorization.
 *
 *  @param SymbolicFactor A symbolic factorization, as returned by a call of the
 *         form SymbolicFactor = SparseFactor(Matrix.structure).
 *
 *  @param Matrix The matrix to factorize.
 *
 *  @returns Factorization of Matrix.                                         */
SparseOpaqueFactorization_Float SparseFactorNumeric_Float(SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Float Matrix);

/*! @abstract Returns the factorization of a sparse matrix of double values
 *            corresponding to the supplied symbolic factorization, using the
 *            specified options.
 *
 *  @param SymbolicFactor A symbolic factorization, as returned by a call of the
 *         form SymbolicFactor = SparseFactor(Matrix.structure).
 *
 *  @param Matrix The matrix to factorize.
 *
 *  @param nfoptions Numeric factor options, for example pivoting parameters.
 *
 *  @returns Factorization of Matrix.                                         */
SparseOpaqueFactorization_Double SparseFactorNumericOpt_Double(SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Double Matrix, SparseNumericFactorOptions nfoptions);

/*! @abstract Returns the factorization of a sparse matrix of float values
 *            corresponding to the supplied symbolic factorization, using the
 *            specified options.
 *
 *  @param SymbolicFactor A symbolic factorization, as returned by a call of the
 *         form SymbolicFactor = SparseFactor(Matrix.structure).
 *
 *  @param Matrix The matrix to factorize.
 *
 *  @param nfoptions Numeric factor options, for example pivoting parameters.
 *
 *  @returns Factorization of Matrix.                                         */
SparseOpaqueFactorization_Float SparseFactorNumericOpt_Float(SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Float Matrix, SparseNumericFactorOptions nfoptions);

/**** Sparse Factor Functions with user-defined workspace *********************/

/*! @abstract Returns the factorization of a sparse matrix of double values
 *            corresponding to the supplied symbolic factorization, using the
 *            specified options and without any internal memory allocations.
 *
 *  @discussion Note that internal memory allocations may occur in the case of
 *  pivoted factorizations that result in delayed pivots. If you require closer
 *  control over memory allocations, supply a sfoptions.malloc() function that
 *  implements the required behaviour, or use an alternative non-pivoted
 *  factorization returns. Note that if sfoptions.malloc() returns NULL the
 *  factorization will abort immediately.
 *
 *  @param SymbolicFactor A symbolic factorization, as returned by a call of the
 *         form SymbolicFactor = SparseFactor(Matrix.structure).
 *
 *  @param Matrix The matrix to factorize.
 *
 *  @param nfoptions Numeric factor options, for example pivoting parameters.
 *
 *  @param factorStorage A pointer to space used to store the factorization
 *         of size at least SymbolicFactor.factorSize_Double bytes. This storage
 *         should not be altered by the user during the lifetime of the return
 *         value. This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 *
 *  @param workspace A pointer to a workspace of size at least
 *         SymbolicFactor.workspaceSize_Double bytes. This workspace may be
 *         reused or destroyed by the user as soon as the function returns.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 *
 *  @returns Factorization of Matrix.                                         */
SparseOpaqueFactorization_Double SparseFactorNumericOptWS_Double(
  SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Double Matrix,
  SparseNumericFactorOptions nfoptions, void *_Nullable factorStorage,
  void *_Nullable workspace);

/*! @abstract Returns the factorization of a sparse matrix of float values
 *            corresponding to the supplied symbolic factorization, using the
 *            specified options and without any internal memory allocations.
 *
 *  @discussion Note that internal memory allocations may occur in the case of
 *  pivoted factorizations that result in delayed pivots. If you require closer
 *  control over memory allocations, supply a sfoptions.malloc() function that
 *  implements the required behaviour, or use an alternative non-pivoted
 *  factorization returns. Note that if sfoptions.malloc() returns NULL the
 *  factorization will abort immediately.
 *
 *  @param SymbolicFactor A symbolic factorization, as returned by a call of the
 *         form SymbolicFactor = SparseFactor(Matrix.structure).
 *
 *  @param Matrix The matrix to factorize.
 *
 *  @param nfoptions Numeric factor options, for example pivoting parameters.
 *
 *  @param factorStorage A pointer to space used to store the factorization
 *         of size at least SymbolicFactor.factorSize_Float bytes. This storage
 *         should not be altered by the user during the lifetime of the return
 *         value. This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 *
 *  @param workspace A pointer to a workspace of size at least
 *         SymbolicFactor.workspaceSize_Float bytes. This workspace may be
 *         reused or destroyed by the user as soon as the function returns.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 *
 *  @returns Factorization of Matrix.                                         */
SparseOpaqueFactorization_Float SparseFactorNumericOptWS_Float(
  SparseOpaqueSymbolicFactorization SymbolicFactor, SparseMatrix_Float Matrix,
  SparseNumericFactorOptions nfoptions, void *_Nullable factorStorage,
  void *_Nullable workspace);

/******************************************************************************
 *  @group Conversion From Other Formats
 ******************************************************************************
 *  @discussion In the conversion functions below, the variables `rowCount`,
 *  `columnCount`, `blockCount`, `row[]` and `col[]` describe a sparse matrix
 *  structure with `blockCount` structurally non-zero entries, each of which
 *  is a "block".
 *
 *  The matrix described has `rowCount*blockSize` rows and
 *  `columnCount*blockSize` columns. For each `i` in `0..<blockCount`, there
 *  is a a structurally non-zero block at block position `(row[i], column[i])`
 *  with numerical values `data[i*blockSize*blockSize:(i+1)*blockSize*blockSize-1]`
 *  interpreted as the elements of a dense column-major matrix with `blockSize`
 *  rows and columns.
 *
 *  If the coordinates `(row[i], column[i])` are invalid (meaning that they
 *  lie outside the ranges `0..<rowCount` or `0..<columnCount`, respectively),
 *  or `attributes.kind` is `SparseTriangular` or `SparseUnitTriangular` and
 *  they lie in the wrong triangle, then that block element is ignored and
 *  not included in the returned matrix.
 *
 *  If `attributes.kind` is `SparseSymmetric`, any entries in the wrong
 *  triangle are transposed and summed into the block at `(column[i], row[i])`
 *  if one is present.
 *
 *  In all cases, if any duplicate coordinates are present, the elements are
 *  summed and replaced with a single entry.
 *
 *  There are two variants of each converter; one which allocates its own
 *  workspace internally, and also allocates space for the resulting sparse
 *  matrix, and one which requires you to pass storage for the new matrix and
 *  a separate workspace for precise control over allocations.                */

/*! @abstract Convert from coordinate format arrays to a SparseMatrix_Double
 *  object, dropping out-of-range entries and summing duplicates.
 *
 *  @discussion For symmetric matrices, entries are accepted in either triangle
 *  (if they are in the "wrong" triangle as specified by attributes.triangle,
 *  they are transposed, and if an entry is already present, are treated as
 *  duplicates and summed).
 *  For triangular matrices, entries in the "wrong" triangle as specified by
 *  attributes.triangle are treated as out-of-range and dropped.
 *
 *  @param rowCount (input) Number of rows in structure.
 *
 *  @param columnCount (input) Number of columns in structure.
 *
 *  @param blockCount (input) Number of blocks in matrix.
 *
 *  @param blockSize (input) Block size for data storage on both input and
 *  ouput.
 *
 *  @param attributes (input) Attributes of matrix to create. The matrix will
 *  be forced to conform to the specified attributes by copying or dropping
 *  elements as needed.
 *
 *  @param row (input) Row indices of matrix structure.
 *
 *  @param column (input) Column indices of matrix structure.
 *
 *  @param data (input) The contents of the structurally non-zero (block)
 *  matrix elements.
 *
 *  @return A new SparseMatrix_Double object. When you are done with this
 *  matrix, release the memory that has been allocated by calling
 *  SparseCleanup( ) on it.                                                   */
SparseMatrix_Double SparseConvertFromCoordinate_Double(int rowCount, int columnCount,
  long blockCount, uint8_t blockSize, SparseAttributes_t attributes,
  const int *row, const int *column, const double *data);

/*! @abstract Convert from coordinate format arrays to a SparseMatrix_Float
 *  object, dropping out-of-range entries and summing duplicates.
 *
 *  @discussion For symmetric matrices, entries are accepted in either triangle
 *  (if they are in the "wrong" triangle as specified by attributes.triangle,
 *  they are transposed, and if an entry is already present, are treated as
 *  duplicates and summed).
 *  For triangular matrices, entries in the "wrong" triangle as specified by
 *  attributes.triangle are treated as out-of-range and dropped.
 *
 *  @param rowCount (input) Number of rows in structure.
 *
 *  @param columnCount (input) Number of columns in structure.
 *
 *  @param blockCount (input) Number of blocks in matrix.
 *
 *  @param blockSize (input) Block size for data storage on both input and
 *  ouput.
 *
 *  @param attributes (input) Attributes of matrix to create. The matrix will
 *  be forced to conform to the specified attributes by copying or dropping
 *  elements as needed.
 *
 *  @param row (input) Row indices of matrix structure.
 *
 *  @param column (input) Column indices of matrix structure.
 *
 *  @param data (input) The contents of the structurally non-zero (block)
 *  matrix elements.
 *
 *  @return A new SparseMatrix_Float object. When you are done with this
 *  matrix, release the memory that has been allocated by calling
 *  SparseCleanup( ) on it.                                                   */
SparseMatrix_Float SparseConvertFromCoordinate_Float(int rowCount, int columnCount, long blockCount,
  uint8_t blockSize, SparseAttributes_t attributes, const int *row,
  const int *column, const float *data);

/*! @abstract Convert from coordinate format arrays to a SparseMatrix_Double
 *  object, dropping out-of-range entries and summing duplicates.
 *
 *  @discussion For symmetric matrices, entries are accepted in either triangle 
 *  (if they are in the "wrong" triangle as specified by attributes.triangle, 
 *  they are transposed, and if an entry is already present, are treated as 
 *  duplicates and summed).
 *  For triangular matrices, entries in the "wrong" triangle as specified by
 *  attributes.triangle are treated as out-of-range and dropped.
 *
 *  @param rowCount (input) Number of rows in structure.
 *
 *  @param columnCount (input) Number of columns in structure.
 *
 *  @param blockCount (input) Number of blocks in matrix.
 *
 *  @param blockSize (input) Block size for data storage on both input and
 *  ouput.
 *
 *  @param attributes (input) Attributes of matrix to create. The matrix will
 *  be forced to conform to the specified attributes by copying or dropping
 *  elements as needed.
 *
 *  @param row (input) Row indices of matrix structure.
 *
 *  @param column (input) Column indices of matrix structure.
 *
 *  @param data (input) The contents of the structurally non-zero (block)
 *  matrix elements.
 *
 *  @param storage (output) A block of memory of size at least:
 *
 *    48 + (columnCount+1)*sizeof(long) + blockCount*sizeof(int)
 *      + blockCount*blockSize*blockSize*sizeof(double)
 *
 *  The returned structures .structure.columnStarts, .structure.rowIndices,
 *  and .data will point into this storage. You are responsible for managing
 *  the allocation and cleanup of this memory.
 *
 *  @param workspace (scratch) Workspace of size rowCount*sizeof(int).
 *
 *  @return A new SparseMatrix_Double object, using the memory you provided in
 *  the `storage` parameter.                                                  */
SparseMatrix_Double SparseConvertFromCoordinateWS_Double(int rowCount, int columnCount, long blockCount,
  uint8_t blockSize, SparseAttributes_t attributes, const int *row,
  const int *column, const double *data, void *storage, void *workspace);

/*! @abstract Convert from coordinate format arrays to a SparseMatrix_Float
 *  object, dropping out-of-range entries and summing duplicates.
 *
 *  @discussion For symmetric matrices, entries are accepted in either triangle
 *  (if they are in the "wrong" triangle as specified by attributes.triangle,
 *  they are transposed, and if an entry is already present, are treated as
 *  duplicates and summed).
 *  For triangular matrices, entries in the "wrong" triangle as specified by
 *  attributes.triangle are treated as out-of-range and dropped.
 *
 *  @param rowCount (input) Number of rows in structure.
 *
 *  @param columnCount (input) Number of columns in structure.
 *
 *  @param blockCount (input) Number of blocks in matrix.
 *
 *  @param blockSize (input) Block size for data storage on both input and
 *  ouput.
 *
 *  @param attributes (input) Attributes of matrix to create. The matrix will
 *  be forced to conform to the specified attributes by copying or dropping
 *  elements as needed.
 *
 *  @param row (input) Row indices of matrix structure.
 *
 *  @param column (input) Column indices of matrix structure.
 *
 *  @param data (input) The contents of the structurally non-zero (block)
 *  matrix elements.
 *
 *  @param storage (output) A block of memory of size at least:
 *
 *    48 + (columnCount+1)*sizeof(long) + blockCount*sizeof(int)
 *      + blockCount*blockSize*blockSize*sizeof(float)
 *
 *  The returned structures .structure.columnStarts, .structure.rowIndices,
 *  and .data will point into this storage. You are responsible for managing
 *  the allocation and cleanup of this memory.
 *
 *  @param workspace (scratch) Workspace of size rowCount*sizeof(int).
 *
 *  @return A new SparseMatrix_Float object, using the memory you provided in
 *  the `storage` parameter.                                                  */
SparseMatrix_Float SparseConvertFromCoordinateWS_Float(int rowCount, int columnCount, long blockCount,
  uint8_t blockSize, SparseAttributes_t attributes, const int *row,
  const int *column, const float *data, void *storage, void *workspace);

/******************************************************************************
 *  @group Sparse Direct Solve Functions (DenseMatrix)
 ******************************************************************************/

/*! @abstract Solves the system AX=B for X, using the supplied factorization
 *            of A, in place.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || X ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_X || AX - B ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAX=B.
 *
 *  @param Factored A factorization of A.
 *
 *  @param XB On entry, the right-hand sides B. On return, the solution vectors
 *         X. If A has dimension m x n, then XB must have dimension k x nrhs,
 *         where k=max(m,n) and nrhs is the number of right-hand sides to find
 *         solutions for.                                                     */
void SparseSolveMatrixInPlace_Double(SparseOpaqueFactorization_Double Factored, DenseMatrix_Double XB);

/*! @abstract Solves the system AX=B for X, using the supplied factorization
 *            of A, in place.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || X ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_X || AX - B ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAX=B.
 *
 *  @param Factored A factorization of A.
 *
 *  @param XB On entry, the right-hand sides B. On return, the solution vectors
 *         X. If A has dimension m x n, then XB must have dimension k x nrhs,
 *         where k=max(m,n) and nrhs is the number of right-hand sides to find
 *         solutions for.                                                     */
void SparseSolveMatrixInPlace_Float(SparseOpaqueFactorization_Float Factored, DenseMatrix_Float XB);

/*! @abstract Solves the system AX=B for X, using the supplied factorization
 *            of A, in place.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || X ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_X || AX - B ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAX=B.
 *
 *  @param Factored A factorization of A.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X Matrix in which to return solutions. If A has dimension m x n, and
 *         B has dimension m x nrhs, then X must have dimension n x nrhs.     */
void SparseSolveMatrix_Double(SparseOpaqueFactorization_Double Factored, DenseMatrix_Double B, DenseMatrix_Double X);

/*! @abstract Solves the system AX=B for X, using the supplied factorization
 *            of A, in place.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || X ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_X || AX - B ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAX=B.
 *
 *  @param Factored A factorization of A.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X Matrix in which to return solutions. If A has dimension m x n, and
 *         B has dimension m x nrhs, then X must have dimension n x nrhs.     */
void SparseSolveMatrix_Float(SparseOpaqueFactorization_Float Factored, DenseMatrix_Float B, DenseMatrix_Float X);

/**** Solving Systems with User Defined Workspace *****************************/

/*! @abstract Solves the system AX=B for X, using the supplied factorization
 *            of A, in place, and without any internal memory allocations.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || X ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_X || AX - B ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAX=B.
 *
 *  @param Factored A factorization of A.
 *
 *  @param XB On entry, the right-hand sides B. On return, the solution vectors
 *         X. If A has dimension m x n, then XB must have dimension k x nrhs,
 *         where k=max(m,n) and nrhs is the number of right-hand sides to find
 *         solutions for.
 *
 *  @param workspace Scratch space of size
 *         Factored.solveWorkspaceRequiredStatic + nrhs * Factored.solveWorkspaceRequiredPerRHS.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 */
void SparseSolveMatrixInPlaceWS_Double(SparseOpaqueFactorization_Double Factored,
  DenseMatrix_Double XB, void *workspace);

/*! @abstract Solves the system AX=B for X, using the supplied factorization
 *            of A, in place, and without any internal memory allocations.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || X ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_X || AX - B ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAX=B.
 *
 *  @param Factored A factorization of A.
 *
 *  @param XB On entry, the right-hand sides B. On return, the solution vectors
 *         X. If A has dimension m x n, then XB must have dimension k x nrhs,
 *         where k=max(m,n) and nrhs is the number of right-hand sides to find
 *         solutions for.
 *
 *  @param workspace Scratch space of size
 *         Factored.solveWorkspaceRequiredStatic + nrhs * Factored.solveWorkspaceRequiredPerRHS.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 */
void SparseSolveMatrixInPlaceWS_Float(SparseOpaqueFactorization_Float Factored, DenseMatrix_Float XB,
  void *workspace);

/*! @abstract Solves the system AX=B for X, using the supplied factorization
 *            of A, and without any internal memory allocations.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || X ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_X || AX - B ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAX=B.
 *
 *  @param Factored A factorization of A.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X Matrix in which to return solutions. If A has dimension m x n, and
 *         B has dimension m x nrhs, then X must have dimension n x nrhs.
 *
 *  @param workspace Scratch space of size
 *         Factored.solveWorkspaceRequiredStatic + nrhs * Factored.solveWorkspaceRequiredPerRHS.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 */
void SparseSolveMatrixWS_Double(SparseOpaqueFactorization_Double Factored,
  DenseMatrix_Double X, DenseMatrix_Double B, void *workspace);

/*! @abstract Solves the system AX=B for X, using the supplied factorization
 *            of A, and without any internal memory allocations.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || X ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_X || AX - B ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAX=B.
 *
 *  @param Factored A factorization of A.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X Matrix in which to return solutions. If A has dimension m x n, and
 *         B has dimension m x nrhs, then X must have dimension n x nrhs.
 *
 *  @param workspace Scratch space of size
 *         Factored.solveWorkspaceRequiredStatic + nrhs * Factored.solveWorkspaceRequiredPerRHS.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 */
void SparseSolveMatrixWS_Float(SparseOpaqueFactorization_Float Factored, DenseMatrix_Float X,
  DenseMatrix_Float B, void *workspace);

/******************************************************************************
 *  @group Sparse Direct Solve Functions (DenseVector)
 ******************************************************************************/

/*! @abstract Solves the system Ax=b for x, using the supplied factorization
 *            of A, in place.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || x ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_x || Ax - b ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAx=b.
 *
 *  @param Factored A factorization of A.
 *
 *  @param xb On entry, the right-hand side b. On return, the solution vector
 *         x. If A has dimension m x n, then xb must have length k, where
 *         k=max(m,n).                                                        */
void SparseSolveInPlace_Double(SparseOpaqueFactorization_Double Factored,
  DenseVector_Double xb);

/*! @abstract Solves the system Ax=b for x, using the supplied factorization
 *            of A, in place.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || x ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_x || Ax - b ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAx=b.
 *
 *  @param Factored A factorization of A.
 *
 *  @param xb On entry, the right-hand side b. On return, the solution vector
 *         x. If A has dimension m x n, then xb must have length k, where
 *         k=max(m,n).                                                        */
void SparseSolveInPlace_Float(SparseOpaqueFactorization_Float Factored, DenseVector_Float xb);

/*! @abstract Solves the system Ax=b for x, using the supplied factorization
 *            of A.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || x ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_x || Ax - b ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAx=b.
 *
 *  @param Factored A factorization of A.
 *
 *  @param b The right-hand side b to solve for. If A has dimension m x n, then
 *         b must have length m.
 *
 *  @param x Vector in which to return solution. If A has dimension m x n, then
 *         x must have length n.                                              */
void SparseSolve_Double(SparseOpaqueFactorization_Double Factored,
  DenseVector_Double b, DenseVector_Double x);

/*! @abstract Solves the system Ax=b for x, using the supplied factorization
 *            of A.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || x ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_x || Ax - b ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAx=b.
 *
 *  @param Factored A factorization of A.
 *
 *  @param b The right-hand side b to solve for. If A has dimension m x n, then
 *         b must have length m.
 *
 *  @param x Vector in which to return solution. If A has dimension m x n, then
 *         x must have length n.                                              */
void SparseSolve_Float(SparseOpaqueFactorization_Float Factored, DenseVector_Float b,
  DenseVector_Float x);

/**** Solving Systems with User Defined Workspace *****************************/

/*! @abstract Solves the system Ax=b for x, using the supplied factorization
 *            of A, in place.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || x ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_x || Ax - b ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAx=b.
 *
 *  @param Factored A factorization of A.
 *
 *  @param xb On entry, the right-hand side b. On return, the solution vector
 *         x. If A has dimension m x n, then xb must have length k, where
 *         k=max(m,n).
 *
 *  @param workspace Scratch space of size
 *         Factored.solveWorkspaceRequiredStatic + 1*Factored.solveWorkspaceRequiredPerRHS.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 */
void SparseSolveInPlaceWS_Double(SparseOpaqueFactorization_Double Factored,
  DenseVector_Double xb, void *workspace);

/*! @abstract Solves the system Ax=b for x, using the supplied factorization
 *            of A, in place.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || x ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_x || Ax - b ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAx=b.
 *
 *  @param Factored A factorization of A.
 *
 *  @param xb On entry, the right-hand side b. On return, the solution vector
 *         x. If A has dimension m x n, then xb must have length k, where
 *         k=max(m,n).
 *
 *  @param workspace Scratch space of size
 *         Factored.solveWorkspaceRequiredStatic + 1*Factored.solveWorkspaceRequiredPerRHS.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 */
void SparseSolveInPlaceWS_Float(SparseOpaqueFactorization_Float Factored, DenseVector_Float xb,
                 void *workspace);

/*! @abstract Solves the system Ax=b for x, using the supplied factorization
 *            of A, in place.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || x ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_x || Ax - b ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAx=b.
 *
 *  @param Factored A factorization of A.
 *
 *  @param b The right-hand side b to solve for. If A has dimension m x n, then
 *         b must have length m.
 *
 *  @param x Vector in which to return solution. If A has dimension m x n, then
 *         x must have length n.
 *
 *  @param workspace Scratch space of size
 *         Factored.solveWorkspaceRequiredStatic + 1*Factored.solveWorkspaceRequiredPerRHS.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 */
void SparseSolveWS_Double(SparseOpaqueFactorization_Double Factored,
  DenseVector_Double x, DenseVector_Double b, void *workspace);

/*! @abstract Solves the system Ax=b for x, using the supplied factorization
 *            of A, in place.
 *
 *  @discussion If the factorization is A=QR and the system is underdetermined,
 *  the solution of minimum norm || x ||_2 is returned.
 *  If the factorization is A=QR and the system is overdetermined, the least
 *  squares solution arg min_x || Ax - b ||_2 is returned.
 *  In the case of a factorization of type=SparseCholeskyAtA, the factorization
 *  is in fact of A^T A, so the solution returned is for the system A^TAx=b.
 *
 *  @param Factored A factorization of A.
 *
 *  @param b The right-hand side b to solve for. If A has dimension m x n, then
 *         b must have length m.
 *
 *  @param x Vector in which to return solution. If A has dimension m x n, then
 *         x must have length n.
 *
 *  @param workspace Scratch space of size
 *         Factored.solveWorkspaceRequiredStatic + 1*Factored.solveWorkspaceRequiredPerRHS.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 */
void SparseSolveWS_Float(SparseOpaqueFactorization_Float Factored, DenseVector_Float x,
  DenseVector_Float b, void *workspace);


/**** Symbolic Factorization Functions ****************************************/

/*! @abstract Returns a symbolic factorization of the requested type for a
 *            matrix with the given structure.
 *
 *  @discussion The resulting symbolic factorization may be used for multiple
 *  numerical factorizations with different numerical values but the same
 *  non-zero structure.
 *
 *  @param type The type of factorization to perform.
 *
 *  @param Matrix The structure of the sparse matrix to be factorized.
 *
 *  @returns The requested symbolic factorization of Matrix.                  */
SparseOpaqueSymbolicFactorization SparseFactorSymbolic(SparseFactorization_t type,
  SparseMatrixStructure Matrix);

/*! @abstract Returns a symbolic factorization of the requested type for a
 *            matrix with the given structure, with the supplied options.
 *
 *  @discussion The resulting symbolic factorization may be used for multiple
 *  numerical factorizations with different numerical values but the same
 *  non-zero structure.
 *
 *  @param type The type of factorization to perform.
 *
 *  @param Matrix The structure of the sparse matrix to be factorized.
 *
 *  @param sfoptions Symbolic factor options, for example the ordering algorithm
 *         to use.
 *
 *  @returns The requested symbolic factorization of Matrix.                  */
SparseOpaqueSymbolicFactorization SparseFactorSymbolicOpt(SparseFactorization_t type,
  SparseMatrixStructure Matrix, SparseSymbolicFactorOptions sfoptions);

/**** Symbolic Refactor Functions *********************************************/


/*! @abstract Reuses supplied factorization object's storage to compute a new
 *            factorization of the supplied matrix.
 *
 *  @discussion Matrix must have the same non-zero structure as that used for
 *  the original factorization.
 *  The same numerical factorization options will be used as in the original
 *  construction of Factorization.
 *  This call provides very similar behavior to that which can be achieved by
 *  reusing explicit storage supplied to SparseFactor() as the argument
 *  factorStorage. However, in addition to providing a simplified call sequence,
 *  this call can also reuse any additional storage allocated to accomodate
 *  delayed pivots.
 *  Note that if the reference count of the underlying object is not
 *  exactly one (i.e. if there are any implict copies as a result of calls to
 *  SparseGetTranspose() or SparseCreateSubfactor() that have not been destroyed
 *  through a call to SparseCleanup()), then new storage will be allocated
 *  regardless.
 *
 *  @param Matrix The matrix to be factorized.
 *
 *  @param Factorization The factorization to be updated.                     */
void SparseRefactor_Double(SparseMatrix_Double Matrix,
  SparseOpaqueFactorization_Double *Factorization);

/*! @abstract Reuses supplied factorization object's storage to compute a new
 *            factorization of the supplied matrix.
 *
 *  @discussion Matrix must have the same non-zero structure as that used for
 *  the original factorization.
 *  The same numerical factorization options will be used as in the original
 *  construction of Factorization.
 *  This call provides very similar behavior to that which can be achieved by
 *  reusing explicit storage supplied to SparseFactor() as the argument
 *  factorStorage. However, in addition to providing a simplified call sequence,
 *  this call can also reuse any additional storage allocated to accomodate
 *  delayed pivots.
 *  Note that if the reference count of the underlying object is not
 *  exactly one (i.e. if there are any implict copies as a result of calls to
 *  SparseGetTranspose() or SparseCreateSubfactor() that have not been destroyed
 *  through a call to SparseCleanup()), then new storage will be allocated
 *  regardless.
 *
 *  @param Matrix The matrix to be factorized.
 *
 *  @param Factorization The factorization to be updated.                     */
void SparseRefactor_Float(SparseMatrix_Float Matrix,
  SparseOpaqueFactorization_Float *Factorization);

/*! @abstract Reuses supplied factorization object's storage to compute a new
 *            factorization of the supplied matrix, using different options.
 *
 *  @discussion Matrix must have the same non-zero structure as that used for
 *  the original factorization.
 *  This call provides very similar behavior to that which can be achieved by
 *  reusing explicit storage supplied to SparseFactor() as the argument
 *  factorStorage. However, in addition to providing a simplified call sequence,
 *  this call can also reuse any additional storage allocated to accomodate
 *  delayed pivots.
 *  Note that if the reference count of the underlying object is not
 *  exactly one (i.e. if there are any implict copies as a result of calls to
 *  SparseGetTranspose() or SparseCreateSubfactor() that have not been destroyed
 *  through a call to SparseCleanup()), then new storage will be allocated
 *  regardless.
 *
 *  @param Matrix The matrix to be factorized.
 *
 *  @param Factorization The factorization to be updated.
 *
 *  @param nfoptions Numeric factor options, for example pivoting parameters. */
void SparseRefactorOpt_Double(SparseMatrix_Double Matrix,
  SparseOpaqueFactorization_Double *Factorization,
  SparseNumericFactorOptions nfoptions);

/*! @abstract Reuses supplied factorization object's storage to compute a new
 *            factorization of the supplied matrix, using different options.
 *
 *  @discussion Matrix must have the same non-zero structure as that used for
 *  the original factorization.
 *  This call provides very similar behavior to that which can be achieved by
 *  reusing explicit storage supplied to SparseFactor() as the argument
 *  factorStorage. However, in addition to providing a simplified call sequence,
 *  this call can also reuse any additional storage allocated to accomodate
 *  delayed pivots.
 *  Note that if the reference count of the underlying object is not
 *  exactly one (i.e. if there are any implict copies as a result of calls to
 *  SparseGetTranspose() or SparseCreateSubfactor() that have not been destroyed
 *  through a call to SparseCleanup()), then new storage will be allocated
 *  regardless.
 *
 *  @param Matrix The matrix to be factorized.
 *
 *  @param Factorization The factorization to be updated.
 *
 *  @param nfoptions Numeric factor options, for example pivoting parameters. */
void SparseRefactorOpt_Float(SparseMatrix_Float Matrix,
  SparseOpaqueFactorization_Float *Factorization,
  SparseNumericFactorOptions nfoptions);

/*! @abstract Reuses supplied factorization object's storage to compute a new
 *            factorization of the supplied matrix, without any internal
 *            allocations.
 *
 *  @discussion Matrix must have the same non-zero structure as that used for
 *  the original factorization.
 *  The same numerical factorization options will be used as in the original
 *  construction of Factorization.
 *  This call provides very similar behavior to that which can be achieved by
 *  reusing explicit storage supplied to SparseFactor() as the argument
 *  factorStorage. However, in addition to providing a simplified call sequence,
 *  this call can also reuse any additional storage allocated to accomodate
 *  delayed pivots.
 *  Note that internal memory allocations may occur in the case of
 *  pivoted factorizations that result in delayed pivots. If you require closer
 *  control over memory allocations, supply a sfoptions.malloc() function that
 *  implements the required behaviour, or use an alternative non-pivoted
 *  factorization returns. Note that if sfoptions.malloc() returns NULL the
 *  factorization will abort immediately.
 *  Note that if the reference count of the underlying object is not
 *  exactly one (i.e. if there are any implict copies as a result of calls to
 *  SparseGetTranspose() or SparseCreateSubfactor() that have not been destroyed
 *  through a call to SparseCleanup()), then new storage will be allocated
 *  regardless.
 *
 *  @param Matrix The matrix to be factorized.
 *
 *  @param Factorization The factorization to be updated.                     */
void SparseRefactorWS_Double(SparseMatrix_Double Matrix,
  SparseOpaqueFactorization_Double *Factorization, void *workspace);

/*! @abstract Reuses supplied factorization object's storage to compute a new
 *            factorization of the supplied matrix, without any internal
 *            allocations.
 *
 *  @discussion Matrix must have the same non-zero structure as that used for
 *  the original factorization.
 *  The same numerical factorization options will be used as in the original
 *  construction of Factorization.
 *  This call provides very similar behavior to that which can be achieved by
 *  reusing explicit storage supplied to SparseFactor() as the argument
 *  factorStorage. However, in addition to providing a simplified call sequence,
 *  this call can also reuse any additional storage allocated to accomodate
 *  delayed pivots.
 *  Note that internal memory allocations may occur in the case of
 *  pivoted factorizations that result in delayed pivots. If you require closer
 *  control over memory allocations, supply a sfoptions.malloc() function that
 *  implements the required behaviour, or use an alternative non-pivoted
 *  factorization returns. Note that if sfoptions.malloc() returns NULL the
 *  factorization will abort immediately.
 *  Note that if the reference count of the underlying object is not
 *  exactly one (i.e. if there are any implict copies as a result of calls to
 *  SparseGetTranspose() or SparseCreateSubfactor() that have not been destroyed
 *  through a call to SparseCleanup()), then new storage will be allocated
 *  regardless.
 *
 *  @param Matrix The matrix to be factorized.
 *
 *  @param Factorization The factorization to be updated.                     */
void SparseRefactorWS_Float(SparseMatrix_Float Matrix,
  SparseOpaqueFactorization_Float *Factorization, void *workspace);

/*! @abstract Reuses supplied factorization object's storage to compute a new
 *            factorization of the supplied matrix, using updated options and
 *            without any internal allocations.
 *
 *  @discussion Matrix must have the same non-zero structure as that used for
 *  the original factorization.

 *  This call provides very similar behavior to that which can be achieved by
 *  reusing explicit storage supplied to SparseFactor() as the argument
 *  factorStorage. However, in addition to providing a simplified call sequence,
 *  this call can also reuse any additional storage allocated to accomodate
 *  delayed pivots.
 *  Note that internal memory allocations may occur in the case of
 *  pivoted factorizations that result in delayed pivots. If you require closer
 *  control over memory allocations, supply a sfoptions.malloc() function that
 *  implements the required behaviour, or use an alternative non-pivoted
 *  factorization returns. Note that if sfoptions.malloc() returns NULL the
 *  factorization will abort immediately.
 *  Note that if the reference count of the underlying object is not
 *  exactly one (i.e. if there are any implict copies as a result of calls to
 *  SparseGetTranspose() or SparseCreateSubfactor() that have not been destroyed
 *  through a call to SparseCleanup()), then new storage will be allocated
 *  regardless.
 *
 *  @param Matrix The matrix to be factorized.
 *
 *  @param Factorization The factorization to be updated.
 *
 *  @param nfoptions Numeric factor options, for example pivoting parameters.
 *
 *  @param workspace A pointer to a workspace of size at least
 *         Factorization->symbolicFactorization.workspaceSize_Double bytes.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 *         This workspace may be reused or destroyed by the user as soon as the
 *         function returns.                                                  */
void SparseRefactorOptWS_Double(SparseMatrix_Double Matrix,
  SparseOpaqueFactorization_Double *Factorization,
  SparseNumericFactorOptions nfoptions, void *workspace);

/*! @abstract Reuses supplied factorization object's storage to compute a new
 *            factorization of the supplied matrix, using updated options and
 *            without any internal allocations.
 *
 *  @discussion Matrix must have the same non-zero structure as that used for
 *  the original factorization.

 *  This call provides very similar behavior to that which can be achieved by
 *  reusing explicit storage supplied to SparseFactor() as the argument
 *  factorStorage. However, in addition to providing a simplified call sequence,
 *  this call can also reuse any additional storage allocated to accomodate
 *  delayed pivots.
 *  Note that internal memory allocations may occur in the case of
 *  pivoted factorizations that result in delayed pivots. If you require closer
 *  control over memory allocations, supply a sfoptions.malloc() function that
 *  implements the required behaviour, or use an alternative non-pivoted
 *  factorization returns. Note that if sfoptions.malloc() returns NULL the
 *  factorization will abort immediately.
 *  Note that if the reference count of the underlying object is not
 *  exactly one (i.e. if there are any implict copies as a result of calls to
 *  SparseGetTranspose() or SparseCreateSubfactor() that have not been destroyed
 *  through a call to SparseCleanup()), then new storage will be allocated
 *  regardless.
 *
 *  @param Matrix The matrix to be factorized.
 *
 *  @param Factorization The factorization to be updated.
 *
 *  @param nfoptions Numeric factor options, for example pivoting parameters.
 *
 *  @param workspace A pointer to a workspace of size at least
 *         Factorization->symbolicFactorization.workspaceSize_Float bytes.
 *         This memory must be 16-byte aligned (any allocation returned
 *         by malloc() has this property).
 *         This workspace may be reused or destroyed by the user as soon as the
 *         function returns.                                                  */
void SparseRefactorOptWS_Float(SparseMatrix_Float Matrix,
  SparseOpaqueFactorization_Float *Factorization,
  SparseNumericFactorOptions nfoptions, void *workspace);

/**** Cleaning up resources ***************************************************/

void SparseCleanupOpaqueSymbolic(SparseOpaqueSymbolicFactorization toFree);
void SparseCleanupOpaqueNumeric_Double(SparseOpaqueFactorization_Double toFree);
void SparseCleanupOpaqueNumeric_Float(SparseOpaqueFactorization_Float toFree);
void SparseCleanupOpaqueSubfactor_Double(SparseOpaqueSubfactor_Double toFree);
void SparseCleanupOpaqueSubfactor_Float(SparseOpaqueSubfactor_Float toFree);
void SparseCleanupSparseMatrix_Double(SparseMatrix_Double toFree);
void SparseCleanupSparseMatrix_Float(SparseMatrix_Float toFree);
void SparseCleanupOpaquePreconditioner_Double(SparseOpaquePreconditioner_Double Preconditioner);
void SparseCleanupOpaquePreconditioner_Float(SparseOpaquePreconditioner_Float Preconditioner);

/**** Iterative Solvers *******************************************************/

/**** Iterative Method Factory Functions **************************************/

/*! @abstract Return a Conjugate Gradient Method with default options.        */
SparseIterativeMethod SparseConjugateGradientDefault(void);

/*! @abstract Return a Conjugate Gradient Method with specified options.
 *
 *  @param options Options for conjugate gradient method.                     */
SparseIterativeMethod SparseConjugateGradientOpt(SparseCGOptions options);

/*! @abstract Return a GMRES Method with default options.                     */
extern SparseIterativeMethod SparseGMRESDefault(void);

/*! @abstract Return a GMRES Method with specified options.
 *
 *  @param options Options for GMRES method.                                  */
SparseIterativeMethod SparseGMRESOpt(SparseGMRESOptions options);

/*! @abstract Return a LSMR Method with default options.                      */
SparseIterativeMethod SparseLSMRDefault(void);

/*! @abstract Return a LSMR Method with specified options
 *
 *  @param options Options for LSMR method.                                   */
SparseIterativeMethod SparseLSMROpt(SparseLSMROptions options);

/**** Solve without preconditioner ********************************************/

/*! @abstract Solve AX=B using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, and B has dimension m x nrhs, then X must have
 *         dimension n x nrhs. If no good initial estimate is available, user
 *         should set the initial guess to be the zero vector.                */
SparseIterativeStatus_t SparseSolveMatrixIterative_Double(SparseIterativeMethod method,
  SparseMatrix_Double A, DenseMatrix_Double B, DenseMatrix_Double X);

/*! @abstract Solve AX=B using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, and B has dimension m x nrhs, then X must have
 *         dimension n x nrhs. If no good initial estimate is available, user
 *         should set the initial guess to be the zero vector.                */
SparseIterativeStatus_t SparseSolveMatrixIterative_Float(SparseIterativeMethod method,
  SparseMatrix_Float A, DenseMatrix_Float B, DenseMatrix_Float X);

/*! @abstract Solve Ax=b using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param b The right-hand side b to solve for. If A has dimension m x n, then
 *         b must have length m.
 *
 *  @param x On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, then x must have length n. If no good initial
 *         estimate is available, user should set the initial guess to be the
 *         zero vector.                                                       */
SparseIterativeStatus_t SparseSolveIterative_Double(SparseIterativeMethod method,
  SparseMatrix_Double A, DenseVector_Double b, DenseVector_Double x);

/*! @abstract Solve Ax=b using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param b The right-hand side b to solve for. If A has dimension m x n, then
 *         b must have length m.
 *
 *  @param x On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, then x must have length n. If no good initial
 *         estimate is available, user should set the initial guess to be the
 *         zero vector.                                                       */
SparseIterativeStatus_t SparseSolveIterative_Float(SparseIterativeMethod method,
  SparseMatrix_Float A, DenseVector_Float b, DenseVector_Float x);

/*! @abstract Solve AX=B using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  \@callback ApplyOperator (block) ApplyOperator(accumulate, trans, X, Y)
 *             should perform the operation Y = op(A)X if accumulate is false,
 *             or Y += op(A)X if accumulate is true.
 *    \@param accumulate (input) Indicates whether to perform Y += op(A)X (if
 *            true) or Y = op(A)X (if false).
 *    \@param trans (input) Indicates whether op(A) is the application of A
 *            (trans=CblasNoTrans) or A^T (trans=CblasTrans).
 *    \@param X The matrix to multiply.
 *    \@param Y The matrix in which to accumulate or store the result.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, and B has dimension m x nrhs, then X must have
 *         dimension n x nrhs. If no good initial estimate is available, user
 *         should set the initial guess to be the zero vector.                */
SparseIterativeStatus_t SparseSolveMatrixIterativeOp_Double(SparseIterativeMethod method,
  void (^_Nonnull ApplyOperator)(bool accumulate, enum CBLAS_TRANSPOSE trans,
  DenseMatrix_Double X, DenseMatrix_Double Y),
  DenseMatrix_Double B, DenseMatrix_Double X);

/*! @abstract Solve AX=B using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  \@callback ApplyOperator (block) ApplyOperator(accumulate, trans, X, Y)
 *             should perform the operation Y = op(A)X if accumulate is false,
 *             or Y += op(A)X if accumulate is true.
 *    \@param accumulate (input) Indicates whether to perform Y += op(A)X (if
 *            true) or Y = op(A)X (if false).
 *    \@param trans (input) Indicates whether op(A) is the application of A
 *            (trans=CblasNoTrans) or A^T (trans=CblasTrans).
 *    \@param X The matrix to multiply.
 *    \@param Y The matrix in which to accumulate or store the result.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, and B has dimension m x nrhs, then X must have
 *         dimension n x nrhs. If no good initial estimate is available, user
 *         should set the initial guess to be the zero vector.                */
SparseIterativeStatus_t SparseSolveMatrixIterativeOp_Float(SparseIterativeMethod method,
  void (^_Nonnull ApplyOperator)(bool accumulate, enum CBLAS_TRANSPOSE trans,
  DenseMatrix_Float X, DenseMatrix_Float Y),
  DenseMatrix_Float B, DenseMatrix_Float X);

/*! @abstract Solve Ax=b using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  \@callback ApplyOperator (block) ApplyOperator(accumulate, trans, x, y)
 *             should perform the operation y = op(A)x if accumulate is false,
 *             or y += op(A)x if accumulate is true.
 *    \@param accumulate (input) Indicates whether to perform y += op(A)x (if
 *            true) or y = op(A)x (if false).
 *    \@param trans (input) Indicates whether op(A) is the application of A
 *            (trans=CblasNoTrans) or A^T (trans=CblasTrans).
 *    \@param x The vector to multiply.
 *    \@param y The vector in which to accumulate or store the result.
 *
 *  @param b The right-hand side b to solve for. If a has dimension m x n, then
 *         b must have length m.
 *
 *  @param x On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, then x must have length n. If no good initial
 *         estimate is available, user should set the initial guess to be the
 *         zero vector.                                                       */
SparseIterativeStatus_t SparseSolveIterativeOp_Double(SparseIterativeMethod method,
  void (^_Nonnull ApplyOperator)(bool accumulate, enum CBLAS_TRANSPOSE trans,
  DenseVector_Double x, DenseVector_Double y),
        DenseVector_Double b, DenseVector_Double x);

/*! @abstract Solve Ax=b using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  \@callback ApplyOperator (block) ApplyOperator(accumulate, trans, x, y)
 *             should perform the operation y = op(A)x if accumulate is false,
 *             or y += op(A)x if accumulate is true.
 *    \@param accumulate (input) Indicates whether to perform y += op(A)x (if
 *            true) or y = op(A)x (if false).
 *    \@param trans (input) Indicates whether op(A) is the application of A
 *            (trans=CblasNoTrans) or A^T (trans=CblasTrans).
 *    \@param x The vector to multiply.
 *    \@param y The vector in which to accumulate or store the result.
 *
 *  @param b The right-hand side b to solve for. If a has dimension m x n, then
 *         b must have length m.
 *
 *  @param x On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, then x must have length n. If no good initial
 *         estimate is available, user should set the initial guess to be the
 *         zero vector.                                                       */
SparseIterativeStatus_t SparseSolveIterativeOp_Float(SparseIterativeMethod method,
  void (^_Nonnull ApplyOperator)(bool accumulate, enum CBLAS_TRANSPOSE trans,
  DenseVector_Float x, DenseVector_Float y),
  DenseVector_Float b, DenseVector_Float x);

/**** Solve with preconditioner ***********************************************/

/*! @abstract Solve AX=B using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, and B has dimension m x nrhs, then X must have
 *         dimension n x nrhs. If no good initial estimate is available, user
 *         should set the initial guess to be the zero vector.
 *
 *  @param Preconditioner Type of preconditioner to create and apply.         */
SparseIterativeStatus_t SparseSolveMatrixIterativePrecond_Double(SparseIterativeMethod method,
  SparseMatrix_Double A, DenseMatrix_Double B, DenseMatrix_Double X,
  SparsePreconditioner_t Preconditioner);

/*! @abstract Solve AX=B using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, and B has dimension m x nrhs, then X must have
 *         dimension n x nrhs. If no good initial estimate is available, user
 *         should set the initial guess to be the zero vector.
 *
 *  @param Preconditioner Type of preconditioner to create and apply.         */
SparseIterativeStatus_t SparseSolveMatrixIterativePrecond_Float(SparseIterativeMethod method,
  SparseMatrix_Float A, DenseMatrix_Float B, DenseMatrix_Float X,
  SparsePreconditioner_t Preconditioner);

/*! @abstract Solve Ax=b using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param b The right-hand side b to solve for. If A has dimension m x n, then
 *         b must have length m.
 *
 *  @param x On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, then x must have length n. If no good initial
 *         estimate is available, user should set the initial guess to be the
 *         zero vector.
 *
 *  @param Preconditioner Type of preconditioner to create and apply.         */
SparseIterativeStatus_t SparseSolveIterativePrecond_Double(SparseIterativeMethod method,
  SparseMatrix_Double A, DenseVector_Double b, DenseVector_Double x,
  SparsePreconditioner_t Preconditioner);

/*! @abstract Solve Ax=b using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param b The right-hand side b to solve for. If A has dimension m x n, then
 *         b must have length m.
 *
 *  @param x On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, then x must have length n. If no good initial
 *         estimate is available, user should set the initial guess to be the
 *         zero vector.
 *
 *  @param Preconditioner Type of preconditioner to create and apply.         */
SparseIterativeStatus_t SparseSolveIterativePrecond_Float(SparseIterativeMethod method,
  SparseMatrix_Float A, DenseVector_Float b, DenseVector_Float x,
  SparsePreconditioner_t Preconditioner);

/*! @abstract Solve AX=B using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, and B has dimension m x nrhs, then X must have
 *         dimension n x nrhs. If no good initial estimate is available, user
 *         should set the initial guess to be the zero vector.
 *
 *  @parameter Preconditioner The preconditioner to apply.                    */
SparseIterativeStatus_t SparseSolveMatrixIterativeOpaquePrecond_Double(SparseIterativeMethod method,
  SparseMatrix_Double A, DenseMatrix_Double B, DenseMatrix_Double X,
  SparseOpaquePreconditioner_Double Preconditioner);

/*! @abstract Solve AX=B using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, and B has dimension m x nrhs, then X must have
 *         dimension n x nrhs. If no good initial estimate is available, user
 *         should set the initial guess to be the zero vector.
 *
 *  @parameter Preconditioner The preconditioner to apply.                    */
SparseIterativeStatus_t SparseSolveMatrixIterativeOpaquePrecond_Float(SparseIterativeMethod method,
  SparseMatrix_Float A, DenseMatrix_Float B, DenseMatrix_Float X,
  SparseOpaquePreconditioner_Float Preconditioner);

/*! @abstract Solve Ax=b using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param b The right-hand side b to solve for. If A has dimension m x n, then
 *         b must have length m.
 *
 *  @param x On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, then x must have length n. If no good initial
 *         estimate is available, user should set the initial guess to be the
 *         zero vector.
 *
 *  @param Preconditioner The preconditioner to apply.                        */
SparseIterativeStatus_t SparseSolveIterativeOpaquePrecond_Double(SparseIterativeMethod method,
  SparseMatrix_Double A, DenseVector_Double b, DenseVector_Double x,
  SparseOpaquePreconditioner_Double Preconditioner);

/*! @abstract Solve Ax=b using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  @param A (input) The matrix A to solve the system for. Only used for
 *         multiplication by A or A^T.
 *
 *  @param b The right-hand side b to solve for. If A has dimension m x n, then
 *         b must have length m.
 *
 *  @param x On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, then x must have length n. If no good initial
 *         estimate is available, user should set the initial guess to be the
 *         zero vector.
 *
 *  @param Preconditioner The preconditioner to apply.                        */
SparseIterativeStatus_t SparseSolveIterativeOpaquePrecond_Float(SparseIterativeMethod method,
  SparseMatrix_Float A, DenseVector_Float b, DenseVector_Float x,
  SparseOpaquePreconditioner_Float Preconditioner);

/*! @abstract Solve AX=B using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  \@callback ApplyOperator (block) ApplyOperator(accumulate, trans, X, Y)
 *             should perform the operation Y = op(A)X if accumulate is false,
 *             or Y += op(A)X if accumulate is true.
 *    \@param accumulate (input) Indicates whether to perform Y += op(A)X (if
 *            true) or Y = op(A)X (if false).
 *    \@param trans (input) Indicates whether op(A) is the application of A
 *            (trans=CblasNoTrans) or A^T (trans=CblasTrans).
 *    \@param X The matrix to multiply.
 *    \@param Y The matrix in which to accumulate or store the result.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, and B has dimension m x nrhs, then X must have
 *         dimension n x nrhs. If no good initial estimate is available, user
 *         should set the initial guess to be the zero vector.
 *
 *  @param Preconditioner (input) The preconditioner to use.                  */
SparseIterativeStatus_t SparseSolveMatrixIterativeOpOpaquePrecond_Double(SparseIterativeMethod method,
  void (^_Nonnull ApplyOperator)(bool accumulate, enum CBLAS_TRANSPOSE trans,
  DenseMatrix_Double X, DenseMatrix_Double Y),
  DenseMatrix_Double B, DenseMatrix_Double X,
  SparseOpaquePreconditioner_Double Preconditioner);

/*! @abstract Solve AX=B using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  \@callback ApplyOperator (block) ApplyOperator(accumulate, trans, X, Y)
 *             should perform the operation Y = op(A)X if accumulate is false,
 *             or Y += op(A)X if accumulate is true.
 *    \@param accumulate (input) Indicates whether to perform Y += op(A)X (if
 *            true) or Y = op(A)X (if false).
 *    \@param trans (input) Indicates whether op(A) is the application of A
 *            (trans=CblasNoTrans) or A^T (trans=CblasTrans).
 *    \@param X The matrix to multiply.
 *    \@param Y The matrix in which to accumulate or store the result.
 *
 *  @param B The right-hand sides B to solve for. If A has dimension m x n, then
 *         B must have dimension m x nrhs, where nrhs is the number of
 *         right-hand sides to find solutions for.
 *
 *  @param X On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, and B has dimension m x nrhs, then X must have
 *         dimension n x nrhs. If no good initial estimate is available, user
 *         should set the initial guess to be the zero vector.
 *
 *  @param Preconditioner (input) The preconditioner to use.                  */
SparseIterativeStatus_t SparseSolveMatrixIterativeOpOpaquePrecond_Float(SparseIterativeMethod method,
  void (^_Nonnull ApplyOperator)(bool accumulate, enum CBLAS_TRANSPOSE trans,
  DenseMatrix_Float X, DenseMatrix_Float Y),
  DenseMatrix_Float B, DenseMatrix_Float X,
  SparseOpaquePreconditioner_Float Preconditioner);

/*! @abstract Solve Ax=b using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  \@callback ApplyOperator (block) ApplyOperator(accumulate, trans, x, y)
 *             should perform the operation y = op(A)x if accumulate is false,
 *             or y += op(A)x if accumulate is true.
 *    \@param accumulate (input) Indicates whether to perform y += op(A)x (if
 *            true) or y = op(A)x (if false).
 *    \@param trans (input) Indicates whether op(A) is the application of A
 *            (trans=CblasNoTrans) or A^T (trans=CblasTrans).
 *    \@param x The vector to multiply.
 *    \@param y The vector in which to accumulate or store the result.
 *
 *  @param b The right-hand side b to solve for. If a has dimension m x n, then
 *         b must have length m.
 *
 *  @param x On entry, initial guess for solution, on return the solution. If A
 *         has dimension m x n, then x must have length n. If no good initial
 *         estimate is available, user should set the initial guess to be the
 *         zero vector.
 *
 *  @param Preconditioner (input) The preconditioner to use.                  */
SparseIterativeStatus_t SparseSolveIterativeOpOpaquePrecond_Double(SparseIterativeMethod method,
  void (^_Nonnull ApplyOperator)(bool accumulate, enum CBLAS_TRANSPOSE trans,
  DenseVector_Double x, DenseVector_Double y),
  DenseVector_Double b, DenseVector_Double x,
  SparseOpaquePreconditioner_Double Preconditioner);

/*! @abstract Solve Ax=b using the specified iterative method.
 *
 *  @param method (input) Iterative method specification, eg return value of
 *         SparseConjugateGradient().
 *
 *  \@callback ApplyOperator (block) ApplyOperator(accumulate, trans, x, y)
 *             should perform the operation y = op(A)x if accumulate is false,
 *             or y += op(A)x if accumulate is true.
 *    \@param accumulate (input) Indicates whether to perform y += op(A)x (if
 *            true) or y = op(A)x (if false).
 *    \@param trans (input) Indicates whether op(A) is the application of A
 *            (trans=CblasNoTrans) or A^T (trans=CblasTrans).
 *    \@param x The vector to multiply.
 *    \@param y The vector in which to accumulate or store the result.
 *
 *  \@param b The right-hand side b to solve for. If a has dimension m x n, then
 *          b must have length m.
 *
 *  \@param x On entry, initial guess for solution, on return the solution. If A
 *          has dimension m x n, then x must have length n. If no good initial
 *          estimate is available, user should set the initial guess to be the
 *          zero vector.
 *
 *  @param Preconditioner (input) The preconditioner to use.                  */
SparseIterativeStatus_t SparseSolveIterativeOpOpaquePrecond_Float(SparseIterativeMethod method,
  void (^_Nonnull ApplyOperator)(bool accumulate, enum CBLAS_TRANSPOSE trans,
  DenseVector_Float x, DenseVector_Float y),
  DenseVector_Float b, DenseVector_Float x,
  SparseOpaquePreconditioner_Float Preconditioner);

/******************************************************************************
 *  @group Matrix and Vector Operations (Sparse BLAS Wrappers)
 ******************************************************************************/

/**** Multiplication **********************************************************/

/*! @abstract Performs the multiplication Y = AX for double values
 *
 *  @param A (input) sparse matrix.
 *
 *  @param X (input) dense matrix. Inner dimensions of A and X must match.
 *
 *  @param Y (output) dense matrix. Dimensions must match the outer dimensions
 *  of A and X. Overwritten with their product.                               */
void SparseMultiplyMatrix_Double(SparseMatrix_Double A, DenseMatrix_Double X, DenseMatrix_Double Y);

/*! @abstract Performs the multiplication Y = AX for float values.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param X (input) dense matrix. Inner dimensions of A and X must match.
 *
 *  @param Y (output) dense matrix. Dimensions must match the outer dimensions
 *  of A and X. Overwritten with their product.                               */
void SparseMultiplyMatrix_Float(SparseMatrix_Float A, DenseMatrix_Float X, DenseMatrix_Float Y);

/*! @abstract Performs the multiplication Y = alpha * AX for double values
 *
 *  @param alpha (input) scale to apply to the result.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param X (input) dense matrix. Inner dimensions of A and X must match.
 *
 *  @param Y (output) dense matrix. Dimensions must match the outer dimensions
 *  of A and X. Overwritten with alpha * AX.                                  */
void SparseMultiplyScalarMatrix_Double(double alpha, SparseMatrix_Double A, DenseMatrix_Double X, DenseMatrix_Double Y);

/*! @abstract Performs the multiplication Y = alpha * AX for float values.
 *
 *  @param alpha (input) scale to apply to the result.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param X (input) dense matrix. Inner dimensions of A and X must match.
 *
 *  @param Y (output) dense matrix. Dimensions must match the outer dimensions
 *  of A and X. Overwritten with alpha * AX.                                  */
void SparseMultiplyScalarMatrix_Float(float alpha, SparseMatrix_Float A, DenseMatrix_Float X, DenseMatrix_Float Y);

/*! @abstract Performs the multiplication y = Ax for double values
 *
 *  @param A (input) sparse matrix.
 *
 *  @param x (input) dense vector.
 *
 *  @param y (output) dense vector.                                           */
void SparseMultiply_Double(SparseMatrix_Double A, DenseVector_Double x, DenseVector_Double y);

/*! @abstract Performs the multiplication y = Ax for float values
 *
 *  @param A (input) sparse matrix.
 *
 *  @param x (input) dense vector.
 *
 *  @param y (output) dense vector.                                           */
void SparseMultiply_Float(SparseMatrix_Float A, DenseVector_Float x, DenseVector_Float y);

/*! @abstract Performs the multiplication y = alpha * Ax for double values
 *
 *  @param alpha (input) scale to apply to the result.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param x (input) dense vector.
 *
 *  @param y (output) dense vector.                                           */
void SparseMultiplyScalar_Double(double alpha, SparseMatrix_Double A, DenseVector_Double x, DenseVector_Double y);

/*! @abstract Performs the multiplication y = alpha * Ax for float values.
 *
 *  @param alpha (input) scale to apply to the result.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param x (input) dense vector.
 *
 *  @param y (output) dense vector.                                           */
void SparseMultiplyScalar_Float(float alpha, SparseMatrix_Float A, DenseVector_Float x, DenseVector_Float y);

/**** Multipy-Add *************************************************************/

/*! @abstract Y += AX for double values
 *
 *  @param A (input) sparse matrix.
 *
 *  @param X (input) dense matrix. Inner dimensions of A and X must match.
 *
 *  @param Y (output) dense matrix. Dimensions must match the outer dimensions
 *  of A and X. Overwritten with their product.                               */
void SparseMultiplyAddMatrix_Double(SparseMatrix_Double A, DenseMatrix_Double X, DenseMatrix_Double Y);

/*! @abstract Y += AX for float values.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param X (input) dense matrix. Inner dimensions of A and X must match.
 *
 *  @param Y (output) dense matrix. Dimensions must match the outer dimensions
 *  of A and X. Overwritten with their product.                               */
void SparseMultiplyAddMatrix_Float(SparseMatrix_Float A, DenseMatrix_Float X, DenseMatrix_Float Y);

/*! @abstract Y += alpha * AX for double values
 *
 *  @param alpha (input) scale to apply to the product of A and X.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param X (input) dense matrix. Inner dimensions of A and X must match.
 *
 *  @param Y (output) dense matrix. Dimensions must match the outer dimensions
 *  of A and X. Overwritten with alpha * AX.                                  */
void SparseMultiplyAddScalarMatrix_Double(double alpha, SparseMatrix_Double A, DenseMatrix_Double X, DenseMatrix_Double Y);

/*! @abstract Y += alpha * AX for float values.
 *
 *  @param alpha (input) scale to apply to the product of A and X.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param X (input) dense matrix. Inner dimensions of A and X must match.
 *
 *  @param Y (output) dense matrix. Dimensions must match the outer dimensions
 *  of A and X. Overwritten with alpha * AX.                                  */
void SparseMultiplyAddScalarMatrix_Float(float alpha, SparseMatrix_Float A, DenseMatrix_Float X, DenseMatrix_Float Y);

/*! @abstract y += Ax for double values
 *
 *  @param A (input) sparse matrix.
 *
 *  @param x (input) dense vector.
 *
 *  @param y (output) dense vector.                                           */
void SparseMultiplyAdd_Double(SparseMatrix_Double A, DenseVector_Double x, DenseVector_Double y);

/*! @abstract y += Ax for float values.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param x (input) dense vector.
 *
 *  @param y (output) dense vector.                                           */
void SparseMultiplyAdd_Float(SparseMatrix_Float A, DenseVector_Float x, DenseVector_Float y);

/*! @abstract y += alpha * Ax for double values
 *
 *  @param alpha (input) scale to apply to the product of A and x.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param x (input) dense vector.
 *
 *  @param y (output) dense vector.                                           */
void SparseMultiplyAddScalar_Double(double alpha, SparseMatrix_Double A, DenseVector_Double x, DenseVector_Double y);

/*! @abstract y += alpha * Ax for float values.
 *
 *  @param alpha (input) scale to apply to the product of A and x.
 *
 *  @param A (input) sparse matrix.
 *
 *  @param x (input) dense vector.
 *
 *  @param y (output) dense vector.                                           */
void SparseMultiplyAddScalar_Float(float alpha, SparseMatrix_Float A, DenseVector_Float x, DenseVector_Float y);


#ifdef __cplusplus
} // extern "C"
#endif

#endif //__ACCELERATE_SYS_C_API_H__
