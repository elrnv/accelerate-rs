use std::convert::TryFrom;
use std::marker::PhantomData;

use accelerate_sys as ffi;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SparseTriangle {
    Upper,
    Lower,
}

impl From<ffi::SparseTriangle_t> for SparseTriangle {
    fn from(st: ffi::SparseTriangle_t) -> Self {
        match st as u32 {
            ffi::SparseTriangle_SparseUpperTriangle => SparseTriangle::Upper,
            ffi::SparseTriangle_SparseLowerTriangle => SparseTriangle::Lower,
            _ => panic!("Invalid triangle type"),
        }
    }
}

impl From<SparseTriangle> for ffi::SparseTriangle_t {
    fn from(st: SparseTriangle) -> Self {
        match st {
            SparseTriangle::Upper => {
                ffi::SparseTriangle_SparseUpperTriangle as ffi::SparseTriangle_t
            }
            SparseTriangle::Lower => {
                ffi::SparseTriangle_SparseLowerTriangle as ffi::SparseTriangle_t
            }
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SparseKind {
    Ordinary,
    Triangular,
    UnitTriangular,
    Symmetric,
}

impl From<ffi::SparseKind_t> for SparseKind {
    fn from(st: ffi::SparseKind_t) -> Self {
        match st {
            ffi::SparseOrdinary => SparseKind::Ordinary,
            ffi::SparseTriangular => SparseKind::Triangular,
            ffi::SparseUnitTriangular => SparseKind::UnitTriangular,
            ffi::SparseSymmetric => SparseKind::Symmetric,
            _ => panic!("Invalid sparse kind"),
        }
    }
}

impl From<SparseKind> for ffi::SparseKind_t {
    fn from(st: SparseKind) -> Self {
        match st {
            SparseKind::Ordinary => ffi::SparseOrdinary,
            SparseKind::Triangular => ffi::SparseTriangular,
            SparseKind::UnitTriangular => ffi::SparseUnitTriangular,
            SparseKind::Symmetric => ffi::SparseSymmetric,
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct SparseAttributes {
    attrs: ffi::SparseAttributes_t,
}

impl From<ffi::SparseAttributes_t> for SparseAttributes {
    fn from(attrs: ffi::SparseAttributes_t) -> Self {
        Self { attrs }
    }
}

impl From<SparseAttributes> for ffi::SparseAttributes_t {
    fn from(attrs: SparseAttributes) -> Self {
        attrs.attrs
    }
}

impl SparseAttributes {
    /// Constructs default attributes.
    ///
    /// The default corresponds to non-transposed, ordinary sparse matrix.
    pub fn new() -> Self {
        SparseAttributes {
            attrs: ffi::SparseAttributes_t {
                _bitfield_align_1: [0; 0],
                _bitfield_1: ffi::SparseAttributes_t::new_bitfield_1(
                    false,
                    SparseTriangle::Lower.into(),
                    ffi::SparseOrdinary,
                    0,
                    true, // Will cause memory leak if not handled carefully.
                ),
                __bindgen_padding_0: 0,
            },
        }
    }
    pub fn transposed(mut self) -> Self {
        self.attrs.set_transpose(true);
        self
    }
    pub fn with_kind(mut self, kind: SparseKind) -> Self {
        self.attrs.set_kind(kind.into());
        self
    }
    pub fn with_triangle(mut self, triangle: SparseTriangle) -> Self {
        self.attrs.set_triangle(triangle.into());
        self
    }
    /// Tells Accelerate that the matrix is allocated by the user.
    pub(crate) fn with_owned(mut self, owned: bool) -> Self {
        self.attrs.set__allocatedBySparse(owned);
        self
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SparseFactorizationType {
    Cholesky,
    LDLT,
    LDLTUnpivoted,
    LDLTSBK,
    LDLTTPP,
    QR,
    CholeskyAtA,
}

impl From<ffi::SparseFactorization_t> for SparseFactorizationType {
    fn from(st: ffi::SparseFactorization_t) -> Self {
        match st as u32 {
            ffi::SparseFactorizationCholesky => SparseFactorizationType::Cholesky,
            ffi::SparseFactorizationCholeskyAtA => SparseFactorizationType::CholeskyAtA,
            ffi::SparseFactorizationLDLT => SparseFactorizationType::LDLT,
            ffi::SparseFactorizationLDLTUnpivoted => SparseFactorizationType::LDLTUnpivoted,
            ffi::SparseFactorizationLDLTSBK => SparseFactorizationType::LDLTSBK,
            ffi::SparseFactorizationLDLTTPP => SparseFactorizationType::LDLTTPP,
            ffi::SparseFactorizationQR => SparseFactorizationType::QR,
            _ => panic!("Invalid factorization type"),
        }
    }
}

impl From<SparseFactorizationType> for ffi::SparseFactorization_t {
    fn from(st: SparseFactorizationType) -> Self {
        match st {
            SparseFactorizationType::Cholesky => {
                ffi::SparseFactorizationCholesky as ffi::SparseFactorization_t
            }
            SparseFactorizationType::CholeskyAtA => {
                ffi::SparseFactorizationCholeskyAtA as ffi::SparseFactorization_t
            }
            SparseFactorizationType::LDLT => {
                ffi::SparseFactorizationLDLT as ffi::SparseFactorization_t
            }
            SparseFactorizationType::LDLTUnpivoted => {
                ffi::SparseFactorizationLDLTUnpivoted as ffi::SparseFactorization_t
            }
            SparseFactorizationType::LDLTSBK => {
                ffi::SparseFactorizationLDLTSBK as ffi::SparseFactorization_t
            }
            SparseFactorizationType::LDLTTPP => {
                ffi::SparseFactorizationLDLTTPP as ffi::SparseFactorization_t
            }
            SparseFactorizationType::QR => ffi::SparseFactorizationQR as ffi::SparseFactorization_t,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SparseOrder {
    Default,
    User,
    AMD,
    Metis,
    COLAMD,
}
impl Default for SparseOrder {
    fn default() -> Self {
        SparseOrder::Default
    }
}
impl From<ffi::SparseOrder_t> for SparseOrder {
    fn from(st: ffi::SparseOrder_t) -> Self {
        match st as u32 {
            ffi::SparseOrderDefault => SparseOrder::Default,
            ffi::SparseOrderUser => SparseOrder::User,
            ffi::SparseOrderAMD => SparseOrder::AMD,
            ffi::SparseOrderMetis => SparseOrder::Metis,
            ffi::SparseOrderCOLAMD => SparseOrder::COLAMD,
            _ => panic!("Invalid sparse order"),
        }
    }
}

impl From<SparseOrder> for ffi::SparseOrder_t {
    fn from(st: SparseOrder) -> Self {
        match st {
            SparseOrder::Default => ffi::SparseOrderDefault as ffi::SparseOrder_t,
            SparseOrder::User => ffi::SparseOrderUser as ffi::SparseOrder_t,
            SparseOrder::AMD => ffi::SparseOrderAMD as ffi::SparseOrder_t,
            SparseOrder::Metis => ffi::SparseOrderMetis as ffi::SparseOrder_t,
            SparseOrder::COLAMD => ffi::SparseOrderCOLAMD as ffi::SparseOrder_t,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SparseScaling {
    Default,
    User,
    EquilibriationInf,
}

impl Default for SparseScaling {
    fn default() -> Self {
        SparseScaling::Default
    }
}
impl From<ffi::SparseScaling_t> for SparseScaling {
    fn from(st: ffi::SparseScaling_t) -> Self {
        match st as u32 {
            ffi::SparseScalingDefault => SparseScaling::Default,
            ffi::SparseScalingUser => SparseScaling::User,
            ffi::SparseScalingEquilibriationInf => SparseScaling::EquilibriationInf,
            _ => panic!("Invalid scaling type"),
        }
    }
}

impl From<SparseScaling> for ffi::SparseScaling_t {
    fn from(st: SparseScaling) -> Self {
        match st {
            SparseScaling::Default => ffi::SparseScalingDefault as ffi::SparseScaling_t,
            SparseScaling::User => ffi::SparseScalingUser as ffi::SparseScaling_t,
            SparseScaling::EquilibriationInf => {
                ffi::SparseScalingEquilibriationInf as ffi::SparseScaling_t
            }
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SparseStatus {
    /// Factorization was successful.
    Ok,
    /// Factorization failed due to a numerical issue.
    FactorizationFailed,
    /// Factorization aborted as matrix is singular.
    MatrixIsSingular,
    /// Factorization encountered an internal error (e.g. failed to allocate memory).
    InternalError,
    /// Error in user supplied parameter.
    ParameterError,
    /// Factorization object has been freed.
    Released,
}

impl From<ffi::SparseStatus_t> for SparseStatus {
    fn from(st: ffi::SparseStatus_t) -> Self {
        match st {
            ffi::SparseStatusOK => SparseStatus::Ok,
            ffi::SparseFactorizationFailed => SparseStatus::FactorizationFailed,
            ffi::SparseMatrixIsSingular => SparseStatus::MatrixIsSingular,
            ffi::SparseInternalError => SparseStatus::InternalError,
            ffi::SparseParameterError => SparseStatus::ParameterError,
            ffi::SparseStatusReleased => SparseStatus::Released,
            _ => panic!("Invalid status"),
        }
    }
}

impl From<SparseStatus> for ffi::SparseStatus_t {
    fn from(st: SparseStatus) -> Self {
        match st {
            SparseStatus::Ok => ffi::SparseStatusOK,
            SparseStatus::FactorizationFailed => ffi::SparseFactorizationFailed,
            SparseStatus::MatrixIsSingular => ffi::SparseMatrixIsSingular,
            SparseStatus::InternalError => ffi::SparseInternalError,
            SparseStatus::ParameterError => ffi::SparseParameterError,
            SparseStatus::Released => ffi::SparseStatusReleased,
        }
    }
}

/// Trait defining scalar types supported by the Accelerate sparse module.
pub trait SupportedScalar {
    type FFISparseMatrixType: Copy;
    type FFISparseFactorizationType: Copy;
    unsafe fn cleanup_matrix(mtx: Self::FFISparseMatrixType);
    unsafe fn cleanup_factorization(fact: Self::FFISparseFactorizationType);
}

impl SupportedScalar for f32 {
    type FFISparseMatrixType = ffi::SparseMatrix_Float;
    type FFISparseFactorizationType = ffi::SparseOpaqueFactorization_Float;
    unsafe fn cleanup_matrix(mtx: Self::FFISparseMatrixType) {
        if mtx.structure.attributes._allocatedBySparse() {
            ffi::SparseCleanupSparseMatrix_Float(mtx);
        }
    }
    unsafe fn cleanup_factorization(fact: Self::FFISparseFactorizationType) {
        ffi::SparseCleanupOpaqueNumeric_Float(fact);
    }
}

impl SupportedScalar for f64 {
    type FFISparseMatrixType = ffi::SparseMatrix_Double;
    type FFISparseFactorizationType = ffi::SparseOpaqueFactorization_Double;
    unsafe fn cleanup_matrix(mtx: Self::FFISparseMatrixType) {
        if mtx.structure.attributes._allocatedBySparse() {
            ffi::SparseCleanupSparseMatrix_Double(mtx);
        }
    }
    unsafe fn cleanup_factorization(fact: Self::FFISparseFactorizationType) {
        ffi::SparseCleanupOpaqueNumeric_Double(fact);
    }
}

/// A block CSC format sparse matrix.
#[derive(Debug)]
pub struct SparseMatrix<'a, T: SupportedScalar> {
    mtx: T::FFISparseMatrixType,
    phantom: PhantomData<&'a ()>,
}

// `SparseMatrix` cannot be copied or cloned.
unsafe impl<'a, T: SupportedScalar> Send for SparseMatrix<'a, T> {}

impl<'a, T: SupportedScalar> Drop for SparseMatrix<'a, T> {
    fn drop(&mut self) {
        // SAFETY: cleanup occurs only if the matrix was allocated internally.
        unsafe {
            T::cleanup_matrix(self.mtx);
        }
    }
}

pub type SparseMatrixF32<'a> = SparseMatrix<'a, f32>;
pub type SparseMatrixF64<'a> = SparseMatrix<'a, f64>;

/// A block CSC format sparse matrix.
#[derive(Debug)]
pub struct SparseFactorization<T: SupportedScalar> {
    fact: T::FFISparseFactorizationType,
}

impl<T: SupportedScalar> Drop for SparseFactorization<T> {
    fn drop(&mut self) {
        unsafe {
            T::cleanup_factorization(self.fact);
        }
    }
}

/// Factorization cannot be copied or cloned.
unsafe impl<T: SupportedScalar> Send for SparseFactorization<T> {}

impl From<ffi::SparseOpaqueFactorization_Float> for SparseFactorization<f32> {
    fn from(f: ffi::SparseOpaqueFactorization_Float) -> Self {
        Self { fact: f }
    }
}
impl From<ffi::SparseOpaqueFactorization_Double> for SparseFactorization<f64> {
    fn from(f: ffi::SparseOpaqueFactorization_Double) -> Self {
        Self { fact: f }
    }
}

impl From<SparseFactorization<f32>> for ffi::SparseOpaqueFactorization_Float {
    fn from(f: SparseFactorization<f32>) -> ffi::SparseOpaqueFactorization_Float {
        f.fact
    }
}

impl From<SparseFactorization<f64>> for ffi::SparseOpaqueFactorization_Double {
    fn from(f: SparseFactorization<f64>) -> ffi::SparseOpaqueFactorization_Double {
        f.fact
    }
}

macro_rules! impl_matrix {
    ($t:ident, $mtx:ident, $convert:ident, $factor:ident, $factor_numeric:ident, $factor_numeric_opt:ident, $cleanup:ident $(,)?) => {
        impl SparseMatrix<'static, $t> {
            /// Constructs a sparse matrix from data arrays.
            ///
            /// Out-of-range entries are dropped and duplicates are summed.
            pub fn from_coordinate(
                num_rows: i32,
                num_cols: i32,
                block_size: u8,
                attributes: SparseAttributes,
                rows: &[i32],
                cols: &[i32],
                values: &[$t],
            ) -> Self {
                let block_count = i64::try_from(values.len() / block_size as usize).unwrap();

                Self {
                    mtx: unsafe {
                        ffi::$convert(
                            num_rows,
                            num_cols,
                            block_count,
                            block_size,
                            attributes.with_owned(true).into(),
                            rows.as_ptr(),
                            cols.as_ptr(),
                            values.as_ptr(),
                        )
                    },
                    phantom: PhantomData,
                }
            }
        }

        impl<'a> SparseMatrix<'a, $t> {
            /// Constructs a sparse BCSC matrix from the given raw indices, offsets and values.
            pub fn from_raw_parts(
                num_rows: usize,
                num_cols: usize,
                block_size: usize,
                attributes: SparseAttributes,
                indices: &mut [i32],
                offsets: &mut [i64],
                values: &mut [$t],
            ) -> Self {
                Self {
                    mtx: ffi::$mtx {
                        // The attributes will be changed to reflect that this matrix structure and
                        // values are owned by the caller.
                        structure: SparseMatrixStructure::from_raw_parts(
                            num_rows,
                            num_cols,
                            block_size,
                            attributes.with_owned(false),
                            indices,
                            offsets,
                        )
                        .into(),
                        data: values.as_mut_ptr(),
                    },
                    phantom: PhantomData,
                }
            }
        }

        impl<'a> SparseMatrix<'a, $t> {
            /// A mutable slice of all structurally non-zero entries in this matrix.
            ///
            /// This can be useful for updating the values of the matrix without changing the
            /// sparsity pattern.
            pub fn data_mut(&mut self) -> &mut [$t] {
                let num_cols = usize::try_from(self.mtx.structure.columnCount).unwrap();
                unsafe {
                    std::slice::from_raw_parts_mut(
                        self.mtx.data,
                        usize::try_from(*self.mtx.structure.columnStarts.add(num_cols)).unwrap(),
                    )
                }
            }
            /// A slice of all structurally non-zero entries in this matrix.
            pub fn data(&self) -> &[$t] {
                let num_cols = usize::try_from(self.mtx.structure.columnCount).unwrap();
                unsafe {
                    std::slice::from_raw_parts(
                        self.mtx.data,
                        usize::try_from(*self.mtx.structure.columnStarts.add(num_cols)).unwrap(),
                    )
                }
            }
            /// A slice of all row indices represented in this BCSC matrix.
            ///
            /// Each index corresponds to an structurally non-zero entry in the matrix.
            pub fn indices(&self) -> &[i32] {
                sparse_matrix_structure_indices(&self.mtx.structure)
            }
            /// A slice of all column offsets represented in this BCSC matrix.
            ///
            /// Each offset (except for one) corresponds to a column in the matrix. The last
            /// offset corresponds to the total number of structural non-zeros in the
            /// matrix.
            pub fn offsets(&self) -> &[i64] {
                sparse_matrix_structure_offsets(&self.mtx.structure)
            }
            /// Factor this matrix according to the given factorization type.
            pub fn factor(
                &self,
                factorization_type: SparseFactorizationType,
            ) -> SparseFactorization<$t> {
                unsafe { ffi::$factor(factorization_type.into(), self.mtx) }.into()
            }
            /// Factor this matrix using the given symbolic factorization.
            pub fn factor_with(
                &self,
                factorization: &SparseSymbolicFactorization,
            ) -> SparseFactorization<$t> {
                unsafe { ffi::$factor_numeric(factorization.factorization.clone(), self.mtx) }
                    .into()
            }
            /// Factor this matrix using the given symbolic factorization and options.
            pub fn factor_with_options(
                &self,
                factorization: &SparseSymbolicFactorization,
                options: SparseNumericFactorOptions<'_, $t>,
            ) -> SparseFactorization<$t> {
                unsafe {
                    ffi::$factor_numeric_opt(
                        factorization.factorization.clone(),
                        self.mtx,
                        options.into(),
                    )
                }
                .into()
            }
            /// Returns a copy of this matrix structure without entries.
            pub fn structure(&self) -> SparseMatrixStructure {
                self.mtx.structure.into()
            }
        }
    };
}

macro_rules! impl_factorization {
    ($t:ident, $solve:ident, $solve_in_place:ident, $dense_vec:ident, $dense_mtx:ident,
     $cleanup:ident, $refactor:ident, $refactor_opt:ident $(,)?) => {
        impl SparseFactorization<$t> {
            /// Returns the status of this factorization state.
            pub fn status(&self) -> SparseStatus {
                self.fact.status.into()
            }
            /// Refactor the given matrix.
            ///
            /// The given matrix must have the exact same sparsity pattern as originally provided.
            pub fn refactor<'a>(&mut self, mtx: &SparseMatrix<'a, $t>) -> SparseStatus {
                unsafe {
                    ffi::$refactor(
                        mtx.mtx,
                        &mut self.fact as *mut <$t as SupportedScalar>::FFISparseFactorizationType,
                    );
                }
                self.status()
            }
            /// Refactor the given matrix with the given options.
            ///
            /// The given matrix must have the exact same sparsity pattern as originally provided.
            pub fn refactor_with_options<'a>(
                &mut self,
                mtx: &SparseMatrix<'a, $t>,
                options: SparseNumericFactorOptions<'a, $t>,
            ) -> SparseStatus {
                unsafe {
                    ffi::$refactor_opt(
                        mtx.mtx,
                        &mut self.fact as *mut <$t as SupportedScalar>::FFISparseFactorizationType,
                        options.into(),
                    );
                }
                self.status()
            }
            /// Solve the system `Ax = b` where `xb` is `b` on the input and `x` on the output.
            pub fn solve_in_place(&self, mut xb: impl AsMut<[$t]>) -> SparseStatus {
                let status = self.status();
                if status != SparseStatus::Ok {
                    return status;
                }
                let xb = xb.as_mut();
                let xb = ffi::$dense_vec {
                    count: i32::try_from(xb.len()).unwrap(),
                    data: xb.as_mut_ptr(),
                };
                unsafe { ffi::$solve_in_place(self.fact, xb) }
                self.status()
            }
            /// Solve the system `Ax = b`.
            pub fn solve(&self, mut b: impl AsMut<[$t]>, mut x: impl AsMut<[$t]>) -> SparseStatus {
                let status = self.status();
                if status != SparseStatus::Ok {
                    return status;
                }
                let b = b.as_mut();
                let b = ffi::$dense_vec {
                    count: i32::try_from(b.len()).unwrap(),
                    data: b.as_mut_ptr(),
                };
                let x = x.as_mut();
                let x = ffi::$dense_vec {
                    count: i32::try_from(x.len()).unwrap(),
                    data: x.as_mut_ptr(),
                };
                unsafe { ffi::$solve(self.fact, b, x) }
                self.status()
            }
        }
    };
}

impl_matrix!(
    f32,
    SparseMatrix_Float,
    SparseConvertFromCoordinate_Float,
    SparseFactor_Float,
    SparseFactorNumeric_Float,
    SparseFactorNumericOpt_Float,
    SparseCleanupSparseMatrix_Float,
);

impl_matrix!(
    f64,
    SparseMatrix_Double,
    SparseConvertFromCoordinate_Double,
    SparseFactor_Double,
    SparseFactorNumeric_Double,
    SparseFactorNumericOpt_Double,
    SparseCleanupSparseMatrix_Double,
);

impl_factorization!(
    f32,
    SparseSolve_Float,
    SparseSolveInPlace_Float,
    DenseVector_Float,
    DenseMatrix_Float,
    SparseCleanupOpaqueNumeric_Float,
    SparseRefactor_Float,
    SparseRefactorOpt_Float,
);

impl_factorization!(
    f64,
    SparseSolve_Double,
    SparseSolveInPlace_Double,
    DenseVector_Double,
    DenseMatrix_Double,
    SparseCleanupOpaqueNumeric_Double,
    SparseRefactor_Double,
    SparseRefactorOpt_Double,
);

#[derive(Debug)]
pub struct SparseSymbolicFactorOptions<'a> {
    pub order_method: SparseOrder,
    pub order: Option<&'a mut [i32]>,
    pub ignore_rows_and_columns: Option<&'a mut [i32]>,
}

impl<'a> Default for SparseSymbolicFactorOptions<'a> {
    fn default() -> Self {
        SparseSymbolicFactorOptions {
            order_method: SparseOrder::Default,
            order: None,
            ignore_rows_and_columns: None,
        }
    }
}

impl<'a> From<SparseSymbolicFactorOptions<'a>> for ffi::SparseSymbolicFactorOptions {
    fn from(opts: SparseSymbolicFactorOptions<'a>) -> Self {
        ffi::SparseSymbolicFactorOptions {
            control: ffi::SparseDefaultControl,
            orderMethod: opts.order_method.into(),
            order: opts
                .order
                .map(|x| x.as_mut_ptr())
                .unwrap_or(std::ptr::null_mut()),
            ignoreRowsAndColumns: opts
                .ignore_rows_and_columns
                .map(|x| x.as_mut_ptr())
                .unwrap_or(std::ptr::null_mut()),
            malloc: None,
            free: None,
            reportError: None,
        }
    }
}

#[derive(Debug)]
pub struct SparseNumericFactorOptions<'a, T> {
    pub scaling_method: SparseScaling,
    pub scaling: Option<&'a mut [T]>,
    pub pivot_tolerance: f64,
    pub zero_tolerance: f64,
}

impl<'a, T> Default for SparseNumericFactorOptions<'a, T> {
    fn default() -> Self {
        SparseNumericFactorOptions {
            scaling_method: SparseScaling::default(),
            scaling: None,
            pivot_tolerance: 0.0,
            zero_tolerance: 0.0,
        }
    }
}

impl<'a, T> From<SparseNumericFactorOptions<'a, T>> for ffi::SparseNumericFactorOptions {
    fn from(opts: SparseNumericFactorOptions<'a, T>) -> Self {
        ffi::SparseNumericFactorOptions {
            control: ffi::SparseDefaultControl,
            scalingMethod: opts.scaling_method.into(),
            scaling: opts
                .scaling
                .map(|x| x.as_mut_ptr())
                .unwrap_or(std::ptr::null_mut()) as *mut std::ffi::c_void,
            pivotTolerance: opts.pivot_tolerance,
            zeroTolerance: opts.zero_tolerance,
        }
    }
}

#[derive(Debug)]
pub struct SparseMatrixStructure<'a> {
    inner: ffi::SparseMatrixStructure,
    phantom: PhantomData<&'a ()>,
}

/// Matrix structure cannot be copied or cloned.
unsafe impl<'a> Send for SparseMatrixStructure<'a> {}

impl From<ffi::SparseMatrixStructure> for SparseMatrixStructure<'_> {
    fn from(s: ffi::SparseMatrixStructure) -> Self {
        SparseMatrixStructure {
            inner: s,
            phantom: PhantomData,
        }
    }
}
impl From<SparseMatrixStructure<'_>> for ffi::SparseMatrixStructure {
    fn from(structure: SparseMatrixStructure) -> Self {
        structure.inner
    }
}

/// A slice of all row indices represented in this BCSC matrix.
///
/// Each index corresponds to an structurally non-zero entry in the matrix.
fn sparse_matrix_structure_indices(structure: &ffi::SparseMatrixStructure) -> &[i32] {
    let num_cols = usize::try_from(structure.columnCount).unwrap();
    unsafe {
        std::slice::from_raw_parts(
            structure.rowIndices,
            usize::try_from(*structure.columnStarts.add(num_cols)).unwrap(),
        )
    }
}
/// A slice of all column offsets represented in this BCSC matrix.
///
/// Each offset (except for one) corresponds to a column in the matrix. The last
/// offset corresponds to the total number of structural non-zeros in the
/// matrix.
fn sparse_matrix_structure_offsets(structure: &ffi::SparseMatrixStructure) -> &[i64] {
    let num_cols = usize::try_from(structure.columnCount).unwrap();
    unsafe { std::slice::from_raw_parts(structure.columnStarts, num_cols + 1) }
}

impl<'a> SparseMatrixStructure<'a> {
    /// Create a sparse BCSC matrix structure from raw indices and offsets.
    pub fn from_raw_parts(
        num_rows: usize,
        num_cols: usize,
        block_size: usize,
        attributes: SparseAttributes,
        indices: &mut [i32],
        offsets: &mut [i64],
    ) -> Self {
        Self {
            inner: ffi::SparseMatrixStructure {
                rowCount: i32::try_from(num_rows).unwrap(),
                columnCount: i32::try_from(num_cols).unwrap(),
                columnStarts: offsets.as_mut_ptr(),
                rowIndices: indices.as_mut_ptr(),
                attributes: attributes.with_owned(false).into(),
                blockSize: u8::try_from(block_size).unwrap(),
            },
            phantom: PhantomData,
        }
    }
    /// A slice of all row indices represented in this BCSC matrix.
    ///
    /// Each index corresponds to an structurally non-zero entry in the matrix.
    pub fn indices(&self) -> &[i32] {
        sparse_matrix_structure_indices(&self.inner)
    }
    /// A slice of all column offsets represented in this BCSC matrix.
    ///
    /// Each offset (except for one) corresponds to a column in the matrix. The last
    /// offset corresponds to the total number of structural non-zeros in the
    /// matrix.
    pub fn offsets(&self) -> &[i64] {
        sparse_matrix_structure_offsets(&self.inner)
    }
    /// Factors the matrix represented by this structure symbolically.
    pub fn symbolic_factor(self, ty: SparseFactorizationType) -> SparseSymbolicFactorization {
        unsafe { ffi::SparseFactorSymbolic(ty.into(), self.into()).into() }
    }
    // Currently opts cause memory issues. TODO: resolve this.
    ///// Factors the matrix represented by this structure symbolically using the given options.
    //pub fn symbolic_factor_with_options(self, ty: SparseFactorizationType, options: SparseSymbolicFactorOptions) -> SparseSymbolicFactorization {
    //    unsafe { ffi::SparseFactorSymbolicOpt(ty.into(), self.into(), options.into()).into() }
    //}
}

#[derive(Debug)]
pub struct SparseSymbolicFactorization {
    factorization: ffi::SparseOpaqueSymbolicFactorization,
}

/// SparseSymbolicFactorization cannot be copied or cloned.
unsafe impl Send for SparseSymbolicFactorization {}

impl From<ffi::SparseOpaqueSymbolicFactorization> for SparseSymbolicFactorization {
    fn from(other: ffi::SparseOpaqueSymbolicFactorization) -> Self {
        Self {
            factorization: other,
        }
    }
}

impl From<SparseSymbolicFactorization> for ffi::SparseOpaqueSymbolicFactorization {
    fn from(other: SparseSymbolicFactorization) -> Self {
        other.factorization
    }
}

impl SparseSymbolicFactorization {
    /// Returns the status of this factorization state.
    pub fn status(&self) -> SparseStatus {
        self.factorization.status.into()
    }
    pub fn new(
        factorization_type: SparseFactorizationType,
        matrix_structure: SparseMatrixStructure,
    ) -> Self {
        matrix_structure.symbolic_factor(factorization_type)
    }
}

/*
 * Iterative method types and implementations.
 */

/// The status of an iterative linear solve.
///
/// This enum is returned by solves implemented in `SparseIterativeMethod`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SparseIterativeStatus {
    /// All solution vectors converged.
    Converged,
    /// One or more solutions failed to converge in the maximum number of iterations.
    MaxIterations,
    /// There was an error with one or more user-supplied parameters.
    ParameterError,
    /// Problem determined to be sufficiently ill conditioned convergence is unlikely.
    IllConditioned,
    /// Some internal failure occurred (e.g. memory allocation failed).
    InternalError,
}

impl From<ffi::SparseIterativeStatus_t> for SparseIterativeStatus {
    fn from(st: ffi::SparseIterativeStatus_t) -> Self {
        match st {
            ffi::SparseIterativeConverged => SparseIterativeStatus::Converged,
            ffi::SparseIterativeMaxIterations => SparseIterativeStatus::MaxIterations,
            ffi::SparseIterativeParameterError => SparseIterativeStatus::ParameterError,
            ffi::SparseIterativeIllConditioned => SparseIterativeStatus::IllConditioned,
            ffi::SparseIterativeInternalError => SparseIterativeStatus::InternalError,
            _ => panic!("Invalid factorization type"),
        }
    }
}

impl From<SparseIterativeStatus> for ffi::SparseIterativeStatus_t {
    fn from(st: SparseIterativeStatus) -> Self {
        match st {
            SparseIterativeStatus::Converged => ffi::SparseIterativeConverged,
            SparseIterativeStatus::MaxIterations => ffi::SparseIterativeMaxIterations,
            SparseIterativeStatus::ParameterError => ffi::SparseIterativeParameterError,
            SparseIterativeStatus::IllConditioned => ffi::SparseIterativeIllConditioned,
            SparseIterativeStatus::InternalError => ffi::SparseIterativeInternalError,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct SparseCGOptions {
    options: ffi::SparseCGOptions,
}

impl From<ffi::SparseCGOptions> for SparseCGOptions {
    fn from(options: ffi::SparseCGOptions) -> Self {
        Self { options }
    }
}

impl From<SparseCGOptions> for ffi::SparseCGOptions {
    fn from(other: SparseCGOptions) -> Self {
        other.options
    }
}

#[derive(Debug, Copy, Clone)]
pub struct SparseGMRESOptions {
    options: ffi::SparseGMRESOptions,
}

impl From<ffi::SparseGMRESOptions> for SparseGMRESOptions {
    fn from(options: ffi::SparseGMRESOptions) -> Self {
        Self { options }
    }
}

impl From<SparseGMRESOptions> for ffi::SparseGMRESOptions {
    fn from(other: SparseGMRESOptions) -> Self {
        other.options
    }
}

#[derive(Debug, Copy, Clone)]
pub enum SparseLSMRConvergenceTest {
    /// Use accelerate's default convergence test:
    ///   || A^Tb - A^TAx ||_2 < rtol * || A^Tb - A^TAx_0 ||_2 + atol.
    Default,
    /// Use the convergence test of Fong and Saunders:
    ///   Either `|| b-Ax ||_2 < btol * || b ||_2 + atol * || A ||_2 || x ||_2`
    ///   or     `|| A^T (b-Ax) ||_2 < atol * || A ||_2 * || A-bx ||_2`
    ///   or     Estimated condition of matrix >= conditionLimit                */
    FongSaunders,
}

impl Default for SparseLSMRConvergenceTest {
    fn default() -> Self {
        SparseLSMRConvergenceTest::Default
    }
}

impl From<ffi::SparseLSMRConvergenceTest_t> for SparseLSMRConvergenceTest {
    fn from(ct: ffi::SparseLSMRConvergenceTest_t) -> Self {
        match ct {
            c if c == ffi::SparseLSMRCTFongSaunders as i32 => {
                SparseLSMRConvergenceTest::FongSaunders
            }
            _ => SparseLSMRConvergenceTest::Default,
        }
    }
}

impl Into<ffi::SparseLSMRConvergenceTest_t> for SparseLSMRConvergenceTest {
    fn into(self) -> ffi::SparseLSMRConvergenceTest_t {
        (match self {
            SparseLSMRConvergenceTest::FongSaunders => ffi::SparseLSMRCTFongSaunders,
            SparseLSMRConvergenceTest::Default => ffi::SparseLSMRCTDefault,
        }) as i32
    }
}

#[derive(Debug, Copy, Clone, Default)]
pub struct SparseLSMROptions {
    /// Damping parameter, if non-zero the actual problem solved is
    ///       min_x || Ax-b ||_2 + lambda || x ||_2.
    /// Using damping can often allow the iteration to converge on ill-conditioned
    /// systems.
    lambda: f64,
    /// Number of vectors used for local orthagonalization.
    /// If nvec<=0, no orthagonalization is performed.
    nvec: u32,
    /// Which convergence test to use. See definition of
    /// SparseLSMRConvergenceTest_t for further information.
    convergence_test: SparseLSMRConvergenceTest,
    /// Either absolute tolerance (default test) or A tolerance (Fong-Saunders
    /// test). In the Fong and Saunders case, it should hold an estimate of the
    /// relative error in the data defining the matrix A. For example, if A is
    /// accurate to about 6 digits, set atol = 1.0e-6. In the Fong and Saunders
    /// case, if atol is 0.0, it is treated as machine epsilon. If using the
    /// default test, a value of 0.0 is treated as 0.0.
    atol: f64,
    /// Relative convergence tolerance (default test only).
    /// If rtol = 0.0, default value of sqrt(epsilon) is used.
    /// If negative, rtol is treated as 0.0 (default is not used).
    rtol: f64,
    /// b tolerance (Fong-Saunders test only). It should hold an estimate of the
    /// relative error in the data defining the rhs b. For example, if b is
    /// accurate to about 6 digits, set btol = 1.0e-6. If btol is zero, it
    /// is treated as machine epsilon.
    btol: f64,
    /// Condition number limit (Fong-Saunders test). Iterations will be terminated
    /// if a computed estimate of cond(Abar) exceeds this value. This is intended
    /// to prevent certain small or zero singular values of A or Abar from coming
    /// into effect and causing unwanted growth in the computed solution.
    /// `condition_limit` and `lambda` may be used separately or together to regularize
    /// ill-conditioned systems.
    /// Normally, conlim should be in the range 1000 to 1/eps.
    /// Suggested value:
    /// conlim = 1/(100*eps)  for compatible systems,
    /// conlim = 1/(10*sqrt(eps)) for least squares.
    /// If `condition_limit` is 0.0, it is treated as 1/eps.
    condition_limit: f64,
    /// Maximum number of iterations to perform. If 0, the default value of 4n
    /// is used.
    /// However, if a good preconditioner is available and/or the matrix is well
    /// conditioned such that singular values are clustered, a value of n/2 may
    /// be more approriate.
    max_iterations: u32,
}

impl From<ffi::SparseLSMROptions> for SparseLSMROptions {
    fn from(options: ffi::SparseLSMROptions) -> Self {
        dbg!(options);
        Self {
            lambda: options.lambda,
            nvec: options.nvec as u32,
            convergence_test: options.convergenceTest.into(),
            atol: options.atol,
            rtol: options.rtol,
            btol: options.btol,
            condition_limit: options.conditionLimit,
            max_iterations: options.maxIterations as u32,
        }
    }
}

impl From<SparseLSMROptions> for ffi::SparseLSMROptions {
    fn from(other: SparseLSMROptions) -> Self {
        dbg!(other);
        ffi::SparseLSMROptions {
            lambda: other.lambda,
            nvec: other.nvec as i32,
            convergenceTest: other.convergence_test.into(),
            atol: other.atol,
            rtol: other.rtol,
            btol: other.btol,
            conditionLimit: other.condition_limit,
            maxIterations: other.max_iterations as i32,
            reportError: None,
            reportStatus: None,
        }
    }
}

#[derive(Copy, Clone)]
pub struct SparseIterativeMethod {
    method: ffi::SparseIterativeMethod,
}

impl From<ffi::SparseIterativeMethod> for SparseIterativeMethod {
    fn from(method: ffi::SparseIterativeMethod) -> Self {
        Self { method }
    }
}

impl From<SparseIterativeMethod> for ffi::SparseIterativeMethod {
    fn from(other: SparseIterativeMethod) -> Self {
        other.method
    }
}

impl SparseIterativeMethod {
    /*
     * Constructors
     */

    /// Constructs a new Conjugate Gradient iterative method with default options.
    pub fn cg() -> Self {
        SparseIterativeMethod {
            // SAFETY: ffi function constructs a Copy type.
            method: unsafe { ffi::SparseConjugateGradientDefault() },
        }
    }
    /// Constructs a new Conjugate Gradient iterative method with the given options.
    pub fn cg_with_options(options: SparseCGOptions) -> Self {
        SparseIterativeMethod {
            // SAFETY: ffi function constructs a Copy type.
            method: unsafe { ffi::SparseConjugateGradientOpt(options.into()) },
        }
    }
    /// Constructs a new GMRES iterative method with default options.
    pub fn gmres() -> Self {
        SparseIterativeMethod {
            // SAFETY: ffi function constructs a Copy type.
            method: unsafe { ffi::SparseGMRESDefault() },
        }
    }
    /// Constructs a new GMRES iterative method with the given options.
    pub fn gmres_with_options(options: SparseGMRESOptions) -> Self {
        SparseIterativeMethod {
            // SAFETY: ffi function constructs a Copy type.
            method: unsafe { ffi::SparseGMRESOpt(options.into()) },
        }
    }
    /// Constructs a new LSMR iterative method with default options.
    pub fn lsmr() -> Self {
        SparseIterativeMethod {
            // SAFETY: ffi function constructs a Copy type.
            method: unsafe { ffi::SparseLSMRDefault() },
        }
    }
    /// Constructs a new LSMR iterative method with the given options.
    pub fn lsmr_with_options(options: SparseLSMROptions) -> Self {
        SparseIterativeMethod {
            // SAFETY: ffi function constructs a Copy type.
            method: unsafe { ffi::SparseLSMROpt(options.into()) },
        }
    }
}

/*
 * Solve functions
 */

#[derive(Debug, Copy, Clone)]
pub enum SparsePreconditioner {
    None,
    User,
    Diagonal,
    DiagScaling,
}

impl Default for SparsePreconditioner {
    fn default() -> Self {
        SparsePreconditioner::None
    }
}

impl From<SparsePreconditioner> for ffi::SparsePreconditioner_t {
    fn from(precond: SparsePreconditioner) -> Self {
        (match precond {
            SparsePreconditioner::None => ffi::SparsePreconditionerNone,
            SparsePreconditioner::User => ffi::SparsePreconditionerUser,
            SparsePreconditioner::Diagonal => ffi::SparsePreconditionerDiagonal,
            SparsePreconditioner::DiagScaling => ffi::SparsePreconditionerDiagScaling,
        }) as i32
    }
}

pub trait IterativeSolve<T, M> {
    fn solve(&self, a: &M, b: impl AsMut<[T]>, x: impl AsMut<[T]>) -> SparseIterativeStatus;
    fn solve_precond(
        &self,
        a: &M,
        b: impl AsMut<[T]>,
        x: impl AsMut<[T]>,
        precond: SparsePreconditioner,
    ) -> SparseIterativeStatus;
    fn solve_op(
        &self,
        op: impl Fn(bool, bool, &[T], &mut [T]),
        b: impl AsMut<[T]>,
        x: impl AsMut<[T]>,
    ) -> SparseIterativeStatus;
}

macro_rules! impl_iterative_solve {
    ($t:ident, $opp:ident, $solve:ident, $solve_op:ident, $solve_precond:ident, $dense_vec:ident, $(,)?) => {
        impl<'a> IterativeSolve<$t, SparseMatrix<'a, $t>> for SparseIterativeMethod {
            fn solve(
                &self,
                a: &SparseMatrix<'a, $t>,
                mut b: impl AsMut<[$t]>,
                mut x: impl AsMut<[$t]>,
            ) -> SparseIterativeStatus {
                let b = b.as_mut();
                let b = ffi::$dense_vec {
                    count: i32::try_from(b.len()).unwrap(),
                    data: b.as_mut_ptr(),
                };
                let x = x.as_mut();
                let x = ffi::$dense_vec {
                    count: i32::try_from(x.len()).unwrap(),
                    data: x.as_mut_ptr(),
                };
                unsafe { ffi::$solve(self.method, a.mtx, b, x) }.into()
            }
            fn solve_precond(
                &self,
                a: &SparseMatrix<'a, $t>,
                mut b: impl AsMut<[$t]>,
                mut x: impl AsMut<[$t]>,
                precond: SparsePreconditioner,
            ) -> SparseIterativeStatus {
                let b = b.as_mut();
                let b = ffi::$dense_vec {
                    count: i32::try_from(b.len()).unwrap(),
                    data: b.as_mut_ptr(),
                };
                let x = x.as_mut();
                let x = ffi::$dense_vec {
                    count: i32::try_from(x.len()).unwrap(),
                    data: x.as_mut_ptr(),
                };
                unsafe { ffi::$solve_precond(self.method, a.mtx, b, x, precond.into()) }.into()
            }
            fn solve_op(
                &self,
                op: impl Fn(bool, bool, &[$t], &mut [$t]),
                mut b: impl AsMut<[$t]>,
                mut x: impl AsMut<[$t]>,
            ) -> SparseIterativeStatus {
                let b = b.as_mut();
                let b = ffi::$dense_vec {
                    count: i32::try_from(b.len()).unwrap(),
                    data: b.as_mut_ptr(),
                };
                let x = x.as_mut();
                let x = ffi::$dense_vec {
                    count: i32::try_from(x.len()).unwrap(),
                    data: x.as_mut_ptr(),
                };

                let ffi_op = move |accumulate: bool,
                                   trans: ffi::CBLAS_TRANSPOSE,
                                   x: ffi::$dense_vec,
                                   y: ffi::$dense_vec| {
                    let x = unsafe { std::slice::from_raw_parts(x.data, x.count.max(0) as usize) };
                    let y =
                        unsafe { std::slice::from_raw_parts_mut(y.data, y.count.max(0) as usize) };
                    op(accumulate, trans == ffi::CBLAS_TRANSPOSE_CblasTrans, x, y);
                };

                let ffi_op_block = ffi::ConcreteBlock::new(ffi_op);
                let status = unsafe { ffi::$solve_op(self.method, &*ffi_op_block, b, x) }.into();
                status
            }
        }
    };
}

impl_iterative_solve!(
    f32,
    ffi_oppf32,
    SparseSolveIterative_Float,
    SparseSolveIterativeOp_Float,
    SparseSolveIterativePrecond_Float,
    DenseVector_Float,
);
impl_iterative_solve!(
    f64,
    ffi_oppf64,
    SparseSolveIterative_Double,
    SparseSolveIterativeOp_Double,
    SparseSolveIterativePrecond_Double,
    DenseVector_Double,
);

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn from_raw() {
        // Define the sparsity pattern for `a`.
        let row_indices = vec![0, 1, 1, 2];
        let col_indices = vec![2, 0, 2, 1];

        let a_values = vec![10.0f32, 20.0, 5.0, 50.0];
        let a = SparseMatrixF32::from_coordinate(
            3,
            3,
            1,
            SparseAttributes::new(),
            &row_indices,
            &col_indices,
            &a_values,
        );

        let mut indices = vec![1, 2, 0, 1];
        let mut offsets = vec![0, 1, 2, 4];
        let mut data = vec![20.0f32, 50.0, 10.0, 5.0];
        assert_eq!(a.data(), data.as_slice());
        assert_eq!(a.indices(), indices.as_slice());
        assert_eq!(a.offsets(), offsets.as_slice());

        let a_raw = SparseMatrixF32::from_raw_parts(
            3,
            3,
            1,
            SparseAttributes::new(),
            &mut indices,
            &mut offsets,
            &mut data,
        );

        assert_eq!(a_raw.data(), a.data());
        assert_eq!(a_raw.indices(), a.indices());
        assert_eq!(a_raw.offsets(), a.offsets());
    }

    #[test]
    fn iterative_solve() {
        // Construct matrix `a` from raw parts.
        let mut row_indices = vec![0, 1, 3, 0, 1, 2, 3, 1, 2];
        let mut values = vec![2.0, -0.2, 2.5, 1.0, 3.2, -0.1, 1.1, 1.4, 0.5];
        let mut column_starts = vec![0, 3, 7, 9];
        let a = SparseMatrixF64::from_raw_parts(
            4,
            3,
            1,
            SparseAttributes::new(),
            &mut row_indices,
            &mut column_starts,
            &mut values,
        );

        // Righ-hand side.
        let b = vec![1.200, 1.013, 0.205, -0.172];

        // Vector of unknowns.
        let mut x = vec![0.0; 3];

        // Expected solution
        let exp_sol = vec![0.1, 0.2, 0.3];

        let lsmr = SparseIterativeMethod::lsmr();
        lsmr.solve(&a, b, x.as_mut_slice());

        for (actual, expected) in x.iter().zip(exp_sol.iter()) {
            assert!((*actual - *expected).abs() < 0.001);
        }
    }

    #[test]
    fn iterative_solve_op() {
        // Righ-hand side.
        let b = vec![1.200, 1.013, 0.205, -0.172];

        // Vector of unknowns.
        let mut x = vec![0.0; 3];

        // Expected solution
        let exp_sol = vec![0.1, 0.2, 0.3];

        // In practice this would be described by some sparse matrix.
        let mtx = [
            [2.0, 1.0, 0.0],
            [-0.2, 3.2, 1.4],
            [0.0, -0.1, 0.5],
            [2.5, 1.1, 0.0],
        ];

        // Define the product operator.
        let op = |accumulate: bool, transpose: bool, x: &[f64], y: &mut [f64]| {
            if !accumulate {
                y.iter_mut().for_each(|a| *a = 0.0);
            }
            if transpose {
                for j in 0..4 {
                    for i in 0..3 {
                        y[i] += mtx[j][i] * x[j];
                    }
                }
            } else {
                for j in 0..4 {
                    for i in 0..3 {
                        y[j] += mtx[j][i] * x[i];
                    }
                }
            }
        };

        let lsmr = SparseIterativeMethod::lsmr();
        dbg!(lsmr.solve_op(op, b, x.as_mut_slice()));

        for (actual, expected) in x.iter().zip(exp_sol.iter()) {
            assert!((*actual - *expected).abs() < 0.001);
        }
    }

    #[test]
    fn factor_and_solve() {
        // Define the sparsity pattern of matrices `A0` and `A1`.
        let row_indices = vec![0, 1, 1, 2];
        let col_indices = vec![2, 0, 2, 1];

        // Create the single-precision coefficient matrix _A0_.
        let a0_values = vec![10.0f32, 20.0, 5.0, 50.0];
        let a0 = SparseMatrixF32::from_coordinate(
            3,
            3,
            1,
            SparseAttributes::new(),
            &row_indices,
            &col_indices,
            &a0_values,
        );
        assert_eq!(a0.data(), &[20.0f32, 50.0, 10.0, 5.0][..]);
        assert_eq!(a0.indices(), &[1, 2, 0, 1][..]);
        assert_eq!(a0.offsets(), &[0, 1, 2, 4][..]);

        // Create the double-precision coefficient matrix _A1_.
        let a1_values = vec![5.0f64, 10.0, 2.5, 25.0];
        let a1 = SparseMatrixF64::from_coordinate(
            3,
            3,
            1,
            SparseAttributes::new(),
            &row_indices,
            &col_indices,
            &a1_values,
        );
        assert_eq!(a1.data(), &[10.0f64, 25.0, 5.0, 2.5][..]);
        assert_eq!(a1.indices(), &[1, 2, 0, 1][..]);
        assert_eq!(a1.offsets(), &[0, 1, 2, 4][..]);

        // Compute the symbolic factorization from the structure of either coefficient matrix.
        let structure = a0.structure();
        let symbolic_factorization = structure.symbolic_factor(SparseFactorizationType::QR);

        // Factorize _A0_ using the symbolic factorization.
        let factorization0 = a0.factor_with(&symbolic_factorization);

        // Solve _A0 · x = b0_ in place.
        let mut b0_values = vec![30.0f32, 35.0, 100.0];
        factorization0.solve_in_place(&mut b0_values);

        // Expected solution
        let exp_sol0 = vec![1.0, 2.0, 3.0];

        for (actual, expected) in b0_values.iter().zip(exp_sol0.iter()) {
            assert_eq!(*actual, *expected);
        }

        // Factorize _A1_ using the symbolic factorization.
        let factorization1 = a1.factor_with(&symbolic_factorization);

        // Solve _A1 · x = b1_ in place.
        let mut b1_values = vec![60.0f64, 70.0, 200.0];
        factorization1.solve_in_place(&mut b1_values);

        // Expected solution
        let exp_sol1 = vec![4.0, 8.0, 12.0];

        for (actual, expected) in b1_values.iter().zip(exp_sol1.iter()) {
            assert_eq!(*actual, *expected);
        }

        // Factorize _A0_ again using the symbolic factorization.
        // This checks that a0 hasn't somehow been compromized with the previous factorization.
        let factorization0 = a0.factor_with(&symbolic_factorization);
        let mut b0_values = vec![30.0f32, 35.0, 100.0];
        factorization0.solve_in_place(&mut b0_values);
        for (actual, expected) in b0_values.iter().zip(exp_sol0.iter()) {
            assert_eq!(*actual, *expected);
        }
    }
    #[test]
    fn solve_singular() {
        // [0 0 .]
        // [. 0 .]
        // [0 . 0]
        let row_indices = vec![0, 1, 1, 2];
        let col_indices = vec![2, 0, 2, 1];
        let a0_values = vec![1.0f32, 0.0, 1.0, 1.0];

        let a0 = SparseMatrixF32::from_coordinate(
            3,
            3,
            1,
            SparseAttributes::new(),
            &row_indices,
            &col_indices,
            &a0_values,
        );
        assert_eq!(a0.data(), &[0.0f32, 1.0, 1.0, 1.0][..]);
        assert_eq!(a0.indices(), &[1, 2, 0, 1][..]);
        assert_eq!(a0.offsets(), &[0, 1, 2, 4][..]);

        // Compute the symbolic factorization from the structure of either coefficient matrix.
        let structure = a0.structure();
        let symbolic_factorization = structure.symbolic_factor(SparseFactorizationType::QR);
        dbg!(symbolic_factorization.status());
        let factorization0 = a0.factor_with(&symbolic_factorization);
        dbg!(factorization0.status());

        let mut b0_values = vec![30.0f32, 35.0, 100.0];
        let status = factorization0.solve_in_place(&mut b0_values);
        dbg!(&status);
        dbg!(&b0_values);
    }
}
