use std::convert::TryFrom;
use std::mem::ManuallyDrop;

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
            SparseTriangle::Upper => ffi::SparseTriangle_SparseUpperTriangle as ffi::SparseTriangle_t,
            SparseTriangle::Lower => ffi::SparseTriangle_SparseLowerTriangle as ffi::SparseTriangle_t,
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
        Self {
            attrs,
        }
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
                    false
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
            SparseFactorizationType::Cholesky => ffi::SparseFactorizationCholesky as ffi::SparseFactorization_t,
            SparseFactorizationType::CholeskyAtA => ffi::SparseFactorizationCholeskyAtA as ffi::SparseFactorization_t,
            SparseFactorizationType::LDLT => ffi::SparseFactorizationLDLT as ffi::SparseFactorization_t,
            SparseFactorizationType::LDLTUnpivoted => ffi::SparseFactorizationLDLTUnpivoted as ffi::SparseFactorization_t,
            SparseFactorizationType::LDLTSBK => ffi::SparseFactorizationLDLTSBK as ffi::SparseFactorization_t,
            SparseFactorizationType::LDLTTPP => ffi::SparseFactorizationLDLTTPP as ffi::SparseFactorization_t,
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
impl From<ffi::SparseOrder_t> for SparseOrder {
    fn from(st: ffi::SparseOrder_t) -> Self {
        match st as u32 {
            ffi::SparseOrderDefault => SparseOrder::Default,
            ffi::SparseOrderUser => SparseOrder::User,
            ffi::SparseOrderAMD => SparseOrder::AMD,
            ffi::SparseOrderMetis => SparseOrder::Metis,
            ffi::SparseOrderCOLAMD => SparseOrder::COLAMD,
            _ => panic!("Invalid factorization type"),
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
            ffi::SparseStatusOK                 => SparseStatus::Ok                 ,
            ffi::SparseFactorizationFailed=> SparseStatus::FactorizationFailed,
            ffi::SparseMatrixIsSingular   => SparseStatus::MatrixIsSingular   ,
            ffi::SparseInternalError      => SparseStatus::InternalError      ,
            ffi::SparseParameterError     => SparseStatus::ParameterError     ,
            ffi::SparseStatusReleased           => SparseStatus::Released           ,
            _ => panic!("Invalid factorization type"),
        }
    }
}

impl From<SparseStatus> for ffi::SparseStatus_t {
    fn from(st: SparseStatus) -> Self {
        match st {
            SparseStatus::Ok                  => ffi::SparseStatusOK                 ,
            SparseStatus::FactorizationFailed => ffi::SparseFactorizationFailed,
            SparseStatus::MatrixIsSingular    => ffi::SparseMatrixIsSingular   ,
            SparseStatus::InternalError       => ffi::SparseInternalError      ,
            SparseStatus::ParameterError      => ffi::SparseParameterError     ,
            SparseStatus::Released            => ffi::SparseStatusReleased           ,
        }
    }
}

macro_rules! impl_matrix {
    ($t:ident, $mtx_rs:ident, $mtx_ffi:ident, $convert:ident, $factor:ident,
     $factorization:ident, $factorization_ffi:ident, $factor_numeric:ident,
     $solve:ident, $solve_in_place:ident, $dense_vec:ident, $dense_mtx:ident,
     $cleanup_mtx:ident, $cleanup_fact:ident) => {
        #[derive(Debug)]
        pub struct $mtx_rs {
            mtx: ManuallyDrop<ffi::$mtx_ffi>,
        }

        impl Drop for $mtx_rs {
            fn drop(&mut self) {
                unsafe {
                    ffi::$cleanup_mtx(ManuallyDrop::take(&mut self.mtx));
                    ManuallyDrop::drop(&mut self.mtx);
                }
            }
        }

        impl $mtx_rs {
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

                $mtx_rs {
                    mtx: ManuallyDrop::new(unsafe { ffi::$convert(
                        num_rows,
                        num_cols,
                        block_count,
                        block_size,
                        attributes.into(),
                        rows.as_ptr(),
                        cols.as_ptr(),
                        values.as_ptr(),
                    ) }),
                }
            }
            pub fn factor(self, factorization_type: SparseFactorizationType) -> $factorization {
                unsafe { ffi::$factor(factorization_type.into(), ManuallyDrop::into_inner(self.mtx)) }.into()
            }
            pub fn factor_with(self, factorization: &SparseSymbolicFactorization) -> $factorization {
                unsafe { ffi::$factor_numeric(factorization.factorization.clone(), ManuallyDrop::into_inner(self.mtx)) }.into()
            }
            pub fn structure(&self) -> SparseMatrixStructure {
                self.mtx.structure.into()
            }
        }
        #[derive(Clone, Debug)]
        pub struct $factorization {
            fact: ManuallyDrop<ffi::$factorization_ffi>,
        }
        impl From<ffi::$factorization_ffi> for $factorization {
            fn from(f: ffi::$factorization_ffi) -> $factorization {
                $factorization { fact: ManuallyDrop::new(f) }
            }
        }
        impl From<$factorization> for ffi::$factorization_ffi {
            fn from(f: $factorization) -> ffi::$factorization_ffi {
                ManuallyDrop::into_inner(f.fact)
            }
        }
        impl Drop for $factorization {
            fn drop(&mut self) {
                unsafe {
                    ffi::$cleanup_fact(ManuallyDrop::take(&mut self.fact));
                    ManuallyDrop::drop(&mut self.fact);
                }
            }
        }

        impl $factorization {
            //pub(crate) fn empty() -> Self {
            //    Self {
            //        fact: ManuallyDrop::new(ffi::$factorization_ffi {
            //            status: SparseStatus::Ok.into(),
            //            attributes: SparseAttributes::new().into(),
            //            symbolicFactorization: SparseSymbolicFactorization::empty().into(),
            //            userFactorStorage: false,
            //            numericFactorization: std::ptr::null_mut(),
            //            solveWorkspaceRequiredStatic: 0,
            //            solveWorkspaceRequiredPerRHS: 0,
            //        })
            //    }
            //}
            pub fn solve_in_place(self, mut xb: impl AsMut<[$t]>) {
                let xb = xb.as_mut();
                let xb = ffi::$dense_vec {
                    count: i32::try_from(xb.len()).unwrap(),
                    data: xb.as_mut_ptr(),
                };
                unsafe { ffi::$solve_in_place(ManuallyDrop::into_inner(self.fact), xb) }
            }
            pub fn solve(self, mut b: impl AsMut<[$t]>, mut x: impl AsMut<[$t]>) {
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
                unsafe { ffi::$solve(ManuallyDrop::into_inner(self.fact), b, x) }
            }
        }
    };
}

impl_matrix!(
    f32,
    SparseMatrixF32,
    SparseMatrix_Float,
    SparseConvertFromCoordinate_Float,
    SparseFactor_Float,
    SparseFactorizationF32,
    SparseOpaqueFactorization_Float,
    SparseFactorNumeric_Float,
    SparseSolve_Float,
    SparseSolveInPlace_Float,
    DenseVector_Float,
    DenseMatrix_Float,
    SparseCleanupSparseMatrix_Float,
    SparseCleanupOpaqueNumeric_Float
);

impl_matrix!(
    f64,
    SparseMatrixF64,
    SparseMatrix_Double,
    SparseConvertFromCoordinate_Double,
    SparseFactor_Double,
    SparseFactorizationF64,
    SparseOpaqueFactorization_Double,
    SparseFactorNumeric_Double,
    SparseSolve_Double,
    SparseSolveInPlace_Double,
    DenseVector_Double,
    DenseMatrix_Double,
    SparseCleanupSparseMatrix_Double,
    SparseCleanupOpaqueNumeric_Double
);

#[derive(Debug)]
pub struct SparseMatrixStructure {
    inner: ffi::SparseMatrixStructure,
}

impl From<ffi::SparseMatrixStructure> for SparseMatrixStructure {
    fn from(s: ffi::SparseMatrixStructure) -> Self {
        SparseMatrixStructure { inner: s }
    }
}
impl From<SparseMatrixStructure> for ffi::SparseMatrixStructure {
    fn from(structure: SparseMatrixStructure) -> Self {
        structure.inner
    }
}

impl SparseMatrixStructure {
    ///// Constructs an empty matrix structure.
    //pub(crate) fn empty() -> Self {
    //    Self {
    //        inner: ffi::SparseMatrixStructure {
    //            rowCount: 0,
    //            columnCount: 0,
    //            columnStarts: std::ptr::null_mut(),
    //            rowIndices: std::ptr::null_mut(),
    //            attributes: SparseAttributes::new().into(),
    //            blockSize: 0
    //        }
    //    }
    //}
    pub fn symbolic_factor(self, ty: SparseFactorizationType) -> SparseSymbolicFactorization {
        unsafe { ffi::SparseFactorSymbolic(ty.into(), self.into()).into() }
    }
}

#[derive(Debug)]
pub struct SparseSymbolicFactorization {
    factorization: ffi::SparseOpaqueSymbolicFactorization,
}

impl From<ffi::SparseOpaqueSymbolicFactorization> for SparseSymbolicFactorization {
    fn from(other: ffi::SparseOpaqueSymbolicFactorization) -> Self {
        Self {
            factorization: other
        }
    }
}

impl From<SparseSymbolicFactorization> for ffi::SparseOpaqueSymbolicFactorization {
    fn from(other: SparseSymbolicFactorization) -> Self {
        other.factorization
    }
}

impl SparseSymbolicFactorization {
    //pub(crate) fn empty() -> Self {
    //    Self {
    //        factorization: ffi::SparseOpaqueSymbolicFactorization {
    //            status: SparseStatus::Ok.into(),
    //            rowCount: 0,
    //            columnCount: 0,
    //            attributes: SparseAttributes::new().into(),
    //            blockSize: 0,
    //            type_: SparseFactorizationType::Cholesky.into(),
    //            factorization: std::ptr::null_mut(),
    //            workspaceSize_Float: 0,
    //            workspaceSize_Double: 0,
    //            factorSize_Float:0,
    //            factorSize_Double:0,
    //        }
    //    }
    //}
    pub fn new(
        factorization_type: SparseFactorizationType,
        matrix_structure: SparseMatrixStructure,
    ) -> Self {
        matrix_structure.symbolic_factor(factorization_type)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

        // Compute the symbolic factorization from the structure of either coefficient matrix.
        let structure = a0.structure();
        let symbolic_factorization = structure.symbolic_factor(SparseFactorizationType::QR);

        // Factorize _A0_ using the symbolic factorization.
        let factorization0 = a0.factor_with(&symbolic_factorization);

        // Solve _A0 · x = b0_ in place.
        let mut b0_values = vec![30.0f32, 35.0, 100.0];
        factorization0.solve_in_place( &mut b0_values);

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
    }
}
