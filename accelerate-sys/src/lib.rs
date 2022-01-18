#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iterative_solve() {
        let mut row_indices = vec![0, 1, 3, 0, 1, 2, 3, 1, 2];
        let mut values = vec![2.0, -0.2, 2.5, 1.0, 3.2, -0.1, 1.1, 1.4, 0.5];
        let mut column_starts = vec![0, 3, 7, 9];

        let attributes = SparseAttributes_t {
            _bitfield_align_1: [0; 0],
            _bitfield_1: SparseAttributes_t::new_bitfield_1(
                false,
                SparseTriangle_SparseLowerTriangle as u8, // ignored
                SparseOrdinary,
                0,
                false,
            ),
            __bindgen_padding_0: 0,
        };

        let structure = SparseMatrixStructure {
            rowCount: 4,
            columnCount: 3,
            columnStarts: column_starts.as_mut_ptr(),
            rowIndices: row_indices.as_mut_ptr(),
            attributes,
            blockSize: 1,
        };
        let A = SparseMatrix_Double {
            structure: structure,
            data: values.as_mut_ptr(),
        };

        let mut bValues = vec![1.200, 1.013, 0.205, -0.172];
        let b = DenseVector_Double {
            count: 4,
            data: bValues.as_mut_ptr(),
        };

        let mut xValues = vec![0.00, 0.00, 0.00];
        let x = DenseVector_Double {
            count: 3,
            data: xValues.as_mut_ptr(),
        };

        let lsmr = unsafe { SparseLSMRDefault() };

        unsafe {
            SparseSolveIterative_Double(lsmr, A, b, x);
        }

        dbg!(&xValues);

        for (actual, expected) in xValues.iter().zip([0.10, 0.20, 0.30]) {
            assert!((actual - expected).abs() < 1e-3);
        }
    }

    #[test]
    fn factor_and_solve() {
        let mut row_indices = vec![0, 1, 3, 1, 2, 3, 2, 3];
        let mut values = vec![10.0, 1.0, 2.5, 12.0, -0.3, 1.1, 9.5, 6.0];
        let mut column_starts = vec![0, 3, 6, 7, 8];

        let attributes = SparseAttributes_t {
            _bitfield_align_1: [0; 0],
            _bitfield_1: SparseAttributes_t::new_bitfield_1(
                false,
                SparseTriangle_SparseLowerTriangle as u8, // ignored
                SparseSymmetric,
                0,
                false,
            ),
            __bindgen_padding_0: 0,
        };

        let structure = SparseMatrixStructure {
            rowCount: 4,
            columnCount: 4,
            columnStarts: column_starts.as_mut_ptr(),
            rowIndices: row_indices.as_mut_ptr(),
            attributes,
            blockSize: 1,
        };
        let A = SparseMatrix_Double {
            structure: structure,
            data: values.as_mut_ptr(),
        };

        let mut bValues = vec![2.20, 2.85, 2.79, 2.87];
        let b = DenseVector_Double {
            count: 4,
            data: bValues.as_mut_ptr(),
        };

        let mut xValues = vec![0.00, 0.00, 0.00, 0.00];
        let x = DenseVector_Double {
            count: 4,
            data: xValues.as_mut_ptr(),
        };

        let LLT = unsafe { SparseFactor_Double(SparseFactorizationCholesky as u8, A) };

        unsafe {
            SparseSolve_Double(LLT, b, x);
        }

        dbg!(&x);

        for (actual, expected) in xValues.iter().zip([0.10, 0.20, 0.30, 0.40]) {
            assert!((actual - expected).abs() < 1e-10);
        }
    }

    #[test]
    fn build_matrix() {
        let mut row_indices = vec![0, 1, 3, 0, 1, 2, 3, 1, 2];
        let mut values = vec![2.0, -0.2, 2.5, 1.0, 3.2, -0.1, 1.1, 1.4, 0.5];
        let mut column_starts = vec![0, 3, 7, 9];

        let attributes = SparseAttributes_t {
            _bitfield_align_1: [0; 0],
            _bitfield_1: SparseAttributes_t::new_bitfield_1(
                false,
                SparseTriangle_SparseUpperTriangle as u8, // ignored
                SparseOrdinary,
                0,
                false,
            ),
            __bindgen_padding_0: 0,
        };

        let structure = SparseMatrixStructure {
            rowCount: 4,
            columnCount: 3,
            columnStarts: column_starts.as_mut_ptr(),
            rowIndices: row_indices.as_mut_ptr(),
            attributes,
            blockSize: 1,
        };
        let _A = SparseMatrix_Double {
            structure: structure,
            data: values.as_mut_ptr(),
        };
    }
}
