from __future__ import absolute_import

from sage.matrix.matrix2 cimport Matrix
from sage.matrix.matrix_sparse cimport Matrix_sparse
from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse
from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse

from cysignals.memory cimport sig_malloc, sig_free

cdef tuple faugere_lachartre(Matrix_sparse A, Matrix_sparse C, list blocks, list cols)

cdef Matrix_sparse set_up_matrix(list rows, list columns)

cdef Matrix_sparse augment(Matrix_sparse)

cdef bint is_zero(Matrix_sparse A)

cdef tuple split_along_rows(Matrix_sparse A, list rows)     

cdef split_along_columns(Matrix_sparse A, list cols)

cdef tuple reorder(Matrix_sparse AB, Matrix_sparse CD, list blocks, list cols)

cdef list get_relevant_matrix_data(Matrix_sparse A, list blocks)

###################################################
# arithmetic
###################################################

cdef Matrix_sparse _matrix_times_matrix_(Matrix_sparse self, Matrix_sparse right)

cdef void trsm(Matrix_sparse A, Matrix_sparse B)

cdef void diff(Matrix_sparse A, Matrix_sparse B)

###################################################
# rational echelon computations
###################################################

cdef Matrix_sparse echelon(Matrix_sparse, bint transformation = *)

cdef void block_echelon(Matrix_sparse self, list blocks)