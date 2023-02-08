from __future__ import absolute_import

from sage.matrix.matrix2 cimport Matrix
from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse
from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse

from cysignals.memory cimport sig_malloc, sig_free

from sage.libs.gmp.mpq cimport *
from sage.modules.vector_rational_sparse cimport *


###################################################
# rational auxiliary
###################################################

cdef Matrix_rational_sparse set_up_matrix(list rows, list columns)

cdef bint is_zero(Matrix_rational_sparse A)

cdef tuple split_along_columns(Matrix_rational_sparse A, list columns)

cdef tuple split_along_rows(Matrix_rational_sparse A, list rows)

cdef Matrix_rational_sparse augment(Matrix_rational_sparse)

cdef tuple reorder(Matrix_rational_sparse AB, Matrix_rational_sparse CD, list blocks, list cols)

cdef bint vector_equals(mpq_vector* v, mpq_vector* w)

###################################################
# rational arithmetic
###################################################

cdef Matrix_rational_sparse _matrix_times_matrix_(Matrix_rational_sparse self, Matrix_rational_sparse right)

cdef void trsm(Matrix_rational_sparse A, Matrix_rational_sparse B)

cdef void diff(Matrix_rational_sparse A, Matrix_rational_sparse B)

###################################################
# rational echelon computations
###################################################

cdef Matrix_rational_sparse rational_echelon(Matrix_rational_sparse, bint transformation = *)

cdef Matrix_rational_sparse block_rational_echelon(Matrix_rational_sparse self, list blocks)

###################################################
# multimodular echelon computations
###################################################

cpdef multimodular_echelon_form(Matrix)

cpdef multimodular_matrix_rational_echelon_form_(Matrix)

cpdef modn_echelon_in_place(Matrix_modn_sparse)