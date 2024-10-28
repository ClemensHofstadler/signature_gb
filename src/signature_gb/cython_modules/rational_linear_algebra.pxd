from __future__ import absolute_import

from sage.matrix.matrix2 cimport Matrix
from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse

from cysignals.memory cimport sig_malloc, sig_free

from sage.libs.gmp.mpq cimport *
from sage.modules.vector_rational_sparse cimport *


###################################################
# rational auxiliary
###################################################

cdef bint is_zero_rational(Matrix_rational_sparse A)

cdef tuple split_along_columns_rational(Matrix_rational_sparse A, list columns)

cdef tuple split_along_rows_rational(Matrix_rational_sparse A, list rows)
 
cdef tuple reorder_rational(Matrix_rational_sparse AB, Matrix_rational_sparse CD,  list blocks, list cols)

cdef list get_relevant_matrix_data_rational(Matrix_rational_sparse A, list blocks)

cdef bint vector_equals(mpq_vector* v, mpq_vector* w)

###################################################
# rational arithmetic
###################################################

cdef Matrix_rational_sparse _matrix_times_matrix_rational_(Matrix_rational_sparse self, Matrix_rational_sparse right)

cdef void trsm_rational(Matrix_rational_sparse A, Matrix_rational_sparse B)
 
###################################################
# rational echelon computations
###################################################

cdef Matrix_rational_sparse rational_echelon(Matrix_rational_sparse, bint transformation = *)

cdef void block_echelon_rational(Matrix_rational_sparse self,  list blocks)

cdef bint mpq_vector_equals(mpq_vector* v, mpq_vector* w)