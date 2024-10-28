from __future__ import absolute_import

from sage.matrix.matrix2 cimport Matrix
from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse

from cysignals.memory cimport sig_malloc, sig_free

from sage.libs.gmp.mpq cimport *
from sage.modules.vector_modn_sparse cimport *


###################################################
# modn auxiliary
###################################################

cdef bint is_zero_modn(Matrix_modn_sparse A)

cdef tuple split_along_rows_modn(Matrix_modn_sparse A, list rows)

cdef tuple split_along_columns_modn(Matrix_modn_sparse A, list columns)
   
cdef tuple reorder_modn(Matrix_modn_sparse AB, Matrix_modn_sparse CD, list blocks, list cols)

cdef list get_relevant_matrix_data_modn(Matrix_modn_sparse A, list blocks)

cdef bint vector_equals_modn(c_vector_modint* v, c_vector_modint* w)

###################################################
# modn arithmetic
###################################################

cdef Matrix_modn_sparse _matrix_times_matrix_modn_(Matrix_modn_sparse self, Matrix_modn_sparse right)

cdef void trsm_modn(Matrix_modn_sparse A, Matrix_modn_sparse B)
 
###################################################
# modn echelon computations
###################################################

cdef void block_echelon_modn(Matrix_modn_sparse self, list blocks)