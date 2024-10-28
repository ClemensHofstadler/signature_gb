# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from sage.all import matrix, primes

from .nc_monomial cimport NCMonomial
from .sig_polynomial cimport SigPoly

from .rational_linear_algebra cimport *
from .modn_linear_algebra cimport *
        
############################################################################
############################################################################
# Linear algebra
############################################################################
############################################################################
cdef tuple faugere_lachartre(Matrix_sparse A, Matrix_sparse C, list blocks, list cols):
    cdef Matrix_sparse B,D
        
    # reorder   
    A,B,C,D,blocks,cols = reorder(A,C,blocks,cols)   
                    
    if not C.is_zero(): 
        # compute A^{-1}B  
        trsm(A,B) 
              
        #compute difference D - CA^{-1}B          
        C = _matrix_times_matrix_(C,B)
        diff(D,C)
                                                    
    return D,blocks,cols

############################################################################
cdef Matrix_sparse set_up_matrix(list rows, list columns):
    cdef Matrix_sparse A
    cdef Py_ssize_t nr, nc, i, j
    cdef dict cols
    cdef list lt_cols
    cdef NCMonomial m
    cdef SigPoly f
    
    K = rows[0]._poly._coeffs[0].parent()
   
    nr = len(rows)
    nc = len(columns)
    A = matrix(K,nr,nc,sparse=True)
    
    cols = {m:i for i,m in enumerate(columns)}
                
    for i,f in enumerate(rows):
        for c,m in zip(f._poly._coeffs, f._poly._mons):
            j = cols[m]
            A.set_unsafe(i,j,c)
    
    return A
############################################################################
cdef Matrix_sparse augment(Matrix_sparse A):
    cdef Matrix_sparse B
    cdef Py_ssize_t nr, nc, i, j
    cdef list A_nonzero_pos = <list> A.nonzero_positions(copy=False, column_order=False)
        
    R = A._base_ring
    one = R.one()
    
    nr = A._nrows
    nc = A._ncols
    B = A.new_matrix(nr, nr + nc)
        
    for (i,j) in A_nonzero_pos:
        v = R(A.get_unsafe(i,j))
        B.set_unsafe(i,j,v)
                
    for i from 0 <= i < nr:
        B.set_unsafe(i, i + nc, one) 
          
    return B
############################################################################    
cdef bint is_zero(Matrix_sparse A):
    if isinstance(A,Matrix_rational_sparse): return is_zero_rational(A)
    else: return is_zero_modn(A)

############################################################################    
cdef tuple split_along_rows(Matrix_sparse A, list rows):
    if isinstance(A,Matrix_rational_sparse): return split_along_rows_rational(A, rows)
    else: return split_along_rows_modn(A, rows)
    
############################################################################        
cdef split_along_columns(Matrix_sparse A, list cols):
    if isinstance(A,Matrix_rational_sparse): return split_along_columns_rational(A, cols)
    else: return split_along_columns_modn(A, cols)

############################################################################
cdef tuple reorder(Matrix_sparse AB, Matrix_sparse CD, list blocks, list cols):
    if isinstance(AB,Matrix_rational_sparse): return reorder_rational(AB, CD, blocks, cols)
    else: return reorder_modn(AB, CD, blocks, cols)
       
############################################################################

cdef Matrix_sparse _matrix_times_matrix_(Matrix_sparse self, Matrix_sparse right):
    if isinstance(self,Matrix_rational_sparse): return _matrix_times_matrix_rational_(self, right)
    else: return _matrix_times_matrix_modn_(self, right)
############################################################################

cdef void trsm(Matrix_sparse A, Matrix_sparse B):
    if isinstance(A,Matrix_rational_sparse): trsm_rational(A,B)
    else: trsm_modn(A,B)    
############################################################################ 
       
cdef list get_relevant_matrix_data(Matrix_sparse A, list blocks):
    if isinstance(A,Matrix_rational_sparse): return get_relevant_matrix_data_rational(A,blocks)
    else: return get_relevant_matrix_data_modn(A,blocks)
############################################################################

cdef void diff(Matrix_sparse A, Matrix_sparse B):
    cdef Py_ssize_t i,j
    
    for (i,j) in B.nonzero_positions(copy=False,column_order=False):
        aij = A.get_unsafe(i,j)
        bij = B.get_unsafe(i,j)
        v = aij - bij
        A.set_unsafe(i,j,v)

###################################################
# echelon computations
###################################################

cdef Matrix_sparse echelon(Matrix_sparse self, bint transformation = False):
        """
        Replace self by its reduction to reduced row echelon form.
        """
        if isinstance(self, Matrix_rational_sparse): rational_echelon(self, transformation)
        else: raise NotImplementedError
############################################################################
cdef void block_echelon(Matrix_sparse self, list blocks):
        """
        Replace self by its reduction to reduced row echelon form.
        """
        if isinstance(self, Matrix_rational_sparse): block_echelon_rational(self, blocks)
        else: block_echelon_modn(self, blocks)