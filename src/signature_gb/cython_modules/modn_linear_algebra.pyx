# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cysignals.signals cimport sig_on, sig_off

from sage.modules.vector_modn_sparse cimport *
from sage.rings.finite_rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract

from .nc_monomial cimport NCMonomial
from .sig_polynomial cimport SigPoly

from signature_gb.cython_modules.linear_algebra cimport augment

from time import time
from bisect import bisect_right

from sage.rings.fast_arith cimport arith_int
cdef arith_int ai
ai = arith_int()
        
############################################################################
############################################################################
# Modn linear algebra
############################################################################
############################################################################


############################################################################    
cdef bint is_zero_modn(Matrix_modn_sparse A):
    cdef Py_ssize_t i
    cdef c_vector_modint v
    for i from 0 <= i < A._nrows:
        v = A.rows[i]
        if v.num_nonzero: return False
    return True

############################################################################    
cdef tuple split_along_rows_modn(Matrix_modn_sparse A, list rows):
    cdef Matrix_modn_sparse B,C
    cdef Py_ssize_t i_a,i_b,i_c,k
            
    B = A.new_matrix(len(rows), A._ncols)
    C = A.new_matrix(A._nrows - len(rows), A._ncols)
            
    for i_b from 0 <= i_b < len(rows):
        B.rows[i_b] = A.rows[rows[i_b]]
    
    i_c = 0
    k = 0
    rows.sort()
    for i_a from 0 <= i_a < A._nrows:
        if k < len(rows) and i_a == rows[k]:
            k += 1
            continue
        C.rows[i_c] = A.rows[i_a]
        i_c += 1
    
    A._nrows = 0
 
    return B,C
############################################################################        
cdef tuple split_along_columns_modn(Matrix_modn_sparse A, list cols):
    cdef Matrix_modn_sparse B,C
    cdef Py_ssize_t i,j,k,l
    cdef c_vector_modint* a
    cdef c_vector_modint* b
      
    B = A.new_matrix(A._nrows, len(cols))
    C = A.new_matrix(A._nrows, A._ncols - len(cols))
     
    cdef set set_cols = set(cols)
    cdef list columns = [i in set_cols for i in range(A._ncols)]
    cdef list new_columns = [None for i in range(A._ncols)]
    for i,j in enumerate(cols): new_columns[j] = i
    new_columns = [i-bisect_right(cols,i) if new_columns[i] == None else new_columns[i] for i in range(A._ncols)]
                
    for i from 0 <= i < A._nrows:
        a = &(A.rows[i])
        b = &(B.rows[i])
    
        j = 0
        while j < a.num_nonzero:
            k = a.positions[j]
            if columns[k]:
                set_entry(b, new_columns[k], a.entries[j])
                a.num_nonzero -= 1
                for l from j <= l < a.num_nonzero:
                    a.entries[l] =  a.entries[l+1]
                    a.positions[l] = a.positions[l+1]   
            else: 
                j += 1
        for l from 0 <= l < a.num_nonzero: a.positions[l] = new_columns[a.positions[l]]
        C.rows[i] = A.rows[i]
                    
    A._nrows = 0
            
    return B,C
############################################################################
cdef tuple reorder_modn(Matrix_modn_sparse AB, Matrix_modn_sparse CD, list blocks, list cols):
    cdef Matrix_modn_sparse A,B,C,D
    cdef Py_ssize_t nr, nc, i, j
    cdef list idxs, rest_cols, pivot_rows, pivot_cols
    cdef c_vector_modint* r 
    
    nr = AB._nrows
    nc = AB._ncols
    
    pivots = set()
    idxs = []
    
    for i from 0 <= i < nr:
        r = &(AB.rows[i])
        # zero row or pivot we already have
        if r.num_nonzero == 0 or r.positions[0] in pivots: continue
        # new pivot
        j = r.positions[0]
        pivots.add(j)
        idxs.append((j,i))
        
    idxs.sort()
    pivot_cols = [t[0] for t in idxs]
    pivot_rows = [t[1] for t in idxs]
        
    AB,_ = split_along_rows_modn(AB,pivot_rows)
        
    A,B = split_along_columns_modn(AB, pivot_cols) 
    C,D = split_along_columns_modn(CD, pivot_cols)
            
    blocks = [i - nr for i in blocks if i > nr]
    cols = [cols[i] for i in range(nc) if not i in pivots]
    
    return A,B,C,D,blocks,cols

############################################################################       
cdef list get_relevant_matrix_data_modn(Matrix_modn_sparse A, list blocks):        
    cdef Py_ssize_t i, j, k
    cdef list data
    cdef bint flag
    cdef c_vector_modint* v
    
    data = []
    R = A._base_ring
        
    j = 0
    for i from 0 <= i < A._nrows:
        if i == blocks[j]: j += 1
        flag = True
        v = &(A.rows[i])
        for k from i < k < blocks[j]:
            if vector_equals_modn(v,&(A.rows[k])):
                flag = False
                break  
        if flag:
            if v.num_nonzero:
                pos,coeffs = zip(*to_list(v))
                coeffs = [R(c) for c in coeffs]
            else: 
                pos = coeffs = None
            data.append((i,pos,coeffs))
    return data
     
############################################################################    
cdef bint vector_equals_modn(c_vector_modint* v, c_vector_modint* w):
    cdef Py_ssize_t i
    
    if v.num_nonzero != w.num_nonzero: return False
    for i from 0 <= i < v.num_nonzero:
        if v.positions[i] != w.positions[i]: return False
        if v.entries[i] != w.entries[i]: return False
    
    return True  

###################################################
# modn arithmetic
###################################################

cdef Matrix_modn_sparse _matrix_times_matrix_modn_(Matrix_modn_sparse self, Matrix_modn_sparse right):
    cdef Matrix_modn_sparse ans

    cdef c_vector_modint* v
    cdef c_vector_modint* w
    
    cdef int_fast64_t mod = self.p

    # Build a table that gives the nonzero positions in each column of right
    cdef list nonzero_positions_in_columns = [set() for _ in range(right._ncols)]
    cdef Py_ssize_t i, j, k
    
    for i from 0 <= i < right._nrows:
        v = &(right.rows[i])
        for j from 0 <= j < v.num_nonzero:
            (<set> nonzero_positions_in_columns[v.positions[j]]).add(i)

    ans = self.new_matrix(self._nrows, right._ncols)

    # Now do the multiplication, getting each row complete before filling it in.
    cdef set relevant_columns
    cdef int_fast64_t x, y, s
    cdef Py_ssize_t pos
    
    for i from 0 <= i < self._nrows:
        v = &(self.rows[i])
        if not v.num_nonzero: continue
                
        # find all relevant columns
        relevant_columns = set()
        for k from 0 <= k < v.num_nonzero:
            w = &(right.rows[v.positions[k]])
            (<set> relevant_columns).update({w.positions[p] for p in range(w.num_nonzero)})
                
        # compute dot product with relevant columns
        for j in list(relevant_columns):    
            s = 0
            for k from 0 <= k < v.num_nonzero:
                if v.positions[k] in nonzero_positions_in_columns[j]:
                    y = get_entry(&right.rows[v.positions[k]], j)
                    x = v.entries[k] * y
                    s += x
                    s = s % mod
                    set_entry(&ans.rows[i], j, s)

    return ans
###########################################################################

cdef void trsm_modn(Matrix_modn_sparse A, Matrix_modn_sparse B):
    cdef Py_ssize_t i, k, nc
    cdef c_vector_modint* a
    cdef c_vector_modint tmp
    cdef int_fast64_t a_minus, mod
    
    mod = A.p
    nr = B._nrows

    for i from nr - 1 > i >= 0:
        a = &(A.rows[i])
        if a.num_nonzero == 1: continue
        for k from 1 <= k < a.num_nonzero:
            a_minus = mod - a.entries[k]
            add_c_vector_modint_init(&tmp, &B.rows[i], &B.rows[a.positions[k]], a_minus)
            clear_c_vector_modint(&B.rows[i])
            B.rows[i] = tmp

###################################################
# modn echelon computations
###################################################

cdef void block_echelon_modn(Matrix_modn_sparse self, list blocks):
        """
        Replace self by its reduction to reduced row echelon form.
        """
        cdef Py_ssize_t i, r, c, nr, nc, pivot, j, next_block 
        cdef int_fast64_t inv, b, b_minus, mod
        cdef c_vector_modint tmp
        cdef c_vector_modint* row
        cdef c_vector_modint* row_i
        
        mod = self.p
        nr = self._nrows
        nc = self._ncols
              
        for c from 0 <= c < nc:          
            # find pivot row
            pivot = -1
            for r from 0 <= r < nr:
                # find row with pivot
                row = &self.rows[r]
                if row.num_nonzero > 0 and row.positions[0] == c:
                    pivot = r
                    break
            
            if pivot == -1: continue
            sig_on()
            # normalize row    
            row = &self.rows[pivot]
            if row.entries[0] != 1:
                inv = ai.c_inverse_mod_int(row.entries[0], mod)
                scale_c_vector_modint(row, inv)
            
            # find index where next block starts
            next_block = blocks[bisect_right(blocks,pivot)]
                        
            # reduce rows below
            for i from next_block <= i < nr:
                row_i = &self.rows[i]
                b = get_entry(row_i, c)
                if b != 0:
                    add_c_vector_modint_init(&tmp, row_i, row, mod - b)
                    clear_c_vector_modint(row_i)
                    self.rows[i] = tmp
            sig_off()
               
        # normalize rows
        for i from 0 <= i < self._nrows:
            row = &self.rows[i]
            if not row.num_nonzero: continue
            if row.entries[0] != 1:
                inv = ai.c_inverse_mod_int(row.entries[0], mod)
                scale_c_vector_modint(row, inv)