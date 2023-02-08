# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cysignals.signals cimport sig_on, sig_off

from sage.matrix.misc import cmp_pivots
from sage.arith.all import previous_prime, CRT_basis

from sage.all import matrix, QQ
from sage.rings.rational cimport Rational

from sage.modules.vector_modn_sparse cimport *
from sage.modules.vector_rational_sparse cimport *

from cython_modules.ncmonomial cimport NCMonomial
from cython_modules.sigpoly cimport SigPoly

from python_modules import global_data

from sage.rings.fast_arith cimport arith_int
cdef arith_int ai
ai = arith_int()

from time import time
from bisect import bisect_right

        
############################################################################
############################################################################
# Linear algebra
############################################################################
############################################################################

###################################################
# rational auxiliaray
###################################################

cdef Matrix_rational_sparse set_up_matrix(list rows, list columns):
    cdef Matrix_rational_sparse A
    cdef Py_ssize_t nr, nc, i, j
    cdef Rational c
    cdef dict cols
    cdef list lt_cols
    cdef NCMonomial m
    cdef SigPoly f
   
    nr = len(rows)
    nc = len(columns)
    A = matrix(QQ,nr,nc,sparse=True)
    
    cols = {m:i for i,m in enumerate(columns)}
                
    for i,f in enumerate(rows):
        for c,m in zip(f._poly._coeffs, f._poly._mons):
            j = cols[m]
            A.set_unsafe(i,j,c)
    
    return A

############################################################################    
cdef bint is_zero(Matrix_rational_sparse A):
    cdef Py_ssize_t i
    cdef mpq_vector* v
    for i from 0 <= i < A._nrows:
        v = &(A._matrix[i])
        if v.num_nonzero: return False
    return True

############################################################################    
cdef tuple split_along_rows(Matrix_rational_sparse A, list rows):
    cdef Matrix_rational_sparse B,C
    cdef Py_ssize_t i_a,i_b,i_c,k
    
    B = A.new_matrix(len(rows), A._ncols)
    C = A.new_matrix(A._nrows - len(rows), A._ncols)
            
    for i_b from 0 <= i_b < len(rows):
        B._matrix[i_b] = A._matrix[rows[i_b]]
    
    i_c = 0
    k = 0
    rows.sort()
    for i_a from 0 <= i_a < A._nrows:
        if k < len(rows) and i_a == rows[k]:
            k += 1
            continue
        C._matrix[i_c] = A._matrix[i_a]
        i_c += 1
    
    A._nrows = 0
 
    return B,C
############################################################################        
cdef tuple split_along_columns(Matrix_rational_sparse A, list cols):
    cdef Matrix_rational_sparse B,C
    cdef Py_ssize_t i,j,k,l
    cdef mpq_vector* a
    cdef mpq_vector* b
    
    B = A.new_matrix(A._nrows, len(cols))
    C = A.new_matrix(A._nrows, A._ncols - len(cols))
     
    cdef set set_cols = set(cols)
    cdef list columns = [i in set_cols for i in range(A._ncols)]
    cdef list new_columns = [None for i in range(A._ncols)]
    for i,j in enumerate(cols): new_columns[j] = i
    new_columns = [i-bisect_right(cols,i) if new_columns[i] == None else new_columns[i] for i in range(A._ncols)]
                
    for i from 0 <= i < A._nrows:
        a = &A._matrix[i]
        b = &B._matrix[i]
    
        j = 0
        while j < a.num_nonzero:
            k = a.positions[j]
            if columns[k]:
                mpq_vector_set_entry(b, new_columns[k], a.entries[j])
                a.num_nonzero -= 1
                for l from j <= l < a.num_nonzero:
                    mpq_set(a.entries[l], a.entries[l+1])
                    a.positions[l] = a.positions[l+1]   
            else: 
                j += 1
        for l from 0 <= l < a.num_nonzero: a.positions[l] = new_columns[a.positions[l]]
        C._matrix[i] = A._matrix[i]
                    
    A._nrows = 0
            
    return B,C
############################################################################
cdef Matrix_rational_sparse augment(Matrix_rational_sparse A):
    cdef Matrix_rational_sparse B
    cdef Py_ssize_t nr, nc, c, r, k
    cdef mpq_vector v,w
    cdef mpq_t one, a
    mpq_init(one)
    mpq_init(a)
    mpq_set_ui(one, 1, 1)
    
    nr = A._nrows
    nc = A._ncols
    B = A.new_matrix(nr, nr + nc)
                
    sig_on()
    for r from 0 <= r < nr:
       mpq_vector_set_entry(&B._matrix[r], nc + r, one)     
       v = A._matrix[r]
       w = B._matrix[r]

       for k from 0 <= k < v.num_nonzero: 
           c = v.positions[k]
           mpq_vector_get_entry(a, &v, c)
           mpq_vector_set_entry(&B._matrix[r], c, a)
    sig_off()
    
    mpq_clear(one)
    mpq_clear(a)
    
    return B
############################################################################
cdef tuple reorder(Matrix_rational_sparse AB, Matrix_rational_sparse CD, list blocks, list cols):
    cdef Matrix_rational_sparse A,B,C,D
    cdef Py_ssize_t nr, nc, i, j
    cdef list pivots, idxs, rest_cols, pivot_rows, pivot_cols
    cdef mpq_vector* r 
    
    nr = AB._nrows
    nc = AB._ncols
    
    pivots = [False] * nc
    idxs = []
    
    for i from 0 <= i < nr:
        r = &(AB._matrix[i])
        # zero row or pivot we already have
        if r.num_nonzero == 0 or pivots[r.positions[0]]: continue
        # new pivot
        j = r.positions[0]
        pivots[j] = True
        idxs.append((j,i))
        
    idxs.sort()
    pivot_cols = [t[0] for t in idxs]
    pivot_rows = [t[1] for t in idxs]
       
    
    
    AB,_ = split_along_rows(AB,pivot_rows)
    
#     rest_cols = [i for i in range(nc) if not pivots[i]]
#     A = AB.matrix_from_rows_and_columns(range(AB._nrows),pivot_cols)
#     B = AB.matrix_from_rows_and_columns(range(AB._nrows),rest_cols)
#     C = CD.matrix_from_rows_and_columns(range(CD._nrows),pivot_cols)
#     D = CD.matrix_from_rows_and_columns(range(CD._nrows),rest_cols)
        
    A,B = split_along_columns(AB, pivot_cols) 
    C,D = split_along_columns(CD, pivot_cols)
            
    blocks = [i - nr for i in blocks if i > nr]
    cols = [cols[i] for i in range(nc) if not pivots[i]]
    
    return A,B,C,D,blocks,cols
############################################################################    
cdef bint vector_equals(mpq_vector* v, mpq_vector* w):
    cdef Py_ssize_t i
    
    if v.num_nonzero != w.num_nonzero: return False
    for i from 0 <= i < v.num_nonzero:
        if v.positions[i] != w.positions[i]: return False
        if not mpq_equal(v.entries[i],w.entries[i]): return False
    
    return True  

###################################################
# rational arithmetic
###################################################

cdef Matrix_rational_sparse _matrix_times_matrix_(Matrix_rational_sparse self, Matrix_rational_sparse right):
    cdef Matrix_rational_sparse ans

    cdef mpq_vector* v
    cdef mpq_vector* w

    # Build a table that gives the nonzero positions in each column of right
    cdef list nonzero_positions_in_columns = [set() for _ in range(right._ncols)]
    cdef Py_ssize_t i, j, k
    
    for i from 0 <= i < right._nrows:
        v = &(right._matrix[i])
        for j from 0 <= j < v.num_nonzero:
            (<set> nonzero_positions_in_columns[v.positions[j]]).add(i)

    ans = self.new_matrix(self._nrows, right._ncols)

    # Now do the multiplication, getting each row completely before filling it in.
    cdef set relevant_columns
    cdef mpq_t x, y, s
    cdef Py_ssize_t pos
    mpq_init(x)
    mpq_init(y)
    mpq_init(s)
    
    for i from 0 <= i < self._nrows:
        v = &(self._matrix[i])
        if not v.num_nonzero: continue
                
        # find all relevant columns
        relevant_columns = set()
        for k from 0 <= k < v.num_nonzero:
            w = &(right._matrix[v.positions[k]])
            (<set> relevant_columns).update({w.positions[p] for p in range(w.num_nonzero)})
                
        # compute dot product with relevant columns
        for j in list(relevant_columns):    
            mpq_set_si(s, 0, 1)
            for k from 0 <= k < v.num_nonzero:
                if v.positions[k] in nonzero_positions_in_columns[j]:
                    mpq_vector_get_entry(y, &right._matrix[v.positions[k]], j)
                    mpq_mul(x, v.entries[k], y)
                    mpq_add(s, s, x)
            mpq_vector_set_entry(&ans._matrix[i], j, s)

    mpq_clear(x)
    mpq_clear(y)
    mpq_clear(s)
    return ans
############################################################################

cdef void trsm(Matrix_rational_sparse A, Matrix_rational_sparse B):
    cdef Py_ssize_t i, k, nc
    cdef mpq_vector* a
    cdef mpq_vector* b
    cdef mpq_vector tmp
    cdef mpq_t a_minus,minus_one
    
    mpq_init(a_minus)
    mpq_init(minus_one)
    mpq_set_si(minus_one,-1,1)
    
    nr = B._nrows

    for i from nr - 1 > i >= 0:
        a = &(A._matrix[i])
        if a.num_nonzero == 1: continue
        for k from 1 <= k < a.num_nonzero:
            mpq_neg(a_minus, a.entries[k])
            add_mpq_vector_init(&tmp, &B._matrix[i], &B._matrix[a.positions[k]], a_minus)
            mpq_vector_clear(&B._matrix[i])
            B._matrix[i] = tmp

    mpq_clear(a_minus)
    mpq_clear(minus_one)

############################################################################

cdef void diff(Matrix_rational_sparse A, Matrix_rational_sparse B):
    cdef Py_ssize_t i,j
    cdef Rational aij, bij, v
    
    for (i,j) in B.nonzero_positions(copy=False,column_order=False):
        aij = A.get_unsafe(i,j)
        bij = B.get_unsafe(i,j)
        v = aij - bij
        A.set_unsafe(i,j,v)

###################################################
# rational echelon computations
###################################################

cdef Matrix_rational_sparse rational_echelon(Matrix_rational_sparse self, bint transformation = False):
        """
        Replace self by its reduction to reduced row echelon form.
        """
        cdef Matrix_rational_sparse H_m, U
        cdef Py_ssize_t i, r, c, nr, nc, pivot
        cdef mpq_t one, a, a_inverse, b, b_minus
        cdef mpq_vector row, row_i, tmp
        
        mpq_init(one)
        mpq_init(a)
        mpq_init(a_inverse)
        mpq_init(b)
        mpq_init(b_minus)
        mpq_set_ui(one, 1, 1)
        
        nr = self._nrows
        nc = self._ncols
        
        # append the identity matrix to obtain the transformation matrix  
        if transformation:
            self = augment(self)
              
        for c from 0 <= c < nc:
            # find pivot row
            pivot = -1
            for r from 0 <= r < nr:
                # find row with pivot
                row = self._matrix[r]
                if row.num_nonzero > 0 and row.positions[0] == c:
                    pivot = r
                    break
            
            if pivot == -1: continue
            sig_on()
            # normalize row    
            row = self._matrix[pivot]
            mpq_set(a,row.entries[0])
            if ~mpq_equal(a,one):
                mpq_inv(a_inverse, a)
                mpq_vector_scale(&row, a_inverse)
            
            # reduce rows below
            for i from pivot < i < nr:
                row_i = self._matrix[i]
                mpq_vector_get_entry(b, &row_i, c)
                if mpq_sgn(b) != 0:
                    mpq_neg(b_minus, b)
                    add_mpq_vector_init(&tmp, &row_i, &row, b_minus)
                    mpq_vector_clear(&row_i)
                    self._matrix[i] = tmp
            
            sig_off()
        
        mpq_clear(one)
        mpq_clear(a)
        mpq_clear(a_inverse)
        mpq_clear(b)
        mpq_clear(b_minus)
            
        return self
############################################################################
cdef Matrix_rational_sparse block_rational_echelon(Matrix_rational_sparse self, list blocks):
        """
        Replace self by its reduction to reduced row echelon form.
        """
        cdef Py_ssize_t i, r, c, nr, nc, pivot, j, next_block 
        cdef mpq_t one, inv, b, b_minus
        cdef mpq_vector tmp
        cdef mpq_vector* row
        cdef mpq_vector* row_i
        
        mpq_init(one)
        mpq_init(inv)
        mpq_init(b)
        mpq_init(b_minus)
        mpq_set_ui(one, 1, 1)
        
        nr = self._nrows
        nc = self._ncols
              
        for c from 0 <= c < nc:          
            # find pivot row
            pivot = -1
            for r from 0 <= r < nr:
                # find row with pivot
                row = &self._matrix[r]
                if row.num_nonzero > 0 and row.positions[0] == c:
                    pivot = r
                    break
            
            if pivot == -1: continue
            sig_on()
            # normalize row    
            row = &self._matrix[pivot]
            if not mpq_equal(row.entries[0],one):
                mpq_inv(inv,row.entries[0])
                mpq_vector_scale(row, inv)
            
            # find index where next block starts
            next_block = blocks[bisect_right(blocks,pivot)]
                        
            # reduce rows below
            for i from next_block <= i < nr:
                row_i = &self._matrix[i]
                mpq_vector_get_entry(b, row_i, c)
                if mpq_sgn(b) != 0:
                    mpq_neg(b_minus, b) 
                    add_mpq_vector_init(&tmp, row_i, row, b_minus)
                    mpq_vector_clear(row_i)
                    self._matrix[i] = tmp
                    
            sig_off()
               
        # normalize rows
        for i from 0 <= i < self._nrows:
            row = &self._matrix[i]
            if not row.num_nonzero: continue
            if not mpq_equal(row.entries[0],one):
                mpq_inv(inv, row.entries[0])
                mpq_vector_scale(row, inv)
        
        mpq_clear(one)
        mpq_clear(inv)
        mpq_clear(b)
        mpq_clear(b_minus)
            
        return self
  
###################################################
# multimodular computations
###################################################

cpdef multimodular_echelon_form(Matrix self):
    cdef Matrix E
    E, pivots = multimodular_matrix_rational_echelon_form_(self)
    E._parent = self.parent()
    return E, pivots
    
############################################################################

cpdef multimodular_matrix_rational_echelon_form_(Matrix self):
    """
    """
    cdef Matrix E
    if self._nrows == 0 or self._ncols == 0:
        return self, ()

    B, _ = self._clear_denom()

    height = self.height()
    height_guess = 10000000*(height+100)

    M = self._ncols * height_guess * height  +  1

    from sage.matrix.matrix_modn_sparse import MAX_MODULUS
    p = MAX_MODULUS + 1
   
    t = None
    X = []
    best_pivots = []
    prod = 1
    problem = 0
    lifts = {}
    while True:
        p = previous_prime(p)
        while prod < M:
            problem = problem + 1

            # We use denoms=False, since we made self integral by calling clear_denom above.
            A = B._mod_int(p)
            modn_echelon_in_place(A)

            # a worthwhile check / shortcut.
            if self._nrows >= self._ncols and self._nrows == len(A.pivots()):
                E = self.parent()(0)
                one = self.base_ring().one()
                for i in range(self._nrows):
                    E.set_unsafe(i, i, one)
                return E, tuple(range(self._nrows))

            c = cmp_pivots(best_pivots, A.pivots())
            if c <= 0:
                best_pivots = A.pivots()
                X.append(A)
                prod = prod * p
            p = previous_prime(p)
        # Find set of best matrices.
        Y = []
        # recompute product, since may drop bad matrices
        prod = 1
        for i in range(len(X)):
            if cmp_pivots(best_pivots, X[i].pivots()) <= 0:
                p = X[i].base_ring().order()
                if p not in lifts:
                    lift = X[i].lift()
                    lifts[p] = (lift, p)
                Y.append(lifts[p])
                prod = prod * X[i].base_ring().order()
        try:
            if not Y:
                raise ValueError("not enough primes")
            a = CRT_basis([w[1] for w in Y])

            # take the linear combination of the lifts of the elements
            # of Y times coefficients in a
            L = a[0]*(Y[0][0])
            assert Y[0][0].is_sparse() == L.is_sparse()
            for j in range(1,len(Y)):
                L += a[j]*(Y[j][0])
            E = L.rational_reconstruction(prod)
            L = 0  # free memory
        except ValueError as msg:
            M = prod * p*p*p
            continue

        d   = E.denominator()
        hdE = int((d*E).height())
        if hdE * self.ncols() * height < prod:
            break
        M = prod * p*p*p
    #end while
    return E, tuple(best_pivots) 
    
############################################################################

cpdef modn_echelon_in_place(Matrix_modn_sparse self):
        """
        Replace self by its reduction to reduced row echelon form.
        """
        self.check_mutability()

        cdef Py_ssize_t i, r, c, min, min_row
        cdef Py_ssize_t a0, a_inverse, b
        cdef c_vector_modint tmp
        pivots = []

        for c from 0 <= c < self._ncols:
            min_row = -1
            for r from 0 <= r < self._nrows:
                if self.rows[r].num_nonzero > 0:
                    if self.rows[r].positions[0] == c:
                        min_row = r
                        break
            if min_row != -1:
                r = min_row
                pivots.append(c)
                a = self.rows[r].entries[0]
                if a != 1:
                    a_inverse = ai.c_inverse_mod_int(a, self.p)
                    scale_c_vector_modint(&self.rows[r], a_inverse)
                sig_on()
                for i from r < i < self._nrows:
                    b = get_entry(&self.rows[i], c)
                    if b != 0:
                        add_c_vector_modint_init(&tmp, &self.rows[i], &self.rows[r], self.p - b)
                        clear_c_vector_modint(&self.rows[i])
                        self.rows[i] = tmp
                sig_off()
        
        self.cache('pivots',tuple(pivots))
        self.cache('in_echelon_form',True)
        