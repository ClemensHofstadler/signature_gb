from __future__ import absolute_import

from sage.all import QQ, matrix
from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse
from sage.rings.rational cimport Rational

from cython_modules.ncmonomial cimport NCMonomial
from cython_modules.ncpoly cimport NCPoly
from cython_modules.linear_algebra cimport augment

############################################################################
############################################################################
# Normal Form
############################################################################
############################################################################
cpdef list interreduce(list F):
    cdef Py_ssize_t i,s
    cdef list G
    cdef NCPoly g
    
    G = [g.copy() for g in F]
    
    i = 0
    s = len(G)
    while i < s:
        if not G[i]:
            i += 1
            continue

        normal_form = reduced_form(G,i)
        # normal form = 0
        if type(normal_form) is int:
            G[i] = None
        # normal form != 0
        elif normal_form:
            G[i] = normal_form
            i = 0
        # normal form = g_i
        else: i += 1

    return [g for g in G if g]
############################################################################    
cpdef reduced_form(list G, Py_ssize_t i, bint intern=True):
    cdef NCPoly f, normal_form
    cdef list F,R,pos
    cdef Matrix_rational_sparse M,T
    cdef Py_ssize_t row_idx,j
    
    f = G[i]
    F = [f]
    R, cols = symbolic_preprocessing_red(F,G,i)
    cols = list(cols)
    cols.sort(reverse=True)
    
    # no reducers found
    if not R:
        if intern: return None
        else: return f

    M = set_up_matrix(R+F,cols)
    M = augment(M)
    M = M.echelon_form()
    T = M[:,len(cols):]
    M = M[:,:len(cols)]

    row_idx = T.nonzero_positions_in_column(0)[-1]
    pos = M.nonzero_positions_in_row(row_idx)

    # reduction to zero
    if len(pos) == 0:
        if intern: return 0
        else: normal_form = NCPoly.zero()
    else:
        coeffs = [M.get_unsafe(row_idx,j) for j in reversed(pos)]
        mons = [cols[j] for j in reversed(pos)]
        normal_form = NCPoly(mons,coeffs)
    
    return normal_form
############################################################################      
cdef tuple symbolic_preprocessing_red(list P, list G, Py_ssize_t idx):
    cdef Py_ssize_t i
    cdef set done, todo
    cdef NCMonomial m, t
    cdef NCPoly g,agb
    cdef list R, lm
    cdef str t_str, a, b
     
    R = []           
    done = set()
    todo = {m for f in P for m in f._mons}
    
    lm = [(i,g._lm) for i,g in enumerate(G) if i != idx and g]
        
    while todo:
        t = todo.pop()
        done.add(t)
        t_str = t._mon
        
        i,m = next(((i,m) for i,m in lm if m._mon in t_str), (0,None))                                      
        if m:
            g = G[i]
            a,b = t_str.split(m._mon,1)
            agb = g.lrmul(a,b)
            R.append(agb)
            todo.update({m for m in agb._mons if m not in done})
        
    return R, done 
############################################################################
cdef Matrix_rational_sparse set_up_matrix(list rows, list columns):
    cdef Matrix_rational_sparse A
    cdef Py_ssize_t nr, nc, i, j
    cdef Rational c
    cdef dict cols
    cdef NCMonomial m
    cdef NCPoly f
   
    nr = len(rows)
    nc = len(columns)
    A = matrix(QQ,nr,nc,sparse=True)
    
    cols = {m:i for i,m in enumerate(columns)}
                
    for i,f in enumerate(rows):
        for c,m in zip(f._coeffs, f._mons):
            j = cols[m]
            A.set_unsafe(i,j,c)
    
    return A   
