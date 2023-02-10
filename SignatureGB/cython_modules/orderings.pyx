# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from python_modules import global_data

###########################################################################
############################################################################
# Orderings 
############################################################################
############################################################################
cdef bint deglex(str a, str b):
    cdef Py_ssize_t la, lb, i
    cdef str x,y
    cdef dict d
        
    la = len(a)
    lb = len(b)
    
    if la != lb: return la < lb
    
    d = <dict>global_data.vars_dict_
    
    for i from 0 <= i < la:
        if <str>a[i] != <str>b[i]:
            x = <str>a[i]
            y = <str>b[i]
            return d[x] < d[y]
############################################################################
cdef bint multilex(str a, str b):
    cdef Py_ssize_t Xa, Xb

    for block in global_data.var_blocks_:
        Xa = sum(x in block for x in a)
        Xb = sum(x in block for x in b)
        if Xa != Xb: return Xa < Xb

    return deglex(a,b)
   
############################################################################
cdef bint test_order(str a, str b):
    
    # compare overall degree
    if len(a) != len(b): return len(a) < len(b)
    
    # compare mons with t = 1
    a_wo_t = a.replace('t','')
    b_wo_t = b.replace('t','')
    
    if a_wo_t != b_wo_t: return deglex(a_wo_t,b_wo_t)
    
    # use lex where t is largest
    return deglex(a,b)
          
############################################################################
cdef bint deg_term_over_position(Sig s, Sig t):
    
    cdef str sa, sb, ta, tb
    cdef Py_ssize_t ei, ej, la, lb
        
    if s._len != t._len: return s._len < t._len
    
    sa = s._a
    ei = s._ei
    sb = s._b
    
    ta = t._a
    ej= t._ei
    tb = t._b
  
    if sa != ta: return global_data.cmp_mon_(sa,ta)
    if sb != tb: return global_data.cmp_mon_(sb,tb)
    return ei < ej
  
############################################################################
cdef bint deg_position_over_term(Sig s, Sig t):
    
    cdef str sa, sb, ta, tb
    cdef Py_ssize_t ei, ej, la, lb
        
    if s._len != t._len: return s._len < t._len
      
    ei = s._ei; ej = t._ei
    if ei != ej: return ei < ej  
    
    sa = s._a; ta = t._a
    if sa != ta: return global_data.cmp_mon_(sa,ta)
    
    sb = s._b; tb = t._b  
    return global_data.cmp_mon_(sb,tb)
    
############################################################################
cdef bint F5_order(Sig s, Sig t):
    
    cdef str sa, sb, ta, tb
    cdef Py_ssize_t ei, ej, la, lb
    
    ei = s._ei; ej = t._ei

    s_len = s._len + global_data.weight[ei]
    t_len = t._len + global_data.weight[ej]
        
    if s_len != t_len: return s_len < t_len
      
    if ei != ej: return ei < ej  
    
    sa = s._a; ta = t._a
    if sa != ta: return global_data.cmp_mon_(sa,ta)
    
    sb = s._b; tb = t._b  
    return global_data.cmp_mon_(sb,tb)