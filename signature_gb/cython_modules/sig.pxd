# cython: linetrace=True

from __future__ import absolute_import

from cpython cimport *

cdef class Sig:
    cdef readonly:
        str _a,_b
        Py_ssize_t _ei, _len, _hash
            
    cdef Sig copy(Sig self)
    
    cpdef tuple aib(Sig self)
    
    cdef str a(Sig self)
    
    cdef Py_ssize_t ei(Sig self)
    
    cdef str b(Sig self)
        
    cdef bint __ceq__(Sig self, Sig other)
        
    cdef Sig lmul(Sig self, str m)
    
    cdef Sig rmul(Sig self, str m)
    
    cpdef Sig lrmul(Sig self, str l, str r)
    
    cpdef tuple divides(Sig self, Sig other)
    
    cpdef to_normal(Sig self, Parent=*)