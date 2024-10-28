# cython: linetrace=True

from __future__ import absolute_import

from cpython cimport *

cdef class Sig:
    cdef readonly:
        bytes _a,_b
        Py_ssize_t _ei, _len
        object _parent
                    
    cdef Sig copy(Sig self)
    
    cpdef tuple aib(Sig self)
    
    cpdef bytes to_bytes(Sig self)
                        
    cdef bool __ceq__(Sig self, Sig other)
    
    cpdef Sig lrmul(Sig self, bytes l, bytes r)
    
    cpdef tuple divides(Sig self, Sig other)
    
    cpdef to_poly(Sig self, list F)