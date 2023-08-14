# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

############################################################################
############################################################################
# NCMonomial
############################################################################
############################################################################
cdef class NCMonomial:
    cdef public str _mon
  
    cpdef NCMonomial copy(NCMonomial self)
    
    cdef str mon(NCMonomial self)
    
    cdef bool __ceq__(NCMonomial self, NCMonomial other)
    
    cdef Py_ssize_t __chash__(NCMonomial self)
    
    cdef void lmul(NCMonomial self, str s)
    
    cdef void rmul(NCMonomial self, str s)
    
    cdef void lrmul(NCMonomial self, str l, str r)
    
    cpdef to_normal(NCMonomial self,Parent=*)