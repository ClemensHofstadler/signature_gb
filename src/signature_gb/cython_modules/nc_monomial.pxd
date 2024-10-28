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
    cdef readonly:
        bytes _mon
        object _parent
            
    cpdef bytes mon(NCMonomial self) 
  
    cpdef NCMonomial copy(NCMonomial self)
    
    cdef void lrmul(NCMonomial self, bytes l, bytes r)