# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from .nc_monomial cimport NCMonomial

############################################################################
############################################################################
# NCPoly
############################################################################
############################################################################
cdef class NCPoly:
    cdef public:
        NCMonomial _lm
        list _mons, _coeffs
   
    @staticmethod
    cdef inline NCPoly _new() noexcept
    
    cpdef NCPoly copy(self)
    
    cpdef Py_ssize_t degree(NCPoly self)
                        
    cpdef NCPoly lrmul(NCPoly self, bytes l, bytes r)
    
    cdef NCPoly mod(NCPoly self, Py_ssize_t p)
    
    @staticmethod
    cdef NCPoly reconstruct_from_images(list imgs, object B)