# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from .nc_monomial cimport NCMonomial
from .nc_polynomial cimport NCPoly
from .sig cimport Sig

############################################################################
############################################################################
# SigPoly
############################################################################
############################################################################
cdef class SigPoly:

    cdef public:
        NCPoly _poly
        Sig _sig
        Py_ssize_t _deg
    
    cdef SigPoly _new(SigPoly self) noexcept
    
    cpdef Py_ssize_t degree(self)
    
    cpdef SigPoly copy(SigPoly self)
            
    cpdef NCMonomial lm(SigPoly self)
    
    cpdef list mons(SigPoly self)
    
    cdef list coeffs(SigPoly self)
    
    cpdef SigPoly lrmul(SigPoly self, bytes l, bytes r)
    
    cdef Py_ssize_t deficiency(SigPoly self)
    
    cpdef SigPoly mod(SigPoly self, Py_ssize_t p)
        
    @staticmethod
    cdef SigPoly reconstruct_from_images(Sig s, list imgs, object B)