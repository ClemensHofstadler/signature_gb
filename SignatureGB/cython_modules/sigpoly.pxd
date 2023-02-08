# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from cython_modules.ncmonomial cimport NCMonomial
from cython_modules.ncpoly cimport NCPoly
from cython_modules.sig cimport Sig

############################################################################
############################################################################
# SigPoly
############################################################################
############################################################################
cdef class SigPoly:

    cdef readonly:
        NCPoly _poly
        Sig _sig
    cdef public int _deg
    
    cpdef SigPoly copy(SigPoly self)
    
    cdef bool __ceq__(SigPoly self, SigPoly other)
    
    cdef Py_ssize_t __chash__(SigPoly self)
    
    cpdef NCMonomial lm(SigPoly self)
    
    cpdef list mons(SigPoly self)
    
    cdef list coeffs(SigPoly self)
    
    cpdef SigPoly lmul(SigPoly self, str m)
    
    cpdef SigPoly rmul(SigPoly self, str m)
    
    cpdef SigPoly lrmul(SigPoly self, str l, str r)