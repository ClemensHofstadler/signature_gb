# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from cython_modules.ncpoly cimport NCPoly
from cython_modules.sig cimport Sig
from cython_modules.sigpoly cimport SigPoly

############################################################################
############################################################################
# Labelled polynomial 
############################################################################
############################################################################
cdef class LabelledPoly(SigPoly):
    cdef public:
        list _module_coeffs, _module_mons
        Sig _pseudo_sig
  
    cpdef LabelledPoly copy(LabelledPoly self)
    
    cpdef LabelledPoly lrmul(LabelledPoly self, str a, str b)
    
    cpdef LabelledPoly module_lrmul(LabelledPoly self, str a, str b)
    
    cpdef to_normal(LabelledPoly self)
        