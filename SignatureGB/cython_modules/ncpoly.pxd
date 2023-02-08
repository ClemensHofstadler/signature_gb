# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from cython_modules.ncmonomial cimport NCMonomial

############################################################################
############################################################################
# NCPoly
############################################################################
############################################################################
cdef class NCPoly:
    cdef public:
        NCMonomial _lm
        list _mons
        list _coeffs
   
    cdef NCPoly copy(self)
    
    cdef NCMonomial lm(NCPoly self)
    
    cdef bool __ceq__(NCPoly self, NCPoly other)
    
    cdef Py_ssize_t __chash__(NCPoly self)
    
    cpdef NCPoly __cmul__(NCPoly self, other)
        
    cpdef make_monic(NCPoly self)
    
    cdef NCPoly lmul(NCPoly self, str m)
    
    cdef rmul(self,m)
    
    cpdef NCPoly lrmul(NCPoly self, str l, str r)
    
    cpdef to_normal(NCPoly self, Parent=*)