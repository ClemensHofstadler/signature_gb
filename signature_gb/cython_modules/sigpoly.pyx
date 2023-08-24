# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

############################################################################
############################################################################
# SigPoly
############################################################################
############################################################################
cdef class SigPoly:

    def __init__(self,f,sig):
        if isinstance(f,NCPoly):
            self._poly = f
            self._sig = sig
        else:
            self._poly = NCPoly(f)
            self._sig = Sig(*sig)
############################################################################    
    cpdef SigPoly copy(SigPoly self):
        return SigPoly(self._poly.copy(),self._sig.copy())
############################################################################    
    def __repr__(self):
        return 'poly: ' + repr(self._poly) + ' sig : ' + repr(self._sig) + '\n'
############################################################################
    cdef bool __ceq__(SigPoly self, SigPoly other):
        return self._poly == other._poly and self._sig == other._sig
############################################################################
    def __eq__(self,other):
        return self.__ceq__(other)
############################################################################
    def __mul__(self, other): 
        self._poly.__cmul__(other)
        return self
############################################################################
    cdef Py_ssize_t __chash__(SigPoly self):
        return hash((self._poly,self._sig))
############################################################################
    def __hash__(self):
        return self.__chash__()
############################################################################
    cpdef NCMonomial lm(SigPoly self):
        return self._poly._lm
############################################################################ 
    cpdef list mons(SigPoly self):
        return self._poly._mons
############################################################################
    cdef list coeffs(SigPoly self):
        return self._poly._coeffs
############################################################################
    cpdef SigPoly lmul(SigPoly self, str m):
        return SigPoly(self._poly.lmul(m),self._sig.lmul(m))
############################################################################
    cpdef SigPoly rmul(SigPoly self, str m):
        return SigPoly(self._poly.rmul(m),self._sig.rmul(m))
############################################################################
    cpdef SigPoly lrmul(SigPoly self, str l, str r):
        return SigPoly(self._poly.lrmul(l,r),self._sig.lrmul(l,r))  
