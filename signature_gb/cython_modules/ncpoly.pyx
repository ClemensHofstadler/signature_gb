# cython: linetrace=True
# cython: boundscheck=False
# cython: auto_pickle=False

from __future__ import absolute_import

from sage.all import ZZ, QQ, FreeAlgebra, copy
from sage.rings.rational cimport Rational

from cython_modules.ncmonomial cimport NCMonomial
from cython_modules.orderings cimport *

from python_modules import global_data

############################################################################
############################################################################
# NCPoly
############################################################################
############################################################################
cdef class NCPoly:

    def __init__(self,f,*args):
        if args:
            self._lm = f[-1]
            self._mons = f
            self._coeffs = args[0]
        else:
            d = f.monomial_coefficients()
            d = [(d[key],NCMonomial(key)) for key in d]
            d.sort(key=lambda p : p[1])
            self._lm = d[-1][1]
            self._mons = [m for c,m in d]
            self._coeffs = [c for c,m in d]      
############################################################################
    cdef NCPoly copy(self):
        cdef NCMonomial m
        cdef Rational c
        return NCPoly([m.copy() for m in self._mons],[c for c in self._coeffs])
############################################################################
    def __repr__(self): return str(self.to_normal())
###########################################################################
    def degree(self): return len(self._lm._mon)
###########################################################################
    cdef NCMonomial lm(NCPoly self): return self._lm
###########################################################################
    @staticmethod
    def zero(): return NCPoly([NCMonomial('')],[QQ.zero()])
###########################################################################
    cdef bool __ceq__(NCPoly self, NCPoly other):
        if len(self._coeffs) != len(other._coeffs): return False
        return self._lm == other._lm and self._coeffs == other._coeffs and self._mons == other._mons
############################################################################
    def __eq__(self,other): return self.__ceq__(other)
############################################################################
    cdef Py_ssize_t __chash__(NCPoly self):
        return hash((self._lm,tuple(self._mons),tuple(self._coeffs)))
############################################################################
    def __hash__(self): return self.__chash__()
############################################################################
    cpdef NCPoly __cmul__(NCPoly self, other):
        self._coeffs = [other * c for c in self._coeffs]
############################################################################
    def __mul__(self, other): return self.__cmul__(other)
############################################################################
    cpdef make_monic(NCPoly self):
        """
        Really update self
        """
        lc = self._coeffs[-1]
        if lc != 1: self._coeffs = [c / lc for c in self._coeffs]
        return lc
############################################################################
    cdef NCPoly lmul(NCPoly self, str m):
        """
        Only update copy
        """
        cdef NCPoly f
        cdef NCMonomial mon
        
        f = self.copy()
        for mon in f._mons: mon.lmul(m)
        return f
############################################################################
    cdef rmul(self,m):
        """
        Only update copy
        """
        cdef NCPoly f 
        cdef NCMonomial mon
        
        f = self.copy()
        for mon in f._mons: mon.rmul(m)
        return f
############################################################################
    cpdef NCPoly lrmul(NCPoly self, str l, str r):
        cdef NCPoly f 
        cdef NCMonomial mon
        
        f = self.copy()
        for mon in f._mons: mon.lrmul(l,r)
        return f
############################################################################
    cpdef to_normal(NCPoly self, Parent=None):
        cdef str s
        cdef Rational c
        cdef NCMonomial m
        
        if not Parent:
            Parent = FreeAlgebra(QQ,len(global_data.vars_), global_data.vars_)
        s = ""
        for c,m in zip(self._coeffs,self._mons):
             s += "".join([str(c), "*", str(m), "+"])
        return Parent(s + "0")