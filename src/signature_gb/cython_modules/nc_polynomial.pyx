# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from collections import defaultdict

from sage.all import GF, QQ
from sage.arith.misc import rational_reconstruction


############################################################################
############################################################################
# NCPoly
############################################################################
############################################################################
cdef class NCPoly:
    def __init__(self,f,*args):
        if args:
            self._mons = args[0]
            self._coeffs = f
            self._lm = self._mons[0]
        else:
            raise NotImplementedError
#             d = f.monomial_coefficients()
#             d = [(d[key],NCMonomial(key)) for key in d]
#             d.sort(key=lambda p : p[1])
#             self._lm = d[0][1]
#             self._mons = [m for c,m in d]
#             self._coeffs = [c for c,m in d]   
############################################################################
    @staticmethod
    cdef inline NCPoly _new() noexcept:
        return <NCPoly>NCPoly.__new__(NCPoly)  
############################################################################
    def coefficients(self): return self._coeffs
    def monomials(self): return self._mons
############################################################################
    cpdef NCPoly copy(self):
        cdef NCMonomial m
        cdef NCPoly f = <NCPoly>NCPoly._new()
        f._coeffs = [c for c in self._coeffs]
        f._mons = [m.copy() for m in self._mons]
        f._lm = f._mons[0]
        return f
############################################################################
    def __repr__(self): return str(self.to_normal())
###########################################################################
    cpdef Py_ssize_t degree(NCPoly self): return len(self._lm._mon)
###########################################################################
    def __eq__(self,other): 
        if len(self._coeffs) != len(other._coeffs): return False
        return  self._coeffs == other._coeffs and \
                self._mons == other._mons
############################################################################
    def __hash__(self): 
        return hash((tuple(self._mons),tuple(self._coeffs)))
############################################################################
    def __mul__(self, other): 
        self._coeffs = [other * c for c in self._coeffs]
############################################################################
    def make_monic(self):
        """
        Really update self
        """
        lc = self._coeffs[0]
        if lc != 1: self._coeffs = [c / lc for c in self._coeffs]
        return lc
############################################################################
    cpdef NCPoly lrmul(NCPoly self, bytes l, bytes r):
        cdef NCPoly f 
        cdef NCMonomial mon
        
        f = <NCPoly>self.copy()
        for mon in f._mons: mon.lrmul(l,r)
        return f
############################################################################
    cdef NCPoly mod(NCPoly self, Py_ssize_t p):
        cdef NCPoly f = <NCPoly>NCPoly._new()
        R = GF(p)
        coeffs = [R(c) for c in self._coeffs]
        f._mons = [m for c,m in zip(coeffs,self._mons) if c != 0]
        f._coeffs = [c for c in coeffs if c != 0]
        f._lm = f._mons[0]
        return f
############################################################################
    def to_normal(NCPoly self):    
        P = self._lm.parent()
        s = ""
        for c,m in zip(self._coeffs,self._mons):
             s += "".join([str(c), "*", str(m), "+"])
        return P._F(s + "0")
###########################################################################
    @staticmethod
    def zero(P): 
        c = [P.zero()]
        m = [NCMonomial('',P)]
        return NCPoly(c,m)
###########################################################################    
    @staticmethod
    cdef NCPoly reconstruct_from_images(list imgs, object B):
        cdef NCPoly g
        cdef Py_ssize_t i
        cdef NCMonomial m
        
        d = defaultdict(lambda : [0] * len(B))
        
        for i,g in enumerate(imgs):
            for c,m in zip(g._coeffs, g._mons):
                d[m][i] = c
            
        poly_data = [(B.crt(c),m) for m,c in d.items()]
        P = B.prod()

        poly_data = [(rational_reconstruction(c,P),m) for c,m in poly_data]
        poly_data = [(c,m) for c,m in poly_data if c != 0]
        poly_data.sort(key=lambda p : p[1], reverse=True)
        
        coeffs,mons = [list(t) for t in zip(*poly_data)]
        return NCPoly(coeffs,mons)
        
            