# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from signature_gb.cython_modules.linear_algebra cimport *

############################################################################
############################################################################
# SigPoly
############################################################################
############################################################################
cdef class SigPoly:

    cdef SigPoly _new(self) noexcept:
        cdef SigPoly s = <SigPoly>SigPoly.__new__(SigPoly)
        return s

    def __init__(self,f,sig):
        if isinstance(f,NCPoly):
            self._poly = f
            self._sig = sig
            self._deg = f.degree()
        else:
            raise NotImplementedError("Input conversion not implemented")
############################################################################    
    cpdef Py_ssize_t degree(self): return self._poly.degree()
############################################################################    
    cpdef SigPoly copy(SigPoly self):
        cdef SigPoly s = self._new()
        s._poly = <NCPoly>self._poly.copy()
        s._sig = <Sig>self._sig.copy()
        s._deg = <Py_ssize_t>self._deg
        return s
############################################################################    
    def __repr__(self):
        return 'poly: ' + str(self._poly) + ', sig : ' + self._sig.__repr__()
############################################################################
    def __eq__(self,other):
        return self._poly == other._poly and self._sig == other._sig
############################################################################
    def __mul__(self, other): 
        self._poly.__mul__(other)
        return self
############################################################################
    def __hash__(self):
        return hash((self._poly,self._sig))
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
    cpdef SigPoly lrmul(SigPoly self, bytes l, bytes r):
        cdef SigPoly s = self._new()
        s._poly = <NCPoly>self._poly.lrmul(l,r)
        s._sig = <Sig>self._sig.lrmul(l,r)
        s._deg = self._deg + len(l) + len(r)
        return s 
############################################################################
    cdef Py_ssize_t deficiency(SigPoly self):
        return self._sig._len - self._deg
############################################################################
    cpdef SigPoly mod(SigPoly self, Py_ssize_t p):
        cdef SigPoly s = self._new()
        s._poly = <NCPoly>self._poly.mod(p)
        s._sig = <Sig>self._sig.copy()
        s._deg = <Py_ssize_t>self._deg
        return s 
############################################################################
    def diff_in_place(self, other):
        
        res = self.copy()
        
        d = {m : c for m,c in zip(res._poly._mons, res._poly._coeffs)}
        for m,c in zip(other._poly._mons, other._poly._coeffs):
            if m in d:
                d[m] -= c
            else:
                d[m] = -c
        
        d = [(c,m) for m,c in d.items() if c]
        d.sort(key=lambda p: p[1],reverse=True)
        res._poly._mons = [m for c,m in d]
        res._poly._coeffs = [c for c,m in d]
        
        return res
        
############################################################################
                    
    @staticmethod
    cdef SigPoly reconstruct_from_images(Sig s, list imgs, object B):
        cdef NCPoly f = NCPoly.reconstruct_from_images(imgs, B)
        return SigPoly(f,s)