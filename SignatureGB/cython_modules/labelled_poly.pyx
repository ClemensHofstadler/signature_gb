# cython: linetrace=True
# cython: boundscheck=False
# cython: auto_pickle=False

from __future__ import absolute_import

from sage.all import QQ, FreeAlgebra, copy
from sage.rings.rational cimport Rational

import itertools
from cython_modules.ncpoly cimport NCPoly
from cython_modules.sig cimport Sig
from cython_modules.sigpoly cimport SigPoly

from python_modules import global_data

############################################################################
############################################################################
# Labelled polynomial 
############################################################################
############################################################################
cdef class LabelledPoly(SigPoly):

    def __init__(self, NCPoly f, Sig s, Sig pseudo, list mons, list coeffs):
        SigPoly.__init__(self,f,s)
        self._module_mons = mons
        self._module_coeffs = coeffs
        self._pseudo_sig = pseudo
############################################################################
    def __repr__(self):
        p,l = self.to_normal()
        return "poly : " + str(p) + "\nlabel : " + str(l)
############################################################################
    def __len__(self):
        return len(self._module_mons)
############################################################################
    cpdef LabelledPoly copy(LabelledPoly self):
        return LabelledPoly(self._poly.copy(), self._sig.copy(), self._pseudo_sig.copy(), self._module_mons, self._module_coeffs)
############################################################################
    cpdef LabelledPoly lrmul(LabelledPoly self, str a, str b):
        cdef NCPoly p 
        cdef Sig s, ps
        
        p = self._poly.lrmul(a,b)
        s = self._sig.lrmul(a,b)
        ps = self._pseudo_sig.lrmul(a,b)
        return LabelledPoly(p, s, ps, self._module_mons, self._module_coeffs)
############################################################################
    cpdef LabelledPoly module_lrmul(LabelledPoly self, str a, str b):
        cdef Sig s,m,ps
        cdef list mons
        
        s = self._sig.lrmul(a,b)
        ps = self._pseudo_sig.lrmul(a,b)
        mons = [m.lrmul(a,b) for m in self._module_mons]
        return LabelledPoly(self._poly, s, ps, mons, copy(self._module_coeffs))
############################################################################
    cpdef to_normal(LabelledPoly self):
        cdef str s
        cdef Rational c
        cdef Sig m
        
        vars = global_data.vars_ + global_data.module_basis_
        Parent = FreeAlgebra(QQ,len(vars),vars)
        p = self._poly.to_normal(Parent)
        s = ""
        for c,m in zip(self._module_coeffs,self._module_mons):
            s += "".join([str(c), "*", str(m), "+"])
        l = Parent(s + "0")      
        return p,l