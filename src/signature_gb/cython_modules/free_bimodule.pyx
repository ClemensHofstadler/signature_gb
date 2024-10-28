
from __future__ import absolute_import
from sage.all import *
from sage.rings.finite_rings.finite_field_base import FiniteField

from signature_gb.auxiliary import simplify_str
from translator import Translator
    
############################################################################
############################################################################
# FreeAlgebra
############################################################################
############################################################################
cdef class FreeBimodule():
    def __init__(self,F,r,signature_order='DoPoT'):
                
        self._F = F
        self._r = r
        self._is_dpot = (signature_order == 'DoPoT')
        
############################################################################    
    def base_ring(self): return self._F.base_ring()
    def F(self): return self._F
    def gens(self): return [str(x) for x in self._F.gens()]
    def rank(self): return self._r
    def signature_order(self): return 'DoPoT' if self._is_dpot else 'DoToP'
############################################################################      
    def change_ring(self, R):
        self._F.change_ring(R)
############################################################################      
    def __eq__(self,other):
        return self._r == other._r and self._F == other._F
    def __ne__(self,other): return not self.__eq__(other)
############################################################################      
    def __repr__(self):
        s = "Free Bimodule of rank " + str(self._r) + " over " + str(self._F)
        s += " with signature order " + self.signature_order()
        return s
############################################################################
    cpdef bool _cmp(self, Sig s, Sig t): 
        
        cdef object F
        cdef bytes sa, sb, ta, tb
        cdef Py_ssize_t ei, ej, la, lb
        
        if s._len != t._len: return s._len < t._len
        
        if self._is_dpot:  
            ei = s._ei; ej = t._ei
            if ei != ej: return ei < ej  
        
        F = self._F
        sa = s._a; ta = t._a
        if sa != ta: return F._cmp(sa, ta)
            
        sb = s._b; tb = t._b  
        if sb != tb: return F._cmp(sb, tb)
        
        if not self._is_dpot:
            ei = s._ei; ej = t._ei
            if ei != ej: return ei < ej
        
        return False