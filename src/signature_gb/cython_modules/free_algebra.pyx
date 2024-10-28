
from __future__ import absolute_import
from sage.all import *
from sage.rings.finite_rings.finite_field_base import FiniteField

from signature_gb.auxiliary import simplify_str
from translator import Translator

import sage.parallel.multiprocessing_sage as MP
    
############################################################################
############################################################################
# FreeAlgebra
############################################################################
############################################################################
cdef class MyFreeAlgebra():
    def __init__(self,K,X):
                
        if K not in {QQ} and not isinstance(K, FiniteField):
            raise TypeError("Coefficient field %s not supported" % str(K))
        if not len(X): raise ValueError("Need at least one variable")
        
        T = Translator(X)
        self._translator = T       
        X = self.set_order(X) 
        self._F = FreeAlgebra(K,X)
        self._gens = [str(x) for x in self._F.gens()]
        
############################################################################    
    def base_ring(self): return self._F.base_ring()
    def F(self): return self._F
    def gens(self): return self._gens
    def is_block_order(self): return self._is_block_order
    def blocks(self): return self._blocks
    def translator(self): return self._translator
############################################################################      
    def change_ring(self, R):
        self._F = FreeAlgebra(R,self._gens)
############################################################################      
    def __eq__(self,other):
        return self._F == other._F and \
            self._is_block_order == other._is_block_order and \
            self._blocks == other._blocks
    def __ne__(self,other): return not self.__eq__(other)
############################################################################      
    def __repr__(self):
        s = str(self._F) + " with "
        if not self._is_block_order:
            for x in self.gens():
                s += x + " < "
        else: 
            idx = -1
            T = self._translator
            block = {T(chr(x),to_internal=False) for x in self._blocks[idx]}
            for x in self.gens():
                if x in block:
                    s += x + " < "
                else:
                    s = s[:-3] + " << " + x + " < "
                    idx -= 1
                    block = {T(chr(x),to_internal=False) for x in self._blocks[idx]}
        return s[:-3]
############################################################################
    def __call__(self,f):
        
        from .nc_monomial import NCMonomial
        from .nc_polynomial import NCPoly
                
        F = self._F
        # from NCPolynomial to FreeAlgebra
        if isinstance(f,NCPoly):
            s = ""
            for c,m in zip(f.coefficients(),f.monomials()):
                s += str(c) + "*" + str(m) + "+"
            return F(s + "0")
        
        # from string to NCPolynomial
        elif isinstance(f,str):
            return self(F(f))
        # from FreeAlgebra to NCPolynomial
        else:
            T = self._translator
            d = f.monomial_coefficients()
            d = [(d[key],NCMonomial(T(simplify_str(str(key))),self)) for key in d]
            d.sort(key=lambda p : p[1],reverse=True)
            mons = [m for c,m in d]
            coeffs = [c for c,m in d]
            if not mons:
                return NCPoly.zero(self)
            return NCPoly(coeffs,mons)
############################################################################   
    def set_order(self,X):
            
        if isinstance(X[0],list):
            T = self._translator
            self._is_block_order = True
            self._blocks = list(reversed([{T(str(x)).encode()[0] for x in block} for block in X if block]))
            return flatten(X)
        else:
            self._is_block_order = False
            self._blocks = None
            return X
        
############################################################################
    cpdef bool _cmp(self, bytes a, bytes b): 
        cdef Py_ssize_t aa, bb
        
        if <bool>self._is_block_order:
            for block in <list>self._blocks:
                aa = sum([x in block for x in a])
                bb = sum([x in block for x in b])
                if aa != bb: return aa < bb
        
        aa = len(a); bb = len(b)
        if aa != bb: return aa < bb

        return <bool>(a < b)
############################################################################   
    def zero(self):
        return self.base_ring().zero() 