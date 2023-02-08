# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from python_modules import global_data
from time import time


############################################################################
############################################################################
# Ambiguities
############################################################################
############################################################################
cdef class Ambiguity:
    """
    Overlap = fiC - Afj
    Inclusion = fi - AfjC
    """
        
    def __init__(self, str ABC, Py_ssize_t A, Py_ssize_t C, Py_ssize_t i, Py_ssize_t j, bint is_overlap):
        self._ABC = ABC
        self._A = A
        self._C = C
        self._i = i
        self._j = j
        self._deg = len(ABC)
        self._is_overlap = is_overlap
############################################################################     
    def __eq__(self,other):
       return self._i == other._i and \
              self._j == other._j and \
              self._A == other._A and \
              self._C == other._C and \
              self._ABC == other._ABC and \
              self._is_overlap == other._is_overlap
#############################################################################
    def __hash__(self):
       return hash( (self._ABC,self._A,self._C,self._i,self._j,self._is_overlap) )
#############################################################################
    def __repr__(self):
            s = "Overlap" if self._is_overlap else "Inclusion"
            A,C = self.AC()
            return s + "(" + self._ABC + ", " + A + ", " + C + ", (" + str(self._i) + ", " + str(self._j) + "))"
############################################################################     
    cdef tuple AC(Ambiguity self):
        cdef str A,C
        A = self._ABC[:self._A]
        C = self._ABC[self._C:]
        return A,C
############################################################################
    @staticmethod
    cdef list generate(tuple a, tuple b):
        
        cdef str v, w
        cdef Py_ssize_t i, j, k, m
        cdef list amb
        
        i,v = a
        j,w = b
        m = min(len(v),len(w))
        
        amb  = [Ambiguity(v + w[k:],len(v[:-k]),len(v),i,j,True) for k in range(1,m) if v[-k:] == w[:k]]
        amb += [Ambiguity(w + v[k:],len(w[:-k]),len(w),j,i,True) for k in range(1,m) if w[-k:] == v[:k]]

        if len(w) > len(v): amb += Ambiguity.generate_incls(b,a)
        else: amb += Ambiguity.generate_incls(a,b)
            
        return amb
############################################################################    
    @staticmethod
    cdef list generate_incls(tuple a, tuple b):
        cdef str v,w
        cdef Py_ssize_t i,j,k
        cdef list amb
        
        amb = []
        
        i,v = a
        j,w = b
        
        if i == j: return amb
           
        k = v.find(w, 0)
        while k >= 0:
            amb.append( Ambiguity(v,k,k+len(w),i,j,False) )
            k = v.find(w,k+1)

        return amb
############################################################################
    cdef SigPoly to_crit_pair(Ambiguity self, SigPoly gi, SigPoly gj):
        cdef str A,C
        cdef Sig s1,s2
        cdef SigPoly g1,g2
        cdef Py_ssize_t d
        
        A,C = self.AC()

        if not self._is_overlap:
            g1 = gi
            g2 = gj.lrmul(A,C)
        else:
            g1 = gi.rmul(C)
            g2 = gj.lmul(A)
            
        s1 = g1._sig
        s2 = g2._sig
                    
        if s1 == s2:
            g1._deg = -2
            return g1
            
        if g1._poly == g2._poly: 
            d = -1
        else:
            d = self._deg  
               
        if s1 < s2:
            g2._deg = d
            return g2
        else:
            g1._deg = d
            return g1
            
############################################################################
    cdef tuple to_crit_pair_poly(Ambiguity self, NCPoly gi, NCPoly gj):
        cdef str A,C
        cdef NCPoly g1,g2
        cdef Py_ssize_t d
        
        A,C = self.AC()

        if not self._is_overlap:
            g1 = gi
            g2 = gj.lrmul(A,C)
        else:
            g1 = gi.rmul(C)
            g2 = gj.lmul(A)
                 
        d = -1 if g1 == g2 else self._deg  
       
        return d,g1,g2
    
############################################################################
    cdef bint is_redundant(Ambiguity self, str V, Py_ssize_t k):
        if self._is_overlap: return self.overlap_test(V)
        else: return self.incl_test(V,k)
############################################################################        
    cdef bint overlap_test(Ambiguity self, str V):
        """
        a) V | m for some m in {A,B,C}, or
        b) ABC = LVR and
            i. |L| = 0 and |R| < |C|, or
            ii. 0 < |L| < |A|, or
            iii. |L| >= |A| and |R| > 0
        """
        cdef str ABC,A,C
        cdef Py_ssize_t L,R
        
        A,C = self.AC()
        if V in A or V in C or V in self._ABC[self._A : self._C]:
            return True
        ABC = self._ABC

        L = ABC.find(V)
        R = len(ABC) - (L + len(V))
        if L == 0 and R < len(ABC) - self._C: return True
        if L > 0 and L < self._A: return True
        if R > 0 and L >= self._A: return True
        return False
############################################################################
    cdef bint incl_test(Ambiguity self, str V, Py_ssize_t k):
        """
        k < j and 
        a) V | m for some m in {A,C}, or
        b) V | B and |AC| > 0, or
        c) B | V and |V| < |ABC|
        """
        cdef str A,B,C
        
        if k >= self._j: return False
        
        A,C = self.AC()
        if V in A or V in C: return True
        
        B = self._ABC[self._A : self._C]
        
        if V in B and A and C: return True
        if B in V and len(V) < self._deg: return True
        return False