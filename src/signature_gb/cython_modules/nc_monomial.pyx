# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import
from sage.all import ZZ
from cpython.object cimport (Py_EQ, Py_NE, Py_LT, Py_LE, Py_GT, Py_GE)

import re
        
############################################################################
############################################################################
# NCMonomial
############################################################################
############################################################################
cdef class NCMonomial:
    def __init__(self,monomial,P):
        cdef bytes m
        cdef Py_ssize_t j
        
        self._parent = P
        
        if isinstance(monomial,str): 
            self._mon = monomial.encode()
        elif str(monomial) == '1': 
            self._mon = b''
        else:
            m = b''
            for v in str(monomial).split('*'):
                 j = v.find('^')
                 if j != -1:
                     m += (v[:j] * ZZ(v[j+1:])).encode()
                 else:
                     m += v.encode()
            self._mon =  m
############################################################################
    cpdef bytes mon(NCMonomial self): return self._mon
    def parent(NCMonomial self): return self._parent
############################################################################        
    def change_parent(self, P): 
        self._parent = P
############################################################################    
    cpdef NCMonomial copy(NCMonomial self): 
        cdef NCMonomial m = <NCMonomial>NCMonomial.__new__(NCMonomial)
        m._mon = <bytes>self._mon
        m._parent = self._parent
        return m
############################################################################
    def __repr__(self):
        P = self._parent
        T = P.translator()        
        if not self._mon: return '1'
        
        m = self._mon.decode()
        m = "*" + "*".join([T(c,to_internal=False) for c in m]) + "*"
        for x in P.gens():
            eq = r"((?<=\*)%s\*){2,}" % x
            l = len(x) + 1
            m = re.sub(eq, lambda y : str(x) + "^" + str((y.end() - y.start())//l) + "*", m)
        return m[1:-1]
############################################################################
    def __richcmp__(NCMonomial self, NCMonomial other, int op):
        
        if op == Py_EQ:
            return self._mon == other._mon
        elif op == Py_LT:
            return self._parent._cmp(self._mon, other._mon)
        elif op == Py_NE:
            return self._mon != other._mon
        elif op == Py_LE:
            return self._mon == other._mon or self._parent._cmp(self._mon, other._mon)
        elif op == Py_GT:
            return other < self
        elif op == Py_GE:
            return other <= self
############################################################################
    def __hash__(self): return self._mon.__hash__()
############################################################################
    cdef void lrmul(NCMonomial self, bytes l, bytes r): 
        self._mon = <bytes>(l + self._mon + r)
############################################################################
    def to_normal(NCMonomial self):
        P = self._parent
        return P._F(str(self))