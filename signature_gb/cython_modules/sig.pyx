# cython: linetrace=True
# cython: boundscheck=False
# cython: auto_pickle=False

from __future__ import absolute_import

from python_modules import global_data
from sage.all import QQ, FreeAlgebra

from time import time

############################################################################
############################################################################
# Sig
############################################################################
############################################################################
cdef class Sig(object):
    def __init__(self, str a, int ei, str b):
        self._a = a
        self._b = b
        self._ei = ei
        if 0 < ei < len(global_data.weight):
            self._len = len(a) + len(b) + global_data.weight[ei]
        else: 
            self._len = 0
        self._hash = hash((a,ei,b))
############################################################################
    cdef Sig copy(Sig self): return Sig(self._a,self._ei,self._b)
############################################################################
    cpdef tuple aib(Sig self): return self._a,self._ei,self._b
############################################################################
    cdef str a(Sig self): return self._a
############################################################################
    def __len__(self): return self._len
############################################################################    
    cdef Py_ssize_t ei(Sig self): return self._ei
############################################################################
    cdef str b(Sig self): return self._b
############################################################################
    def __repr__(self): 
        a = '*'.join([self._a[i:i+1] for i in range(len(self._a))])
        b = '*'.join([self._b[i:i+1] for i in range(len(self._b))])
        if a: a += '*'
        if b: b = '*' + b
        return ''.join([a, "e", str(self._ei), b])
############################################################################
    def __hash__(self): return self._hash
############################################################################
    cdef bint __ceq__(Sig self, Sig other): return <bint>(self._ei == other._ei and self._a == other._a and self._b == other._b)
############################################################################
    def __eq__(self, other): return self.__ceq__(other)
############################################################################
    def __lt__(self, other): 
        # this is a strict order
        if self == other: return False
            
        # handle signature 0
        if not self._ei or not other._ei: return self._ei == 0
        
        # standard case
        return global_data.cmp_sig_(self,other)
############################################################################
    def my_str(self): return ''.join([self._a + "*", str(self._ei), "*" + self._b])
############################################################################
    cdef Sig lmul(Sig self, str m):
        if not self._ei: return self
        return <Sig>Sig(m + self._a, self._ei, self._b) 
############################################################################
    cdef Sig rmul(Sig self, str m):
        if not self._ei: return self
        return <Sig>Sig(self._a, self._ei, self._b + m) 
############################################################################
    cpdef Sig lrmul(Sig self, str l, str r):
        if not self._ei: return self
        return <Sig>Sig(l + self._a, self._ei, self._b + r) 
############################################################################
    cpdef tuple divides(Sig self, Sig other):
        cdef str sa, sb
        cdef Py_ssize_t i,j
        
        if self._ei != other._ei: return None
        
        sa = self._a
        sb = self._b
        
        if len(other._a) < len(sa) or len(other._b) < len(sb): return None
        
        i = len(other._a) - len(sa)
        j = len(sb)
        if other._a[i:] == sa and other._b[:j] == sb:
            return other._a[:i] , other._b[j:]
        return None
############################################################################
    cpdef to_normal(Sig self, Parent=None):
        if not Parent:
            vars = global_data.__vars__ + global_data.__module_basis__
            Parent = FreeAlgebra(QQ,len(vars),vars)
        
        a = '*'.join([self._a[i:i+1] for i in range(len(self._a))])
        b = '*'.join([self._b[i:i+1] for i in range(len(self._b))])
        if a: a += '*'
        if b: b = '*' + b
        return Parent(a + "e" + str(self._ei) + b)