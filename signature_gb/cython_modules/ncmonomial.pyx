# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import
from sage.all import ZZ, QQ, FreeAlgebra

from python_modules import global_data

from time import time
        

############################################################################
############################################################################
# NCMonomial
############################################################################
############################################################################
cdef class NCMonomial:
    def __init__(self,monomial):
        cdef str m
        cdef Py_ssize_t j
        
        if isinstance(monomial,str): self._mon = <str>monomial
        elif str(monomial) == '1': self._mon = ''
        else:
            m = ''
            for v in str(monomial).split('*'):
                 j = v.find('^')
                 if j != -1:
                     m += v[:j] * ZZ(v[j+1:])
                 else:
                     m += v
            self._mon =  m
############################################################################
    cpdef NCMonomial copy(NCMonomial self): return NCMonomial(self._mon)
############################################################################
    cdef str mon(NCMonomial self): return self._mon
############################################################################
    def __lt__(self,other): 
        # this is a strict order
        if self == other: return False

        return global_data.cmp_mon_(self._mon,other._mon)
############################################################################
    def __repr__(self):
        s = '*'.join([self._mon[i:i+1] for i in range(len(self._mon))])
        return s if s else '1'
############################################################################
    cdef bool __ceq__(NCMonomial self, NCMonomial other): return self._mon == other._mon
############################################################################
    def __eq__(self,other): return self.__ceq__(other)
############################################################################
    cdef Py_ssize_t __chash__(NCMonomial self): return hash(self._mon)
############################################################################
    def __hash__(self): return self.__chash__()
############################################################################
    cdef void lmul(NCMonomial self, str s): self._mon = s + self._mon
############################################################################
    cdef void rmul(NCMonomial self, str s): self._mon += s
############################################################################
    cdef void lrmul(NCMonomial self, str l, str r): self._mon = <str>''.join([l, self._mon, r])
############################################################################
    cpdef to_normal(NCMonomial self,Parent=None):
        if not Parent:
            Parent = FreeAlgebra(QQ,len(global_data.__vars__), global_data.__vars__)
        s = '*'.join([self._mon[i:i+1] for i in range(len(self._mon))])
        if not s: s = '1'
        return Parent(s)