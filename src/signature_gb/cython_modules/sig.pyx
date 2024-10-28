# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython.object cimport (Py_EQ, Py_NE, Py_LT, Py_LE, Py_GT, Py_GE)

cdef bytes SIG_SEPARATOR = b'\x7f' #  = chr(127)

############################################################################
############################################################################
# Sig
############################################################################
############################################################################
cdef class Sig:

    def __init__(self, bytes a, Py_ssize_t ei, bytes b, object P):
        self._a = a
        self._b = b
        self._ei = ei
        self._len = 0
        self._parent = P
############################################################################
    cdef Sig copy(Sig self):
        cdef Sig s = Sig.__new__(Sig)
        
        s._a = <bytes>self._a
        s._b = <bytes>self._b
        s._ei = <Py_ssize_t>self._ei
        s._len = <Py_ssize_t>self._len
        s._parent = self._parent
        return s
############################################################################
    cpdef tuple aib(Sig self): return self._a,self._ei,self._b
############################################################################
    def __len__(self): return self._len
############################################################################
    def set_len(self, l): self._len = l
############################################################################
    def __repr__(self): 
        a = '*'.join([str(c) for c in self._a])
        b = '*'.join([str(c) for c in self._b])
        if a: a += '*'
        if b: b = '*' + b
        return ''.join([a, "e", str(self._ei), b])
############################################################################
    def __hash__(self): return hash((self._a,self._ei,self._b))
############################################################################
    def __richcmp__(Sig self, Sig other, int op):
        
        if op == Py_EQ:
            return self.__ceq__(other)
        elif op == Py_LT:
            return self._parent._cmp(self, other)
        elif op == Py_NE:
            return not self.__ceq__(other)
        elif op == Py_LE:
            return self.__ceq__(other) or self._parent._cmp(self, other)
        elif op == Py_GT:
            return other < self
        elif op == Py_GE:
            return other <= self
############################################################################
    cdef bool __ceq__(Sig self, Sig other):
        return self._ei == other._ei and \
                self._a == other._a and \
                self._b == other._b
############################################################################
    cpdef bytes to_bytes(Sig self):
        return <bytes>(self._a + SIG_SEPARATOR + chr(self._ei).encode() + SIG_SEPARATOR + self._b)
############################################################################
    cpdef Sig lrmul(Sig self, bytes l, bytes r):
        cdef Sig s = Sig.__new__(Sig)
        
        s._a = <bytes>(l + self._a)
        s._ei = <Py_ssize_t>self._ei
        s._b = <bytes>(self._b + r)
        s._len = <Py_ssize_t>(self._len + len(l) + len(r))
        s._parent = self._parent
        return s
############################################################################
    cpdef tuple divides(Sig self, Sig other):
        cdef bytes sa, sb
        cdef Py_ssize_t i,j
        
        if self._ei != other._ei: return None
        
        sa = <bytes>self._a
        sb = <bytes>self._b
                
        if len(other._a) < len(sa) or len(other._b) < len(sb): return None
        
        i = len(other._a) - len(sa)
        j = len(sb)
        if other._a[i:] == sa and other._b[:j] == sb:
            return other._a[:i] , other._b[j:]
        
        return None
############################################################################
    cpdef to_poly(Sig self, list F):
        cdef bytes a,b
        cdef Py_ssize_t i
        
        i = self._ei - 1
        a = self._a; b = self._b
        return F[i].lrmul(a,b)

############################################################################
    def to_normal(Sig self, P):
        
        a = '*'.join([chr(c) for c in self._a])
        b = '*'.join([chr(c) for c in self._b])
        if a: a += '*'
        if b: b = '*' + b
        return P(a + "e" + str(self._ei) + b)
