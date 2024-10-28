# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from .sig cimport Sig


from cpython cimport *

###########################################################################
############################################################################
# Ambiguities
############################################################################
############################################################################
cdef class FreeBimodule:
    cdef public:
        object _F
        Py_ssize_t _r
        bool _is_dpot
    
    cpdef bool _cmp(self, Sig a, Sig b)
        