# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

###########################################################################
############################################################################
# Ambiguities
############################################################################
############################################################################
cdef class MyFreeAlgebra:
    cdef public:
        object _F, _translator
        list _gens, _blocks
        bool _is_block_order
    
    cpdef bool _cmp(self, bytes a, bytes b)
        