# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from .nc_monomial cimport NCMonomial
from .nc_polynomial cimport NCPoly
from .sig cimport Sig
from .sig_polynomial cimport SigPoly

###########################################################################
############################################################################
# Ambiguities
############################################################################
############################################################################
cdef class Ambiguity:
    cdef readonly:
        bytes _ABC
        int _Ai, _Ci, _Aj, _Cj, _i, _j
        Sig _sig
        NCMonomial _lm
       
    cdef tuple AC(Ambiguity self)
    
    @staticmethod
    cdef list generate_incls(int i, bytes v, int j , bytes w)
    
    @staticmethod
    cdef list generate_with_tries(prefix_trie, suffix_trie, words, oldlen)
            
    cdef void to_crit_pair(Ambiguity self, SigPoly gi, SigPoly gj)    
           
    
