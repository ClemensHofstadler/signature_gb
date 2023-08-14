# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from cython_modules.ncpoly cimport NCPoly
from cython_modules.sig cimport Sig
from cython_modules.sigpoly cimport SigPoly

from python_modules.auxiliary import flatten

###########################################################################
############################################################################
# Ambiguities
############################################################################
############################################################################
cdef class Ambiguity:
    """
    Overlap = fiC - Afj
    Inclusion = fi - AfjC
    """
    cdef readonly:
        str _ABC
        Py_ssize_t _A, _C, _i, _j, _deg
        bint _is_overlap
     
    cdef tuple AC(Ambiguity self)
    
    @staticmethod
    cdef list generate(tuple a, tuple b)
    
    @staticmethod
    cdef list generate_incls(tuple a, tuple b)
    
    cdef SigPoly to_crit_pair(Ambiguity self, SigPoly gi, SigPoly gj)
    
    cdef tuple to_crit_pair_poly(Ambiguity self, NCPoly gi, NCPoly gj)
    
    cdef bint is_redundant(Ambiguity self, str V, Py_ssize_t k)
           
    cdef bint overlap_test(Ambiguity self, str V)
    
    cdef bint incl_test(Ambiguity self, str V, Py_ssize_t k)