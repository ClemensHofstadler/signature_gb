# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from sage.matrix.matrix2 cimport Matrix

from cython_modules.ncmonomial cimport NCMonomial
from cython_modules.ncpoly cimport NCPoly
from cython_modules.sig cimport Sig
from cython_modules.sigpoly cimport SigPoly
from cython_modules.labelled_poly cimport LabelledPoly
from cython_modules.ambiguity cimport Ambiguity
from cython_modules.algorithm cimport Algorithm


############################################################################
############################################################################
# GVW Algorithm
############################################################################
############################################################################
cdef class GVW(Algorithm):
    
    cdef tuple c_compute_basis(GVW self)
                
    cdef bint cover_criterion(GVW self, SigPoly p)
    