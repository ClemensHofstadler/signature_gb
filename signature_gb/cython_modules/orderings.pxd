# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *
from cython_modules.sig cimport Sig

###########################################################################
############################################################################
# Orderings 
############################################################################
############################################################################
cdef bint deglex(str a, str b)

cdef bint multilex(str a, str b)

cdef bint deg_term_over_position(Sig s, Sig t)

cdef bint deg_position_over_term(Sig s, Sig t)