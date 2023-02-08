# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from sage.matrix.matrix2 cimport Matrix
from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse

from cython_modules.ncmonomial cimport NCMonomial
from cython_modules.ncpoly cimport NCPoly
from cython_modules.sig cimport Sig
from cython_modules.sigpoly cimport SigPoly
from cython_modules.labelled_poly cimport LabelledPoly
from cython_modules.ambiguity cimport Ambiguity


############################################################################
############################################################################
# Parent class for all algorithms
############################################################################
############################################################################
cdef class Algorithm:
    cdef public:
        object syz_automaton, lm_automaton, suffix_trie
        list gens, G, H, labGB, F5_rules, degs_gens, quotient
        Py_ssize_t maxiter, maxdeg, count_interval
        Sig sig_bound
  
    cdef tuple compute_crit_pairs(Algorithm self)
    
    cdef list generate_with_tries(Algorithm self, list words, tuple v)

    cdef tuple c_compute_basis(Algorithm self)
    
    cdef void add_F5_rule(Algorithm self, SigPoly g)
    
    cdef bint F5_criterion(Algorithm self, SigPoly p)
            
    cdef tuple reduction(Algorithm self, SigPoly p)
    
    cdef tuple reconstruct_polynomial(Algorithm self, Matrix_rational_sparse A, list rows, list columns)
                
    cdef tuple symbolic_preprocessing(Algorithm self, SigPoly p)
    
    cdef SigPoly find_reducer(Algorithm self, str t, Sig sigma)
    
    cdef list reconstruct_labelled_basis(Algorithm self)
        
    cdef tuple reconstruct_labelled_poly(Algorithm self, Matrix_rational_sparse A, Matrix_rational_sparse T, list rows, list columns)
    
    cpdef void rewrite_cofactors(Algorithm self)