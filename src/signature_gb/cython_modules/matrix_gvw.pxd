# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from sage.matrix.matrix_sparse cimport Matrix_sparse

from .nc_monomial cimport NCMonomial
from .nc_polynomial cimport NCPoly
from .sig cimport Sig
from .sig_polynomial cimport SigPoly
from .ambiguity cimport Ambiguity

############################################################################
############################################################################
# Matrix GVW 
############################################################################
############################################################################
cdef class Matrix_GVW():
    cdef public:
        object _simplify, _syz_automaton, _lm_automaton, _suffix_trie, _sig_automaton
        list _gens, _G, _H
        Py_ssize_t _maxiter, _maxdeg, _block_size, _verbose
        Sig _sig_bound
            
    cdef tuple c_compute_basis(Matrix_GVW self)
            
    cdef list compute_crit_pairs(Matrix_GVW self, Py_ssize_t oldlen=*)
    
    cdef tuple reduction(Matrix_GVW self, list P, bint trace=*)
    
    cdef tuple reduce_matrix(Matrix_GVW self, Matrix_sparse A, list rows, list cols)
        
    cdef tuple reconstruct_polynomials(Matrix_GVW self, Matrix_sparse A, list rows, list columns, list blocks)
    
    cdef tuple symbolic_preprocessing(Matrix_GVW self, list P)
    
    cdef SigPoly find_reducer(Matrix_GVW self, bytes t)
    
    cdef void update_poly_data(Matrix_GVW self, list P)
        
    cdef void update_syz_data(Matrix_GVW self, list S)
    
    cdef list criteria(Matrix_GVW self, list P)
         
    cdef bint F5_criterion(Matrix_GVW self, object p)
               
    cdef bint cover_criterion(Matrix_GVW self, object p)
    
    cdef void minimize_syz(Matrix_GVW self)
