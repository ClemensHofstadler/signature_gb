# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cpython cimport *

from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse

from cython_modules.ncmonomial cimport NCMonomial
from cython_modules.ncpoly cimport NCPoly
from cython_modules.sig cimport Sig
from cython_modules.sigpoly cimport SigPoly
from cython_modules.labelled_poly cimport LabelledPoly
from cython_modules.ambiguity cimport Ambiguity
from cython_modules.algorithm cimport Algorithm


############################################################################
############################################################################
# Matrix GVW Algorithm
############################################################################
############################################################################
cdef class Matrix_GVW(Algorithm):
    cdef public:
        object simplify
        Py_ssize_t block_size
            
    cdef tuple c_compute_basis(Matrix_GVW self)
            
    cdef tuple compute_crit_pairs_matrix(Matrix_GVW self, Py_ssize_t oldlen=*)
    
    cdef tuple reduction_matrix(Matrix_GVW self, list P)
    
    cdef list update(Matrix_GVW self, list P)
    
    cdef tuple reconstruct_polynomials(Matrix_GVW self, Matrix_rational_sparse A, list rows, list columns, list blocks)
    
    cdef tuple symbolic_preprocessing_matrix(Matrix_GVW self, list P)
    
    cdef SigPoly find_reducer_matrix(Matrix_GVW self, str t)
    
    cdef void update_poly_data(Matrix_GVW self, list P)
        
    cdef void update_syz_data(Matrix_GVW self, list S)
    
    cdef list criteria(Matrix_GVW self, list P)
               
    cdef bint cover_criterion(Matrix_GVW self, SigPoly p)
    
    cdef list reconstruct_syzygies(Matrix_GVW self)
