from __future__ import absolute_import

from cython_modules.labelled_poly cimport LabelledPoly

from python_modules import global_data
from python_modules.auxiliary import flatten, rebalance

from collections import defaultdict


import ahocorasick
from time import time

############################################################################
############################################################################
# Labelled Module
############################################################################
############################################################################
cdef class LabelledModule:
    cdef public:
        object Parent
        list gens, G, labGB, H
        str _monomial_order, _signature_order
        
    cpdef tuple convert_label(LabelledModule self, LabelledPoly lp)
    
    cpdef list reconstruct_syzygies(LabelledModule self)
