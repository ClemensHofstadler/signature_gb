from __future__ import absolute_import

############################################################################
############################################################################
# Labelled Module
############################################################################
############################################################################
cdef class LabelledModule:
    cdef public:
        object _parent
        list _original_gens, _gens, _G, _labGB, _H,
        
    
    cdef void verify_basis(LabelledModule self, list G, list H, args, probabilistic=*)
            

    