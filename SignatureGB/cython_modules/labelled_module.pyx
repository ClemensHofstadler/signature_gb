from __future__ import absolute_import

from sage.all import QQ, FreeAlgebra

from python_modules import global_data
from python_modules.auxiliary import flatten

from cython_modules.algorithm cimport Algorithm
from cython_modules.sig cimport Sig
from cython_modules.sigpoly cimport SigPoly
from cython_modules.sig_gb_algorithm cimport SigGB
from cython_modules.gvw cimport GVW
from cython_modules.matrix_gvw cimport Matrix_GVW
from cython_modules.orderings cimport *

############################################################################
############################################################################
# Labelled Module
############################################################################
############################################################################
cdef class LabelledModule:
    """
    """
    def __init__(self, F, *X, signature_order='dpot'):
                        
        self.Parent = F[0].parent()
        
        self.set_monomial_order(X)
        
        self.translate_polynomials(F)

        self.set_signature_order(signature_order,len(F))
              
        self.reset()
        
########################################################################## 
    def translate_polynomials(self,F):
        self.gens = []
        for i,f in enumerate(F):
            ei = ('',i+1,'')
            self.gens.append(SigPoly(f,ei))
            
##########################################################################     
    def reset(self):
                
        self.G = []
        self.labGB = []
        self.H = []
############################################################################
    def __repr__(self):
        
        s = "Labelled module generated by\n"
        for i,f in enumerate(self.gens):
            s += str(f._poly) + "^[e_" + str(i+1) + "]\n"
        s += "Monomial order: "
        if self._monomial_order == 'deglex':
            for v in global_data.vars_: s += v + " < "
            s = s[:-3] + "\n"
        else:
            for block in reversed(global_data.var_blocks_):
                for x in block: s += x + " < "
                s = s [:-3]+ " << "
            s = s[:-3]
            
        s += "Signature order: " + self._signature_order + "\n"   
        return s
##########################################################################         
    def set_monomial_order(self, *X):
        
        s = ""
        X = X[0]
                         
        # deglex order
        if len(X) == 1:
            X = X[0]
            self._monomial_order = 'deglex'
            global_data.cmp_mon_ = deglex
                        
            global_data.vars_ = [str(x) for x in X]
            global_data.vars_dict_ = {x : i for i,x in enumerate(global_data.vars_)}
             
        else:
            self._monomial_order = 'multilex'
            global_data.cmp_mon_ = multilex
            global_data.var_blocks_ = []
            for block in X:
                if not block: continue
                global_data.var_blocks_.append([str(x) for x in block])
                        
            global_data.vars_ = flatten(global_data.var_blocks_)
            global_data.vars_dict_ = {x : i for i,x in enumerate(global_data.vars_)}
            global_data.var_blocks_.reverse()
            for i,block in enumerate(global_data.var_blocks_): global_data.var_blocks_[i] = set(block) 
             
############################################################################
    def set_signature_order(self, order, nr_gens):
        
        self._signature_order = order
        
        global_data.module_basis_ = ["e" + str(i+1) for i in range(nr_gens)]

        if self._signature_order == 'dtop':
            global_data.cmp_sig_ = deg_term_over_position
        elif self._signature_order == 'dpot':
            global_data.cmp_sig_ = deg_position_over_term
        elif self._signature_order == 'F5':
            global_data.weight = [None] + [len(f.lm()._mon) for f in self.gens]
            global_data.cmp_sig_ = F5_order
        else:
            raise NotImplementedError("Signature order not implemented")
############################################################################
    def signature_GB(self, maxiter=10, maxdeg=-1, count_interval = 10, sig_bound=None, algorithm='matrix gvw', block_size=2000, quotient=[]):
        self.reset()
        
        if algorithm == 'standard':
            algo = SigGB(self,maxiter,maxdeg,count_interval,sig_bound)
        elif algorithm == 'gvw':
            algo = GVW(self,maxiter,maxdeg,count_interval,sig_bound,quotient)
        elif algorithm == 'matrix gvw':
            algo = Matrix_GVW(self,maxiter,maxdeg,count_interval,sig_bound,block_size,quotient)
        else:
            raise ValueError, "Algorithm " + str(algorithm) + "not implemented"
            
        self.G, self.H = algo.compute_basis()
                       
        return self.G, self.H
############################################################################
    def reconstruct_labelled_basis(self):
        algo = Algorithm(self)
        self.labGB = algo.reconstruct_labelled_basis()
        return self.labGB
############################################################################
    cpdef list reconstruct_syzygies(LabelledModule self):
        algo = Matrix_GVW(self)
        self.H = algo.reconstruct_syzygies()
        return self.H
############################################################################    
    cpdef tuple convert_label(LabelledModule self, LabelledPoly lp):
        
        vars = global_data.vars_ + global_data.module_basis_
        P = FreeAlgebra(QQ,len(vars), vars)
        p,l = lp.to_normal()
        d = {P(ei) : P(str(self.gens[i]._poly)) for i,ei in enumerate(global_data.module_basis_)}
        return self.Parent(str(p)), self.Parent(str(l.subs(d)))
############################################################################
############################################################################
# Membership test
############################################################################
############################################################################
    def membership_test(self, f):
       algo = Algorithm(self)
       return algo.membership_test(f)