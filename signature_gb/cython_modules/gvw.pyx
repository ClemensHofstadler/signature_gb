# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cysignals.signals cimport sig_check, sig_on, sig_off

from sage.all import ZZ, QQ, FreeAlgebra, matrix, copy
from sage.modules.vector_modn_sparse cimport *
from sage.modules.vector_rational_sparse cimport *
from sage.libs.gmp.mpq cimport *
from sage.matrix.matrix2 cimport Matrix
from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse
from sage.rings.rational cimport Rational

from cython_modules.linear_algebra cimport *
from cython_modules.orderings cimport *

from python_modules import global_data
from python_modules.auxiliary import flatten, rebalance

from collections import defaultdict

import ahocorasick
from time import time

############################################################################
############################################################################
# GVW Algorithm
############################################################################
############################################################################
cdef class GVW:
    """
    """
    def __init__(self, M, maxiter=10, maxdeg=-1, count_interval=10, sig_bound=None, quotient=[]):
        
        super().__init__(M,maxiter,maxdeg,count_interval,sig_bound,quotient)

############################################################################
    cdef tuple c_compute_basis(GVW self):        
        cdef list G, H, pairs, H_new, new_syzs
        cdef Py_ssize_t count, i
        cdef Ambiguity amb
        cdef SigPoly new_poly,g,p
        cdef Sig new_syz, syz, s, s2, prev_sig, sig_bound
                
        pairs = copy(self.gens)
        for i,f in enumerate(pairs):
            f._deg = len(f.lm()._mon)
            self.degs_gens[i+1] = f._deg
        
        pairs.sort(key = lambda p : p._sig, reverse=True)
        pairs = [pairs]
                 
        count = 1
        reductions = 0
        zero_reductions = 0
        
        sig_bound = self.sig_bound
        
        while count <= self.maxiter and len(pairs) > 0:
           
            # sort critical pairs
            pairs.sort(key = lambda l : l[-1]._sig, reverse=True)
                                                         
            p = pairs[-1].pop()
            if not pairs[-1]: pairs.pop()
                                                   
            if self.F5_criterion(p): continue
            if self.cover_criterion(p): continue
            
            reductions += 1
                         
            new_poly,new_syz = self.reduction(p)
                                                                                                           
            if new_syz:
                zero_reductions += 1
                new_syz_str = new_syz.my_str()
                self.syz_automaton.add_word(new_syz_str,1)
                self.H.append(new_syz)
                pairs = [[p for p in l if not new_syz_str in p._sig.my_str()] for l in pairs]
                pairs = [l for l in pairs if l]
                continue
                                                         
            self.update_poly_data(new_poly)                                  
            new_pairs, new_syzs = self.compute_crit_pairs()
            
            # don't apply new syzygies to new pairs - not worth it   
            for new_syz in new_syzs:
                if new_syz < sig_bound:
                    self.syz_automaton.add_word(new_syz.my_str(),1)
                    self.H.append(new_syz)
                
            # delete pairs with too large signature
            if sig_bound: new_pairs = [p for p in new_pairs if p._sig < sig_bound]
                                                                                        
            # apply old syzygies to new pairs
            new_pairs = self.syzygy_criterion(new_pairs) 
            new_pairs.sort(key = lambda p: p._sig, reverse=True) 
                                    
            # append new pairs
            if new_pairs: pairs.append(new_pairs)
                                                          
            # rebalance pairs if necessary
            if len(pairs) > 2: pairs = rebalance(pairs)   
                                                           
            if global_data.VerboseLevel > 0 and count % self.count_interval == 0:
                print("Iteration " + str(count) + " finished. G has now " + str(len(self.G)) + " elements.\n")
            count += 1
            
        if global_data.VerboseLevel > 0 and count < self.maxiter:
            print("All critical pairs were reduced to 0.")
        
        print("Reductions = %d" % reductions)
        print("Reductions to 0 = %d" % zero_reductions)
                 
        return self.G, self.H
############################################################################
############################################################################
# Additional stuff
############################################################################
############################################################################        
    cdef bint cover_criterion(GVW self, SigPoly p):
        cdef Sig sigma
        cdef SigPoly g
        cdef tuple ab
        cdef NCMonomial m
        cdef str a,b
                     
        sigma = p._sig
        m = p.lm()
        for g in self.G:
            ab = g._sig.divides(sigma)
            if ab:
                a,b = ab
                if g._poly.lrmul(a,b).lm() < m: return True
        return False