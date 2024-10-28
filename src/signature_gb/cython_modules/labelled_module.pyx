from __future__ import absolute_import

from sage.all import flatten, GF, FreeAlgebra, rational_reconstruction, ZZ, QQ, copy

import os

from .free_algebra cimport MyFreeAlgebra
from .sig cimport Sig
from .nc_polynomial cimport NCPoly
from .sig_polynomial cimport SigPoly
from .matrix_gvw cimport Matrix_GVW
from .free_bimodule cimport FreeBimodule
import signature_gb.auxiliary as aux

from sage.arith.misc import random_prime, primes
from sage.arith.multi_modular cimport MultiModularBasis
import sage.parallel.multiprocessing_sage as MP
from sage.parallel.decorate import normalize_input

from collections import defaultdict, Counter

from sage.modules.vector_modn_sparse cimport *
from signature_gb.cython_modules.linear_algebra cimport *
from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse


import random
import time

MAX_MOD = 46341 #due to Sage internal stuff

############################################################################
############################################################################
# Labelled Module
############################################################################
############################################################################
cdef class LabelledModule:
    """
    """
    def __init__(self, gens, monomial_order=None, signature_order='DoPoT'):
        
        if signature_order not in {'DoPoT', 'DoToP'}:
            raise NotImplementedError("Signature order not implemented")
            
        F = gens[0].parent()
        R = F.base_ring()
        X = monomial_order if monomial_order else F.gens()
        
        self._original_gens = gens
        F = MyFreeAlgebra(R,X)
        self._parent = FreeBimodule(F,len(gens),signature_order=signature_order)
        
        self.translate_polynomials(gens)                      
        self.reset()
############################################################################
    def __repr__(self):
        
        s = "Labelled module generated by\n"
        for i,f in enumerate(self._gens):
            s += str(f._poly) + "^[e_" + str(i+1) + "]\n"
        s += "over " + str(self._parent)
        return s
##########################################################################     
    def reset(self):      
        self._G = []
        self._labGB = []
        self._H = []
##########################################################################         
    def translate_polynomials(self,gens):
        self._gens = []
        B = self._parent
        F = self._parent.F()
                  
        for i,f in enumerate(gens):
            g = F(f)
            ei = Sig(b'',i+1,b'',B)
            ei.set_len(g.degree())
            s = SigPoly(g, ei)
            self._gens.append(s)                   
############################################################################
    def characteristic(self):
        return self._parent.base_ring().characteristic()
############################################################################
    def signature_basis(self, maxiter=10, maxdeg=-1, sig_bound=None, verbose=0):
        
        self.reset()
        algo = Matrix_GVW(self,maxiter,maxdeg,sig_bound,verbose)
        self._G, self._H = algo.compute_basis()
                       
        return self._G, self._H
############################################################################
    def modular_signature_basis(self, maxiter=10, maxdeg=-1, sig_bound=None,\
        verbose=0, threads=None, num_primes=None, verification='rigorous'):
            
        if not threads:
            threads = os.cpu_count()
        if not num_primes:
            num_primes = max(8,threads)
            
        original_gens = self._original_gens
        X = self._parent.gens()
        order = self._parent.signature_order()
            
        if num_primes > 2892:
            raise ValueError("Not enough primes. At most 2892 primes are supported.")
        P = []; P_new = []
        while len(P_new) < num_primes: 
            p = random_prime(2**14, MAX_MOD)
            if p not in P_new: P_new.append(p)
            
        Gs = []; Hs = []
        while True:
            P += P_new
            # prepare input for different primes       
            args = []
            for p in P_new:
                R = FreeAlgebra(GF(p), X)
                gens = [R(f) for f in original_gens]
                M = LabelledModule(gens, signature_order=order)
                arg = normalize_input((M,maxiter,maxdeg,sig_bound,verbose))
                args.append(arg)
        
            # compute sig-GBs in parallel 
            out = [(i[0][0],o) for (i,o) in MP.parallel_iter(threads, aux.compute_basis, args)]   
            Gs += [(p,GH[0]) for p,GH in out]
            Hs += [tuple(GH[1]) for p,GH in out]
        
            # reconstruct bases from images  
            try:                
                G = self.reconstruct_basis_from_primes(Gs)
                H = list(Counter(Hs).most_common(1)[0][0])
                break
            except ArithmeticError:
                P_new = []
                while len(P_new) < threads: 
                    p = random_prime(2**14, MAX_MOD)
                    if p not in P_new: P_new.append(p)
                if verbose > 0:
                    print("Basis could not be reconstructed with %d primes " % len(P))
                    print("Increasing to %d" % (len(P) + threads))
                
        # final verification 
        args = (maxdeg,sig_bound,verbose) 
        if verification == 'probabilistic':
            self.verify_basis(G,H,args,probabilistic=True)
        elif verification == 'rigorous':
            self.verify_basis(G,H,args,probabilistic=False)
            
        return G,H
############################################################################
    def reconstruct_basis_from_primes(self, images):
        
        primes = [im[0].characteristic() for im in images]
        imgs = [im[1] for im in images]
             
        sigs = defaultdict(lambda : [None] * len(primes))
             
        for i,Gp in enumerate(imgs):
            for g in Gp:
                s = g._sig
                sigs[s][i] = g._poly
                
        B = MultiModularBasis(primes)
        G = [SigPoly.reconstruct_from_images(s,g,B) for s,g in sigs.items()]
        
        # change base ring of monomials
        P = self._parent.F() 
        for g in G:
            for m in g._poly._mons: m.change_parent(P)
        
        return G
############################################################################    
    cdef void verify_basis(LabelledModule self, list G, list H, args, probabilistic=False):
        
        verbose = args[-1]
        if verbose > 0:
            print("Verifying reconstructed basis...")
        
        # verify that signature data is correct
        start = time.time()
        self.verify_signature_data(G,H,probabilistic=probabilistic)
        if verbose > 0:
            print("Verifying signature data took %.3f" % (time.time() - start))
            
        # check cover criterion
        start = time.time()
        algo = Matrix_GVW(self, maxdeg=args[0], sig_bound=args[1])
        algo.update_poly_data(G)
        algo.update_syz_data(H)
                      
        # create all critical pairs
        P = algo.compute_crit_pairs(oldlen=0)
        
        # remove pairs with too large signature
        s = algo._sig_bound
        if s: P = [p for p in P if p._sig < s]
        
        # check cover criterion
        P = algo.criteria(P)
        if P:
            raise ValueError("Cover criterion failed for signature basis")
        if verbose > 0:
            print("Verifying cover criterion took %.3f" % (time.time() - start))
     
############################################################################
    def verify_signature_data(self, G, H, probabilistic=False):
                    
        F = self._gens

        algo = Matrix_GVW(self)
        G_verified = []
        
        if probabilistic:
            p = random_prime(MAX_MOD, lbound=10**4)
            G = [g.mod(p) for g in G]
        
        D = max(max(len(g._sig) for g in G), max((len(s) for s in H),default=0))
        
        for d in range(D):
            Gd = [g for g in G if len(g._sig) == d]
            Hd = [s for s in H if len(s) == d]
            if not (Gd or Hd): continue
            
            R = []
            for g in Gd + Hd:
                flag = False
                s = g._sig if isinstance(g,SigPoly) else g
                for gg in reversed(G_verified):
                    ab = gg._sig.divides(s)
                    if ab:
                        x,y = ab
                        h = gg.lrmul(x,y)
                        flag = True
                        break
                if not flag: 
                    h = s.to_poly(F)
                    if probabilistic: h = h.mod(p)

                if isinstance(g,SigPoly) and h._poly == g._poly: continue
                R.append(h)
            
            rows, cols, M = algo.reduction(Gd + R,trace=True)
            row_idxs = [i for i,r in enumerate(rows) if r in R]
            
            H,M = split_along_rows(M,row_idxs)
            blocks = list(range(M.nrows()))
            D,_,_ = faugere_lachartre(M,H,blocks,cols)

            if not D.is_zero():
                raise ValueError("Polynomial-signature relation failed for signature basis")
            
            G_verified += Gd
            algo.update_poly_data(Gd)
                                 
############################################################################
