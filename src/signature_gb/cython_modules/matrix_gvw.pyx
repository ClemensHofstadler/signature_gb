# cython: linetrace=True

from __future__ import absolute_import

from cysignals.signals cimport sig_check, sig_on, sig_off
from sage.all import copy, ZZ, primes


from sage.matrix.matrix2 cimport Matrix
from .linear_algebra cimport *

import ahocorasick
from time import time
from bisect import bisect_right
from collections import defaultdict

import sys

############################################################################
############################################################################
# Matrix GVW Algorithm
############################################################################
############################################################################
cdef class Matrix_GVW():
    """
    """
    def __init__(self, M, maxiter=10, maxdeg=-1, sig_bound=None, verbose=0, block_size=2000):
        
        self._gens = [f.copy() for f in M._gens]
        
        self._G = copy(M._G)
        self._H = copy(M._H)
                
        self._syz_automaton = ahocorasick.Automaton()
        self._lm_automaton = ahocorasick.Automaton()
        self._suffix_trie = ahocorasick.Automaton()
        self._sig_automaton = ahocorasick.Automaton()
           
        self._verbose = verbose
        self._maxiter = maxiter
        self._maxdeg = maxdeg
        if sig_bound in ZZ:
            s = self._gens[0]._sig
            x = chr(0).encode()
            sig_bound -= min(len(f.lm()._mon) for f in self._gens)
            sig_bound = s.lrmul(b'',x*sig_bound)            
        self._sig_bound = sig_bound                
        self._simplify = ahocorasick.Automaton()
        self._block_size = block_size

############################################################################
    def compute_basis(self): return self.c_compute_basis() 
    
    cdef tuple c_compute_basis(Matrix_GVW self):
        cdef list P, pairs, new_syz, new_poly
        cdef Py_ssize_t count
        cdef Sig s, sig_bound
                                 
        sig_bound = self._sig_bound
        maxiter = self._maxiter
        verbose = self._verbose
        oldlen = 0
        count = 1
        G = self._G
        
        pairs = self._gens
                        
        while count <= maxiter and len(pairs) > 0:
        
            P,pairs = self.select_pairs(pairs)
                                                                               
            # apply criteria
            P = self.criteria(P)
            
            if P:  
                # turn crit pairs into actual polynomials
                P = [p.to_poly(G[p._i]) if isinstance(p,Ambiguity) else p for p in P] 
                                                    
                # reduce pairs
                new_poly,new_syz = self.reduction(P)
                                                                                                                                                                                                           
                self.update_poly_data(new_poly)
                self.update_syz_data(new_syz)
                            
                if new_poly:
                    new_pairs = self.compute_crit_pairs(oldlen)
                    
                    # delete pairs with too large signature and by applying criteria
                    if sig_bound: new_pairs = [p for p in new_pairs if p._sig < sig_bound]
                    new_pairs = self.criteria(new_pairs)
                    
                    # append new pairs
                    pairs += new_pairs
                                        
                    oldlen = len(self._G)
                                                                                                                                   
            if verbose > 0 and count % 10 == 0:
                print("Iteration " + str(count) + " finished. G has now " + str(len(self._G)) + " elements.\n")
            count += 1
            
        if verbose > 0 and len(pairs) == 0:
            print("All critical pairs were reduced to 0.")
                
        self.minimize_syz()
                 
        return self._G, self._H
############################################################################
############################################################################
# Compute crit pairs
############################################################################
############################################################################         
    cdef list compute_crit_pairs(Matrix_GVW self, Py_ssize_t oldlen=0):
        cdef list words, crit_pairs, syz
        cdef tuple v
        cdef Ambiguity a
        cdef SigPoly p
                                             
        G = self._G
        maxdeg = self._maxdeg
        sig_bound = self._sig_bound
        verbose = self._verbose
        prefix_trie = self._lm_automaton
        suffix_trie = self._suffix_trie
        
        start = time()
        words = [g.lm()._mon for g in G]
        amb = Ambiguity.generate_with_tries(prefix_trie, suffix_trie, words, oldlen)
        if maxdeg > 0: amb = [a for a in amb if a.degree() <= maxdeg]
        if sig_bound: amb = [a for a in amb if a.degree() <= len(sig_bound)]
                
        if verbose > 1:
            print(str(len(amb)) + " ambiguities in total (computation took %.5f)" % (time()-start))
                                         
        start = time()
        for a in amb: a.to_crit_pair(G[a._i],G[a._j])
                        
        # get regular pairs
        crit_pairs = [a for a in amb if a._sig]
                
        if verbose > 1:
            print(str(len(crit_pairs)) + " critical pairs were generated (computation took %.5f)" % (time()-start))
                                            
        return crit_pairs   
############################################################################
############################################################################
# Reduction & Symbolic Preprocessing
############################################################################
############################################################################
    cdef tuple reduction(Matrix_GVW self, list P, bint  trace=False):
        cdef list rows, cols, new_poly, new_syz
        cdef Matrix_sparse A
        cdef SigPoly f
        
        verbose = self._verbose
        
        start = time()
        rows,cols = self.symbolic_preprocessing(P)         
        cols.sort(reverse=True)        
        rows.sort(key=lambda f: f._sig)
        if verbose > 1:
            print("Symbolic preprocessing finished (computation took %.5f)" % (time()-start))
                
        A = set_up_matrix(rows,cols)
        if verbose > 1:
            print("Matrix dimension " + str(A.dimensions()))
                
        if trace:
            return rows, cols, A
        
        start = time()
        new_poly, new_syz = self.reduce_matrix(A, rows, cols)        
        if verbose > 1:
            print("Reduction finished (computation took %.5f)" % (time()-start))
        
        return new_poly, new_syz
###########################################################################
    cdef tuple reduce_matrix(Matrix_GVW self, Matrix_sparse A, list rows, list cols):
        cdef Matrix_sparse B,C,D
        cdef list new_poly, new_syz, new_poly_A, new_syz_A
                              
        block_size = self._block_size
        nc = A._ncols
        nr = A._nrows
                                  
        new_poly = []
        new_syz = [] 
        
        blocks = []
        prev_sig = rows[0]._sig
        for i,s in enumerate([f._sig for f in rows]):
            if s == prev_sig: continue
            blocks.append(i)
            prev_sig = s
        blocks.append(len(rows))
        
        block_idxs = [blocks[bisect_right(blocks, min(nr-1, n * block_size))] for n in range(1,nr//block_size + 1)]        
        if not block_idxs or block_idxs[-1] != blocks[-1]: block_idxs.append(blocks[-1])
        block_idxs.reverse()   
                
        while block_idxs:
                        
            end = block_idxs.pop()
            block_idxs = [n - end for n in block_idxs]
                                                
            # removes the last A._nrows - end rows of A and puts them into C          
            A,C = split_along_rows(A,list(range(end)))
                                                                  
            # reduce A
            block_echelon(A,blocks)
            new_poly_A, new_syz_A = self.reconstruct_polynomials(A,rows,cols,blocks)
            new_poly += new_poly_A
            new_syz += new_syz_A
            rows = rows[end:]
                                                                        
            if not C._nrows: break
            
            # Faugere-Lachartre
            A,blocks,cols = faugere_lachartre(A,C,blocks,cols)
                    
        return new_poly, new_syz   

###########################################################################
    cdef tuple reconstruct_polynomials(Matrix_GVW self, Matrix_sparse A, list rows, list columns, list blocks):
        cdef Py_ssize_t i, j, k
        cdef list mons, poly, syz
        cdef Sig s
        cdef SigPoly p
        cdef bytes m
        
        lm = self._lm_automaton
        G = self._G
        simplify = self._simplify
                                            
        poly = []
        syz = []
        
        data = get_relevant_matrix_data(A,blocks)
        
        for (i,pos,coeffs) in data:            
            s = rows[i]._sig
            # non syzygy
            if pos:
                mons = [columns[k] for k in pos]
                coeffs = list(coeffs)                
                p = SigPoly(NCPoly(coeffs,mons),s)
                m = p.lm()._mon
                redundant = False
                if lm:
                    divisors = list(lm.iter(m))
                    redundant = any( s >= G[i]._sig.lrmul(m[:k-len(mm)+1],m[k+1:]) for k,(i,mm) in divisors )
                if redundant:
                    simplify.add_word(m, (len(m),p))
                else:
                    poly.append(p)
                    
            # syzygy
            else:
                syz.append(s)
        
                out = []

        simplify.make_automaton()
    
        return poly, syz

###########################################################################
    cdef tuple symbolic_preprocessing(Matrix_GVW self, list P):
        cdef Py_ssize_t i
        cdef set T, done
        cdef SigPoly agb
        cdef NCMonomial m, t
        cdef Sig sigma
                
        sigma = max(p._sig for p in P)
                     
        done = set()
        T = {m for p in P for m in p.mons()}
        R = P
                
        if not len(self._lm_automaton): return R, list(T)
        
        while T:
            t = T.pop()
            done.add(t)
            agb = self.find_reducer(<bytes>t._mon)                                
            if agb and agb._sig < sigma:
                R.append(agb)
                T.update({m for m in agb.mons() if m not in done})
        return R, list(done)
###########################################################################
    cdef SigPoly find_reducer(Matrix_GVW self, bytes t):
        
        cdef Py_ssize_t i,k
        cdef bytes m
        cdef list reducers
        cdef SigPoly p
        
        G = self._G
        simplify = self._simplify
        lm = self._lm_automaton
        
        reducers = []
        
        if len(simplify):
            reducers = [ p.lrmul(t[:k-i+1],t[k+1:]) for k,(i,p) in simplify.iter(t) ]
                        
        reducers += [ G[i].lrmul(t[:k-len(m)+1],t[k+1:]) for k,(i,m) in lm.iter(t) ]
        
        if reducers:         
            return min(reducers, key=lambda p : p._sig)

        return None  
############################################################################
############################################################################
# Elimination criteria
############################################################################
############################################################################    
    cdef list criteria(Matrix_GVW self, list P):
                
        P = self.syzygy_criterion(P)  
        P = [p for p in P if not (self.F5_criterion(p) or self.cover_criterion(p))]             
                  
        return P
############################################################################            
    cdef bint cover_criterion(Matrix_GVW self, object p):
        cdef SigPoly g
        cdef NCMonomial m, lm
        cdef bytes a,b, sigma
        
        sig = self._sig_automaton
        if not sig: return False
        
        G = self._G
        sigma = p._sig.to_bytes()
        m = p.lm()
         
        for k,(i,s) in sig.iter(sigma):
            g = G[i]
            
            a = sigma[:k-s+1]
            b = sigma[k+1:]
                        
            lm = g.lm().copy()
            lm.lrmul(a,b)
            if lm < m: return True
        
        return False
            
###########################################################################
    def syzygy_criterion(self, pairs):
        
        syz = self._syz_automaton
    
        if not len(syz): return pairs
        
        syz.make_automaton()
        out = [p for p in pairs if not next(syz.iter(p._sig.to_bytes()),False)]
                         
        return out   
###########################################################################
    cdef bint F5_criterion(Matrix_GVW self, object p):
        cdef bytes a,b,alpha
        cdef Sig beta, sigma, mu
        cdef SigPoly g,h
                
        lm = self._lm_automaton
        sig = self._sig_automaton
        if not lm: return False
        G = self._G
        
        alpha = p._sig.to_bytes()
                
        for k,(i,s) in sig.iter(alpha):
            g = G[i]   
            beta = g._sig
            lm_g = g.lm()._mon
            
            a = alpha[:k-s+1]
            b = alpha[k+1:]
            
            for l,(j,m) in lm.iter(a):
                h = G[j]
                sigma = beta.lrmul(a[l-len(m)+1:], b'')
                mu = h._sig.lrmul(b'', a[l+1:] + lm_g)
                if sigma > mu: return True
                
            for l,(j,m) in lm.iter(b):
                h = G[j]
                sigma =  beta.lrmul(b'', b[:l+1])
                mu = h._sig.lrmul(lm_g + b[:l-len(m)+1], b'')
                if sigma > mu: return True

        return False           
############################################################################
############################################################################
# Additional stuff
############################################################################
############################################################################ 
    def select_pairs(self,pairs):
    
        # TODO: make pairs a dict { deg : [pairs...]}

        pairs.sort(key = lambda p : p._sig)
        d = len(pairs[0]._sig)
        P = [p for p in pairs if len(p._sig) == d]
        P.sort(key = lambda p : p.lm(), reverse=True)
        P = list(({ p._sig : p for p in P}).values())
        pairs = pairs[len(P):]
        
        return P,pairs 
############################################################################     
    cdef void update_poly_data(Matrix_GVW self, list P):
        cdef SigPoly p
        cdef bytes mon,s
        
        G = self._G
        lm = self._lm_automaton
        sig = self._sig_automaton
        suffix = self._suffix_trie
        
        for i,p in enumerate(P):
            mon = p.lm()._mon
            s = p._sig.to_bytes()
            n = len(G)+i
            
            lm.add_word(mon,(n,mon))
            suffix.add_word(mon[::-1],(n,mon))  
            sig.add_word(s,(n,len(s)))
                    
        lm.make_automaton()
        sig.make_automaton()
        self._G += P
############################################################################     
    cdef void update_syz_data(Matrix_GVW self, list S):
        cdef Sig s
        
        bound = self._sig_bound
        syz = self._syz_automaton
        H = self._H
        
        for s in S:
            if bound and not s < bound: continue
            t = s.to_bytes()
            syz.add_word(t,1)
            H.append(s)
        syz.make_automaton()
############################################################################             
    cdef void minimize_syz(Matrix_GVW self):
        
        syz = self._syz_automaton
        H = self._H
        
        H.sort()
        
        B = [s.to_bytes() for s in H]
        
        to_delete = set()
        for i,s in enumerate(B):
            for j,t in enumerate(B[i+1:]):
                if s in t: to_delete.add(i+j+1)
        
        self._H = [s for i,s in enumerate(H) if i not in to_delete]