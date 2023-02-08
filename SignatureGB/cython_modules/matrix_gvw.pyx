# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cysignals.signals cimport sig_check, sig_on, sig_off

from sage.all import ZZ, QQ, FreeAlgebra, matrix, copy, ceil,N
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
import itertools

import ahocorasick
from time import time

from bisect import bisect_right

    
############################################################################
############################################################################
# Matrix GVW Algorithm
############################################################################
############################################################################
cdef class Matrix_GVW(Algorithm):
    """
    """
    def __init__(self, M, maxiter=10, maxdeg=-1, count_interval=10, sig_bound=None, block_size=2000, quotient=[]):
        
        super().__init__(M,maxiter,maxdeg,count_interval,sig_bound,quotient)
        
        self.simplify = ahocorasick.Automaton()
        self.block_size = block_size

############################################################################
    cdef tuple c_compute_basis(Matrix_GVW self):        
        cdef list G, H, pairs, H_new, new_syz, new_poly
        cdef Py_ssize_t count, i
        cdef Ambiguity amb
        cdef SigPoly g,p
        cdef Sig syz, s, s2, prev_sig, sig_bound
                
        self.update_poly_data([SigPoly(f,('',0,'')) for f in self.quotient])

        pairs = copy(self.gens)
        for i,f in enumerate(pairs):
            f._deg = len(f.lm()._mon)
            self.degs_gens[i+1] = f._deg
                         
        count = 1
        sig_bound = self.sig_bound
        oldlen = len(self.quotient)
                
        while count <= self.maxiter and len(pairs) > 0:
                       
            # sort critical pairs
            pairs.sort(key = lambda p : p._sig._len)
                        
            # select pairs
            d = pairs[0]._sig._len
            P = [p for p in pairs if p._sig._len == d]
            pairs = pairs[len(P):]
                                                                                
            # apply criteria
            P = self.criteria(P)

            if P:
                                        
                # reduce pairs         
                new_poly,new_syz = self.reduction_matrix(P)
                                                                                                                                                                                                                                                      
                self.update_poly_data(new_poly)
                self.update_syz_data(new_syz)
            
                if new_poly:
                    new_pairs, new_syz = self.compute_crit_pairs_matrix(oldlen)
                    
                    self.update_syz_data(new_syz)
                    oldlen = len(self.G)
                                    
                    # delete pairs with too large signature
                    if sig_bound: new_pairs = [p for p in new_pairs if p._sig < sig_bound]
                                                                                                                            
                    # append new pairs
                    pairs += new_pairs
                                                                                                              
            if global_data.VerboseLevel > 0 and count % self.count_interval == 0:
                print("Iteration " + str(count) + " finished. G has now " + str(len(self.G)) + " elements.\n")
            count += 1
            
        if global_data.VerboseLevel > 0 and count < self.maxiter:
            print("All critical pairs were reduced to 0.")
        
        self.H = list(set(self.H))
                 
        return self.G, self.H
############################################################################
############################################################################
# Compute crit pairs
############################################################################
############################################################################         
    cdef tuple compute_crit_pairs_matrix(Matrix_GVW self, Py_ssize_t oldlen=0):
        cdef list words, new_words, amb, crit_pairs
        cdef tuple v
        cdef Py_ssize_t i
        cdef Ambiguity a
        cdef SigPoly g,p
                 
        words = [(i,g.lm()._mon) for i,g in enumerate(self.G)]
        new_words = words[oldlen:]
        
        start = time()
        amb = [self.generate_with_tries(words,v) for v in new_words]        
        amb = flatten(amb)
        if self.maxdeg > 0: amb = [a for a in amb if a._deg <= self.maxdeg]
        
        amb = list(set(amb))
        
        if global_data.VerboseLevel > 0:
            print(str(len(amb)) + " ambiguities in total (computation took %.5f)" % (time()-start))
        
        crit_pairs = [a.to_crit_pair(self.G[a._i],self.G[a._j]) for a in amb]
                
        # get syzygies
        syz = [p._sig for p in crit_pairs if p._deg == -1]
        
        # get regular non-zero S-polynomials
        crit_pairs = [p for p in crit_pairs if p._deg > -1]
        
        if global_data.VerboseLevel > 0:
            print(str(len(crit_pairs)) + " critical pairs were generated.")
                                
        return crit_pairs, syz      
############################################################################
############################################################################
# Reduction & Symbolic Preprocessing
############################################################################
############################################################################
    cdef tuple reduction_matrix(Matrix_GVW self, list P):
        cdef list rows, cols, new_poly, new_poly_A, new_syz, new_syz_A, block_idxs
        cdef Matrix_rational_sparse A,B,C,D
        cdef SigPoly f
        cdef Sig s
        cdef Py_ssize_t end, nr, n
          
        rows,cols = self.symbolic_preprocessing_matrix(P)         
        cols.sort(reverse=True)
        rows.sort(key=lambda f: f._sig)
                                
        blocks = []
        prev_sig = rows[0]._sig
        for i,s in enumerate([f._sig for f in rows]):
            if s == prev_sig: continue
            blocks.append(i)
            prev_sig = s
        blocks.append(len(rows))
                                  
        A = set_up_matrix(rows,cols)
        nr = A._nrows
                            
        new_poly = []
        new_syz = []
            
        block_idxs = [blocks[bisect_right(blocks, min(nr-1, n * self.block_size))] for n in range(1,nr//self.block_size + 1)]        
        if not block_idxs or block_idxs[-1] != blocks[-1]: block_idxs.append(blocks[-1])
        block_idxs.reverse()

        while block_idxs:
            
            end = block_idxs.pop()
            block_idxs = [n - end for n in block_idxs]
                                                
            # removes the last A._nrows - end rows of A and puts them into C          
            A,C = split_along_rows(A,list(range(end)))
                                          
            # reduce A
            block_rational_echelon(A,blocks)
            new_poly_A, new_syz_A = self.reconstruct_polynomials(A,rows,cols,blocks)
            new_poly += new_poly_A
            new_syz += new_syz_A
            rows = rows[end:]
                                                                        
            if not C._nrows: break
            # Faugere
            
            # reorder         
            A,B,C,D,blocks,cols = reorder(A,C,blocks,cols)
                  
            if not is_zero(C): 
                # compute A^{-1}B
                trsm(A,B)              
                #compute difference D - CA^{-1}B        
                C = _matrix_times_matrix_(C,B)               
                diff(D,C) 
                                            
            A = D
                
        new_poly = self.update(new_poly)
                         
        return new_poly, new_syz
###########################################################################
    cdef list update(Matrix_GVW self, list P):
        cdef list out
        cdef SigPoly p
        cdef Sig s
        cdef Py_ssize_t i,k
        cdef str a,b,lm, m
        
        if not len(self.lm_automaton): return P
        
        out = []
        
        for p in P:
            lm = p.lm()._mon
            s = p._sig
            divisors = list(self.lm_automaton.iter(lm))
            if all( s < self.G[i]._sig.lrmul(lm[:k-len(m)+1],lm[k+1:]) for k,(i,m) in divisors ):
                out.append(p)
            else:
                m = p.lm()._mon
                self.simplify.add_word(m, (len(m),p))
        
        self.simplify.make_automaton()                
        return out
###########################################################################
    cdef tuple reconstruct_polynomials(Matrix_GVW self, Matrix_rational_sparse A, list rows, list columns, list blocks):
        cdef Py_ssize_t i, j, k
        cdef list mons, poly, syz
        cdef mpq_vector* v
        cdef bint flag
        cdef Sig s
                     
        poly = []
        syz = []
        
        j = 0
        for i from 0 <= i < A._nrows:
            if i == blocks[j]: j += 1
            flag = True
            v = &(A._matrix[i])
            for k from i < k < blocks[j]:
                if vector_equals(v,&(A._matrix[k])):
                    flag = False
                    break
            if flag:
                s = rows[i]._sig
                # non syzygy
                if v.num_nonzero:
                    pos,coeffs = zip(*mpq_vector_to_list(v))
                    mons = [columns[k] for k in pos]
                    coeffs = list(reversed(list(coeffs)))
                    mons.reverse()
                    poly.append(SigPoly(NCPoly(mons,coeffs),s))
                # syzygy
                else:
                    syz.append(s)
    
        return poly, syz
###########################################################################
    cdef tuple symbolic_preprocessing_matrix(Matrix_GVW self, list P):
        cdef Py_ssize_t i
        cdef set T, R, done
        cdef SigPoly agb
        cdef NCMonomial m, t
        cdef Sig sigma
        
        sigma = max(P,key=lambda p : p._sig)._sig
                     
        done = set()
        T = {m for p in P for m in p.mons()}
        R = set(P)
                
        if not len(self.lm_automaton): return list(R), list(T)
        
        while len(T):
            t = T.pop()
            done.add(t)
            agb = self.find_reducer_matrix(t._mon)                                      
            if agb and agb._sig < sigma:
                R.add(agb)
                T.update({m for m in agb.mons() if m not in done})
        return list(R), list(done)
###########################################################################
    cdef SigPoly find_reducer_matrix(Matrix_GVW self, str t):
        
        cdef Py_ssize_t i,k
        cdef str m
        cdef list reducers
        cdef SigPoly p
        
        reducers = []
        
        if len(self.simplify):
            reducers = [ p.lrmul(t[:k-i+1],t[k+1:]) for k,(i,p) in self.simplify.iter(t) ]
                        
        reducers += [ self.G[i].lrmul(t[:k-len(m)+1],t[k+1:]) for k,(i,m) in self.lm_automaton.iter(t) ]
        
        if reducers: 
            return min(reducers,key=lambda p : p._sig)

        return None                  
############################################################################
############################################################################
# Additional stuff
############################################################################
############################################################################ 
    cdef void update_poly_data(Matrix_GVW self, list P):
        cdef SigPoly p
        cdef str mon
        
        for i,p in enumerate(P):
            mon = p.lm()._mon  
            self.lm_automaton.add_word(mon,(len(self.G)+i,mon))
            self.suffix_trie.add_word(mon[::-1],(len(self.G)+i,mon))   
            self.add_F5_rule(p) 
        self.lm_automaton.make_automaton()
        self.G += P
############################################################################     
    cdef void update_syz_data(Matrix_GVW self, list S):
        cdef Sig s
        for s in S:
            if self.sig_bound and not s < self.sig_bound: continue
            self.syz_automaton.add_word(s.my_str(),1)
            self.H.append(s)
        self.syz_automaton.make_automaton()
############################################################################ 
    cdef list criteria(Matrix_GVW self, list P):
        cdef SigPoly p
        
        P = self.syzygy_criterion(P)        
        P = [p for p in P if not (self.F5_criterion(p) or self.cover_criterion(p))]
        
        return P
############################################################################            
    cdef bint cover_criterion(Matrix_GVW self, SigPoly p):
        cdef Sig sigma
        cdef SigPoly g
        cdef tuple ab
        cdef NCMonomial m, lm
        cdef str a,b
                     
        sigma = p._sig
        m = p.lm()
        for g in self.G:
            ab = g._sig.divides(sigma)
            if ab:
                a,b = ab
                lm = g.lm().copy()
                lm.lrmul(a,b)
                if lm < m: return True
        return False
        
###########################################################################
    cdef list reconstruct_syzygies(Matrix_GVW self):
        cdef list pivot_rows, rows, cols, syz, coeffs, mons, sigGB, F
        cdef Matrix_rational_sparse A, T, B
        cdef SigPoly g
        cdef Sig s, m
        cdef str a, b
        cdef Py_ssize_t i, j, k 
        cdef LabelledPoly r
        cdef Rational c, c2
        cdef mpq_t prod, v
        
        mpq_init(prod)
        mpq_init(v)
     
        d = defaultdict(Rational)
        
        sigGB = self.G
        self.G = self.labGB
        
        self.lm_automaton.clear()
        for j,g in enumerate(self.G):
            mon = g.lm()._mon
            self.lm_automaton.add_word(mon,(j,mon))
        self.lm_automaton.make_automaton()
        
        syz = self.H
        self.H = []
        F = [None] + copy(self.gens)    
        P = [F[s._ei].lrmul(s._a,s._b) for s in syz]
        
            
        rows,cols = self.symbolic_preprocessing_matrix(P)            
        set_P = set(P)
        pivot_rows = [f for f in rows if f not in set_P]        
        cols.sort(reverse=True)
        pivot_rows.sort(key=lambda g: g._sig)
                                        
        blocks = [len(pivot_rows),len(rows)]                                          
        A = set_up_matrix(pivot_rows + P,cols)
        A = augment(A)
        block_rational_echelon(A,blocks)
        T = A.matrix_from_rows_and_columns(range(len(pivot_rows),A._nrows),range(A._ncols-A._nrows,A._ncols))
        
        cdef mpq_vector* w
        nonzero_column = [False for _ in range(len(pivot_rows))]
        for i from 0 <= i < T._nrows:
            w = &(T._matrix[i])
            for j from 0 <= j < w.num_nonzero-1:
                nonzero_column[w.positions[j]] = True
                
        for j from 0 <= j < len(pivot_rows):
            if not nonzero_column[j]: continue
            r = pivot_rows[j]
            a,_,b = r._pseudo_sig.aib()
            r._module_mons = [m.lrmul(a,b) for m in r._module_mons]
        
        for i from 0 <= i < T._nrows:
            d.clear()
            pos,values = zip(*mpq_vector_to_list(&T._matrix[i]))
            d[syz[i]] = values[-1]
            for j from 0 <= j < len(pos)-1:
                r = pivot_rows[pos[j]]
                mpq_set(v, (<Rational>values[j]).value)
                for k from 0 <= k < len(r._module_coeffs):
                    c = <Rational>r._module_coeffs[k]
                    m = <Sig>r._module_mons[k]
                    c2 = d[m]
                    sig_on()
                    mpq_mul(prod, c.value, v)
                    mpq_add(c2.value, prod, c2.value)
                    sig_off()
                    d[m] = c2
            
            coeffs = [c for c in d.values() if c]
            mons = [m for m in d if d[m]]
   
            r = LabelledPoly(NCPoly.zero(), syz[i], syz[i], mons, coeffs)
            self.H.append(r)
        
        mpq_clear(prod)
        mpq_clear(v)
                                   
        self.G = sigGB  
        return self.H