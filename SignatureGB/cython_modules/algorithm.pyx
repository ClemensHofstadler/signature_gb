# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from cysignals.signals cimport sig_check, sig_on, sig_off

from sage.all import ZZ, QQ, FreeAlgebra, matrix, copy, vector,N
from sage.modules.vector_modn_sparse cimport *
from sage.modules.vector_rational_sparse cimport *
from sage.libs.gmp.mpq cimport *
from sage.rings.rational cimport Rational
from sage.structure.element cimport Vector


from cython_modules.linear_algebra cimport *
from cython_modules.orderings cimport *

from python_modules import global_data
from python_modules.auxiliary import flatten, rebalance

from collections import defaultdict

import ahocorasick
from time import time

############################################################################
############################################################################
# Algorithm
############################################################################
############################################################################
cdef class Algorithm:
    """
    """
    def __init__(self, M, maxiter=10, maxdeg=-1, count_interval=10, sig_bound=None, quotient=[]):
        
        self.gens = M.gens
        self.quotient = quotient
        
        self.G = M.G
        self.H = M.H
        self.labGB = M.labGB
                
        self.syz_automaton = ahocorasick.Automaton()
        self.lm_automaton = ahocorasick.Automaton()
        self.suffix_trie = ahocorasick.Automaton()
                
        self.F5_rules = [[] for i in range(len(self.gens)+1)]
        self.degs_gens = [0 for i in range(len(self.gens)+1)]
        
        self.maxiter = maxiter
        self.maxdeg = maxdeg
        self.count_interval = count_interval
        self.sig_bound = sig_bound
############################################################################
    def compute_basis(self): return self.c_compute_basis() 
############################################################################
    cdef tuple c_compute_basis(Algorithm self): pass
############################################################################
    cdef tuple compute_crit_pairs(Algorithm self):
        cdef list words, amb, crit_pairs, syz
        cdef tuple v, w
        cdef Ambiguity a
        cdef SigPoly p

        words = list(self.lm_automaton.values())
        v = (len(self.G)-1, self.G[-1].lm()._mon)
        
        # naive way
        # amb = flatten([Ambiguity.generate(v,w) for w in words])
                
        amb = self.generate_with_tries(words, v) 
        if self.maxdeg > 0: amb = [a for a in amb if a._deg <= self.maxdeg]
                          
        amb = list(set(amb))
                     
        crit_pairs = [a.to_crit_pair(self.G[a._i],self.G[a._j]) for a in amb]
                
        # get syzygies
        syz = [p._sig for p in crit_pairs if p._deg == -1]
        
        # get regular non-zero S-polynomials
        crit_pairs = [p for p in crit_pairs if p._deg > -1]
                
        return crit_pairs, syz        
############################################################################        
    cdef list generate_with_tries(Algorithm self, list words, tuple v):
        cdef list amb
        cdef Py_ssize_t i,j,k,len_m
        cdef str m,m_rev,A,B,C,AB,BC
        cdef tuple w
          
        prefix_trie = self.lm_automaton
        suffix_trie = self.suffix_trie
        
        amb = []
        i,m = v
        m_rev = m[::-1]
        len_m = len(m)

        for k from 1 <= k < len_m:
            # overlap ambiguities with m = AB
            A = m[:-k]
            B = m[-k:]
            amb += [Ambiguity(A+BC,len(A),len_m,i,j,True) for j,BC in prefix_trie.values(B) if len(BC) > k]

            # overlap ambiguities with m = BC
            B = m_rev[-k:]
            C = m[k:]
            amb += [Ambiguity(AB+C,len(AB)-k,len(AB),j,i,True) for j,AB in suffix_trie.values(B) if len(AB) > k]
   
        # inclusion ambiguities with m = ABC        
        amb += [Ambiguity(m,k-len(B)+1,k+1,i,j,False) for k,(j,B) in prefix_trie.iter(m) if i != j]
        # inclusion ambiguities with m = B
        amb += flatten([Ambiguity.generate_incls(w,v) for w in words if len(w[1]) > len_m])  
        
        return amb

############################################################################
############################################################################
# Elimination criteria
############################################################################
############################################################################    
    cdef void add_F5_rule(Algorithm self, SigPoly g):    
        cdef str a,b, mon
        cdef Py_ssize_t i, j, lg, lab, lf
        
        a,j,b = g._sig.aib()
        mon = g.lm()._mon
        lg = len(mon)
        lab = len(a) + len(b)
               
        for i,lf in enumerate(self.degs_gens):
            if not i: continue
            if lg > lab + lf: self.F5_rules[i].append(mon)
            if lg == lab + lf and j < i: self.F5_rules[i].append(mon) 
###########################################################################
    cdef bint F5_criterion(Algorithm self, SigPoly p):
        cdef str a,b,m
        cdef Py_ssize_t i 
        
        a,i,b = p._sig.aib()
        for m in self.F5_rules[i]:
            if m in a or m in b: return True
        return False
###########################################################################
    def syzygy_criterion(self, pairs):
    
        if not len(self.syz_automaton): return pairs
        
        self.syz_automaton.make_automaton()
        out = [p for p in pairs if not next(self.syz_automaton.iter(p._sig.my_str()),False)]
                 
        return out
############################################################################
############################################################################
# Reduction & Symbolic Preprocessing
############################################################################
############################################################################
    cdef tuple reduction(Algorithm self, SigPoly p):
        cdef list rows, cols
        cdef Matrix_rational_sparse A
        cdef SigPoly f
                
        rows,cols = self.symbolic_preprocessing(p)
    
        cols.sort(reverse=True)
        rows.sort(key=lambda f: f._sig)
            
        A = set_up_matrix(rows,cols)                                   
        rational_echelon(A) 
           
        new_poly,new_syz = self.reconstruct_polynomial(A,rows,cols)
            
        return new_poly,new_syz
###########################################################################
    cdef tuple reconstruct_polynomial(Algorithm self, Matrix_rational_sparse A, list rows, list columns):
        cdef Py_ssize_t i, j
        cdef list pos, coeffs, mons
        cdef SigPoly new_poly
        cdef Sig s
        
        i = A.nrows()-1
        s = rows[i]._sig
        
        pos = mpq_vector_to_list(&A._matrix[i])
        
        if not pos: return None, s
        
        p,c = zip(*pos)
        mons = [columns[j] for j in p]
        coeffs = list(c)
        mons.reverse()
        coeffs.reverse()
        
        new_poly = SigPoly(NCPoly(mons,coeffs),s)

        return new_poly, None
###########################################################################
    cdef tuple symbolic_preprocessing(Algorithm self, SigPoly p):
        cdef Py_ssize_t i
        cdef set T, R, done
        cdef Sig sigma
        cdef SigPoly agb
        cdef NCMonomial m, t
             
        sigma = p._sig
        done = set()
        T = set(p.mons())
        R = {p}
        
        if not len(self.lm_automaton): return list(R), list(T)
        
        while len(T) > 0:
            t = T.pop()
            done.add(t)
            agb = self.find_reducer(t._mon,sigma)                                      
            if agb:
                R.add(agb)
                T.update({m for m in agb.mons() if m not in done})
        return list(R), list(done)

###########################################################################
    cdef SigPoly find_reducer(Algorithm self, str t, Sig sigma):
        
        cdef Py_ssize_t i,k
        cdef str a,b,m
        cdef Sig min_sig
        cdef SigPoly g, agb, result
                     
        result = None
        min_sig = sigma
                                  
        for k,(i,m) in self.lm_automaton.iter(t):
            g = self.G[i]
            a = t[:k-len(m)+1]
            b = t[k+1:]
                    
            agb = g.lrmul(a,b)
                        
            if min_sig == None or agb._sig < min_sig:
                result = agb
                min_sig = result._sig             
        return result
############################################################################
    cdef list reconstruct_labelled_basis(Algorithm self):
        cdef list rows, cols, sigs, coeffs, mons, F, Algorithm
        cdef Matrix_rational_sparse A, T, A2, T2
        cdef SigPoly g
        cdef Sig s
        cdef str a, b
        cdef Py_ssize_t i, j
        cdef NCPoly poly
        cdef LabelledPoly lp
        
        sigGB = self.G
        self.G = []
        self.lm_automaton.clear()
        F = [None] + copy(self.gens)
              
        sigs = [g._sig for g in sigGB]
        sigs.sort() 
        
        for j,s in enumerate(sigs):
            a,i,b = s.aib()
            g = F[i].lrmul(a,b)
            lp = LabelledPoly(g._poly, s, s, [], [])
            
            rows,cols = self.symbolic_preprocessing(lp)           
            cols.sort(reverse=True)
            rows.sort(key = lp.__eq__)         
                
            A = set_up_matrix(rows,cols) 
            A = augment(A)     
            A = block_rational_echelon(A,[A._nrows-1,A._nrows])
               
            T = A.matrix_from_rows_and_columns(range(A._nrows),range(A._ncols - A._nrows,A._ncols))
            A = A.matrix_from_rows_and_columns(range(A._nrows),range(A._ncols - A._nrows))
                    
            poly, mons, coeffs = self.reconstruct_labelled_poly(A,T,rows,cols)
            lp = LabelledPoly(poly, s, Sig('',j+1,''), mons, coeffs)
            
            self.G.append(lp)
            mon = lp.lm()._mon
    
            self.lm_automaton.add_word(mon,(j,mon))
            self.lm_automaton.make_automaton()
        
        self.rewrite_cofactors()
        
        self.labGB = self.G
        self.G = sigGB
                   
        return self.labGB
###########################################################################
    cdef tuple reconstruct_labelled_poly(Algorithm self, Matrix_rational_sparse A, Matrix_rational_sparse T, list rows, list columns):
        cdef Py_ssize_t i, j
        cdef list pos, coeffs, mons
        cdef NCPoly poly
        
        i = A._nrows-1     
        p = mpq_vector_to_list(&A._matrix[i])
        if p: 
            p,c = zip(*p)
            coeffs = list(c)
            mons = [columns[j] for j in p]
            coeffs.reverse()
            mons.reverse()
            poly = NCPoly(mons,coeffs)
        else: 
            poly = NCPoly.zero()
        
        p,c = zip(*mpq_vector_to_list(&T._matrix[i]))
        coeffs = list(c)
        mons = [rows[j]._pseudo_sig for j in p]
        
        return poly, mons, coeffs
############################################################################
    cpdef void rewrite_cofactors(Algorithm self):
        cdef LabelledPoly g
        cdef Sig s, s2
        cdef str a,b
        cdef Py_ssize_t i, j
        cdef list mons, coeffs, new_mons, new_coeffs, multiples_coeffs, multiples_mons
        cdef list needed_multiples_coeffs, needed_multiples_mons
        cdef Rational c, c2, c3
        cdef Vector v
        
        d = defaultdict(Rational)
        
        # analyze which multiples are needed
        needed_multiples_coeffs = [set() for i in range(len(self.G))]
        needed_multiples_mons = [set() for i in range(len(self.G))]
        for j,g in enumerate(self.G):
            for c,s in zip(g._module_coeffs[:-1],g._module_mons[:-1]):
                a,i,b = s.aib()
                if i-1 == j: continue
                needed_multiples_coeffs[i-1].add(abs(c))
                needed_multiples_mons[i-1].add((a,b))

            
        multiples_coeffs = [None for i in range(len(self.G))] 
        multiples_mons = [None for i in range(len(self.G))] 
                                                 
        for j,g in enumerate(self.G):
            d.clear()
            mons = g._module_mons
            coeffs = g._module_coeffs
            for c,s in zip(coeffs[:-1],mons[:-1]):
                a,i,b = s.aib()
                if i-1 == j:
                    c2 = d[s]
                    mpq_add(c2.value, c.value, c2.value)
                else:
                    new_coeffs = multiples_coeffs[i-1][abs(c)]
                    if c < 0: new_coeffs = [-c2 for c2 in new_coeffs]
                    new_mons = multiples_mons[i-1][(a,b)]
                                                            
                    for c2,s in zip(new_coeffs, new_mons):
                        c3 = d[s]
                        mpq_add(c3.value, c2.value, c3.value)
                                  
            c = coeffs[-1]
            c2 = d[mons[-1]]
            mpq_add(c2.value, c.value, c2.value)
                 
            g._module_mons = [s for s in d if d[s]]
            g._module_coeffs = [c for c in d.values() if c]
            
            v = vector(g._module_coeffs)
            multiples_coeffs[j] = {c : list(c * v) if c != 1 else list(v) for c in needed_multiples_coeffs[j]} 
            multiples_mons[j] = {(a,b) : [s.lrmul(a,b) for s in g._module_mons] for a,b in needed_multiples_mons[j]}
        
############################################################################
############################################################################
# Membership & correctness checks
############################################################################
############################################################################
    def membership_test(Algorithm self, f):
        cdef list rows, cols, coeffs, mons, Algorithm, new_coeffs, new_mons
        cdef Matrix_rational_sparse A, T
        cdef SigPoly g
        cdef Sig m, m2
        cdef str a, b
        cdef Py_ssize_t i
        cdef NCPoly p, poly
        cdef LabelledPoly r
        cdef Rational c, c2
        
        G = self.G
        self.G = self.labGB
        
        self.lm_automaton.clear()
        for j,g in enumerate(self.G):
            mon = g.lm()._mon
            self.lm_automaton.add_word(mon,(j,mon))
            self.lm_automaton.make_automaton()
     
        p = NCPoly(f)
        lp = LabelledPoly(p, None, None, [], [])
            
        rows,cols = self.symbolic_preprocessing(lp)
                        
        cols.sort(reverse=True)
        rows.remove(lp)         
        rows.sort(key=lambda g: g._sig)
        rows.append(lp)
                
        A = set_up_matrix(rows,cols)         
        nc = A._ncols
                       
        A = rational_echelon(A,transformation=True)
                            
        T = A[:,nc:]
        A = A[:,:nc]
                    
        poly, mons, coeffs = self.reconstruct_labelled_poly(A,T,rows,cols)
        
        if poly != NCPoly.zero():
            print("Membership test FAILED.")
            print("Polynomial has nonzero remainder %s." % str(poly))
            return False
        
        c = - coeffs[-1]
        coeffs = [c2 / c for c2 in coeffs[:-1]]
        mons = mons[:-1]
        
        # rewrite mons in terms of generators
        d = defaultdict(Rational)
        for c,m in zip(coeffs,mons):
            a,i,b = m.aib()
            for c2,m2 in zip(self.G[i-1]._module_coeffs,self.G[i-1]._module_mons):
                d[m2.lrmul(a,b)] += c * c2
             
        mons = [m for m in d if d[m]]
        coeffs = [c for c in d.values() if c]
        lp = LabelledPoly(p, max(mons), None, mons, coeffs)
        
        self.G = G
        
        print("Membership test SUCCESSFUL.")
        return lp
############################################################################
############################################################################
# Additional stuff
############################################################################
############################################################################   
    def update_poly_data(self,new_poly):
        self.G.append(new_poly)
        
        mon = new_poly.lm()._mon  
        self.lm_automaton.add_word(mon,(len(self.G)-1,mon))
        self.lm_automaton.make_automaton()
        self.suffix_trie.add_word(mon[::-1],(len(self.G)-1,mon))
        
        self.add_F5_rule(new_poly)
