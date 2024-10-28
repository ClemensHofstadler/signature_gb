# cython: linetrace=True
# cython: boundscheck=False

from __future__ import absolute_import

from sage.all import flatten

from .nc_monomial cimport NCMonomial

############################################################################
############################################################################
# Ambiguities
############################################################################
############################################################################
cdef class Ambiguity:
    def __init__(self, bytes ABC, int Ai, int Ci, int Aj, int Cj, int i, int j):
        self._ABC = ABC
        self._Ai = Ai
        self._Ci = Ci
        self._Aj = Aj
        self._Cj = Cj
        self._i = i
        self._j = j   
        self._sig = None
        self._lm = None
############################################################################
    def ABC(self): return self._ABC
    def lm(self): return self._lm 
    def sig(self): return self._sig  
    cdef tuple AC(Ambiguity self): return self._ABC[:self._Ai], self._ABC[self._Ci:],self._ABC[:self._Aj], self._ABC[self._Cj:]
    def ij(self): return self._i, self._j
    def i(self): return self._i
    def j(self): return self._j
    def degree(self): return len(self._ABC)     
############################################################################     
    def _eq_(self,other):
       return self._i == other._i and \
              self._j == other._j and \
              self._Ai == other._Ai and \
              self._Aj == other._Aj and \
              self._Ci == other._Ci and \
              self._Cj == other._Cj and \
              self._ABC == other._ABC
############################################################################                  
    def __ne__(self,other): return not (self == other)
#############################################################################
    def __hash__(self):
       return hash( (self._ABC,self._Ai,self._Ci,self._Aj,self._Cj,self._i,self._j) )
#############################################################################
    def __repr__(self):
            
            Ai, Ci, Aj, Cj = self.AC()
            
            ABC = str(self._ABC)
            Ai = str(Ai)
            Ci = str(Ci)
            Aj = str(Aj)
            Cj = str(Cj)         
            
            return "(" + ABC + ", " + Ai + ", " + Ci + ", " + Aj + ", " + \
            Cj + ", (" + str(self._i) + ", " + str(self._j) + "))"
############################################################################
    def to_poly(self, g):
        cdef bytes Ai,Ci
        Ai,Ci,_,_ = self.AC()      
        return g.lrmul(Ai,Ci)
############################################################################ 
    cdef void to_crit_pair(Ambiguity self, SigPoly gi, SigPoly gj):

        cdef bytes Ai,Ci,Aj,Cj
        cdef Sig s1,s2
        
        Ai,Ci,Aj,Cj = self.AC()
        s1 = gi._sig.lrmul(Ai,Ci)
        s2 = gj._sig.lrmul(Aj,Cj)
                 
        # singular critical pair -> redundant
        if s1 == s2: return
            
        # regular crit pair -> store relevant data in i-part
        self._lm = gi.lm().copy()
        self._lm.lrmul(Ai,Ci)
        
        if s1 < s2:
            self._sig = s2
            self._i = self._j
            self._Ai = self._Aj
            self._Ci = self._Cj
        else:
            self._sig = s1

############################################################################ 
    def __truediv__(self,other):
        r"""
        Divide ``self`` by ``other`` if possible.
        
        Am ambiguity `(ABC, A_i, C_i, A_j, C_j, i ,j)` is divisible by 
        another ambiguity `(A'B'C', A_i', C_i', A_j', C_j', i', j')` if
        there exist `L` and `R` such that `A_j = L A_j' and `C_j = C_j' R`.
        
        INPUT:
        
        - ``other`` -- Ambiguity
        
        OUTPUT:
        
        - `0` if ``other`` does not divides ``self``
        
        - `1` if ``other`` divides ``self`` and `L = R = 1`
        
        - `2` if if ``other`` divides ``self`` and `LR \neq 1`, i.e, ``other``
        properly divides ``self`` 
        
        TESTS::
        
            sage: a = Ambiguity()
        
        
        """  
        
        sAj, sCj = self._Aj, self._Cj
        oAj, oCj = other._Aj, other._Cj
        sdeg = self.degree()
        odeg = other.degree()
        
        start = sAj - oAj
        end = start + odeg
        
        if start < 0 or end > sdeg: return 0
        elif self._ABC[start:end] == other._ABC:
            if sdeg == odeg: return 1
            else: return 2
############################################################################               
    @staticmethod
    cdef list generate_incls(int i, bytes v, int j , bytes w):        
        cdef list amb = []
        cdef int k
        
        if i != j:
            k = v.find(w, 0)
            while k >= 0:
                amb.append( Ambiguity(v,0,len(v),k,k+len(w),i,j) )
                k = v.find(w,k+1)

        return amb
            
############################################################################
    @staticmethod
    cdef list generate_with_tries(prefix_trie, suffix_trie, words, oldlen):
        cdef list amb
        cdef int i,j,k,l,len_m
        cdef bytes m,m_rev,A,B,C,AB,BC
        
        amb = []
        l = len(words)
        
        for i from oldlen <= i < l:
            m = words[i]
            m_rev = m[::-1]
            len_m = len(m)
            
            for k from 1 <= k < len_m:
                # overlap ambiguities with m = AB
                A = m[:-k]
                B = m[-k:]
                amb += [Ambiguity(A+BC,len(A),len(A+BC),0,len_m,j,i) for j,BC in prefix_trie.values(B) if len(BC) > k and j <= i]
    
                # overlap ambiguities with m = BC
                B = m_rev[-k:]
                C = m[k:]            
                amb += [Ambiguity(AB+C,0,len(AB),len(AB)-k,len(AB+C),j,i) for j,AB in suffix_trie.values(B) if len(AB) > k and j < i]
                                    
            # inclusion ambiguities with m = ABC       
            amb += [Ambiguity(m,k-len(B)+1,k+1,0,len(m),j,i) for k,(j,B) in prefix_trie.iter(m) if j < i]
            # inclusion ambiguities with m = B
            amb += flatten([Ambiguity.generate_incls(j,w,i,m) for j,w in enumerate(words[:i]) if len(w) > len_m])          
        
        return amb
  
###########################################################################
    @staticmethod
    def gebauer_moeller(amb):
        amb = sorted(amb,key=lambda a : (a.degree(), a._i, a._Ai))
        idx = 0
                                    
        while idx < len(amb)-1:
            a = amb[idx]
            idx += 1
            amb_tmp = amb[:idx]
            j = a._i
            for aa in amb[idx:]:
                vw = aa / a
                if vw:
                    # divisible and vw != 1, or
                    # divisible and i > j
                    # => throw away aa
                    i = aa._i
                    if vw > 1 or i > j: continue
                
                    # divisible, i == j, vw = 1 -> look at cofactors of ambiguities
                    # if aa_wi > a_wi -> throw away
                    if i == j and vw == 1 and aa._Ai > a._Ai: continue
             
                # not in one of the above cases -> have to keep aa
                amb_tmp.append(aa)
            amb = amb_tmp
                    
        return amb
    
############################################################################
    @staticmethod
    def chain_criterion(amb_old,amb_s,lm_s,s):
                
        amb = []
        set_amb_s = set(amb_s)
        len_lm_s = len(lm_s)
        for a in amb_old:
            keep = True   
            ABC = a.ABC()
            if lm_s in ABC:
                i,j = a.ij()
                Ai,Ci,Aj,Cj = a._Ai, a._Ci, a._Aj, a._Cj
                k = ABC.find(lm_s)
                len_ABC = len(ABC)
                while k > -1:
                    As,Cs = k, k+len_lm_s
                    # check if new ambiguities are trivial
                    is_trivial_i = Ai >= As + len_lm_s or Ci <= Cs - len_lm_s
                    is_trivial_j = Aj >= As + len_lm_s or Cj <= Cs - len_lm_s
                    # compute new ambiguities and shrink them
                    a_is = Ambiguity(ABC,Ai,Ci,As,Cs,i,s).shrink() if not is_trivial_i else None
                    a_js = Ambiguity(ABC,Aj,Cj,As,Cs,j,s).shrink() if not is_trivial_j else None
                    # check if new ambiguities are trivial or contained
                    if (is_trivial_i or a_is in set_amb_s) and (is_trivial_j or a_js in set_amb_s):
                        keep = False
                        break
                    k = ABC.find(lm_s,k+1)
            if keep: amb.append(a)          
        return amb 
############################################################################       
    def shrink(self):
    
        ABC = self.ABC()
        i,j = self.ij()
        Ai,Ci,Aj,Cj = self._Ai, self._Ci, self._Aj, self._Cj
        k = min(Ai,Aj)
        l = max(Ci,Cj)
        
        return Ambiguity(ABC[k:l], Ai-k, Ci-k, Aj-k, Cj-k,i,j) 