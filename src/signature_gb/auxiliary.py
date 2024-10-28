r"""
Useful auxiliary functions
"""
from __future__ import absolute_import

from sage.all import ZZ,prod
import re
from itertools import chain

import sage.parallel.multiprocessing_sage as MP

############################################################################  
def my_pretty_print(string,A):
    r"""
    Pretty-print a string representing a monomial.

    """
    m = string        
    if not m: return '1'
    
    T = A.translator()
    
    X = set(m)
    m = "*".join([m[i:i+1] for i in range(len(m))]) + "*"
    for x in X:
        eq = r"(%s\*){2,}" % x
        m = re.sub(eq, lambda y : str(x) + "^" + str((y.end() - y.start())//2) + "*", m)
    return T(m[:-1],to_internal=False)
############################################################################
def simplify_str(string):
    r"""
    Simplify a string representing a monomial
    
    Remove `*` characters and expand powers `x^n` as `x\dots x`.

    INPUT: 
    
    - ``string`` -- string
    
    OUTPUT: ``string`` with all `*` characters removed and all powers
     `x^n` expanded as `x\dots x`. The string '1' gets simplified to ''.
            
    TESTS::
     
        sage: from OperatorGB import *
        sage: simplify_string('1')
        ''
        sage: simplify_string('a*b*c')
        'abc'
        sage: simplify_string('a^2*b^3*c')
        'aabbbc'
    """
    m = ''
    if string == '1': return m
    else:
        for v in str(string).split('*'):
            j = v.find('^')
            if j != -1:
                m += v[:j] * ZZ(v[j+1:])
            else:
                m += v
        return m
############################################################################
    
def compute_basis(M, maxiter, maxdeg, sig_bound, verbose):
    return M.signature_basis(maxiter=maxiter, maxdeg=maxdeg, sig_bound=sig_bound, verbose=verbose)
        


