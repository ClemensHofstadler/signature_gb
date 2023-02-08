from __future__ import absolute_import

import itertools
from statistics import median

def flatten(l): return list(itertools.chain(*l))
    
def rebalance(pairs):
    lens = [len(l) for l in pairs]
            
    # no need for rebalancing
    if 2*median(lens) > min(lens) + max(lens): return pairs
        
    C = (min(lens) + max(lens)) // 2
    short_lists = [l for l in pairs if len(l) <= C]
    if len(short_lists) > 1:
        pairs = [l for l in pairs if len(l) > C]
        longer_list = list(itertools.chain(*short_lists))
        longer_list.sort(key = lambda p: p._sig, reverse=True)
        pairs.append(longer_list)
                              
    return pairs
    
