#!/usr/bin/env python

from __future__ import division
'''
Created on Jul 11, 2012
'''
from ga_optimization_ternary.database import MAX_GOOD, MAX_CAND

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jul 11, 2012"

from scipy.stats import hypergeom

from scipy.misc import comb

"""
Hypergeometric distribution

Models drawing objects from a bin. M is total number of objects, n is total number of Type I objects. RV counts number of Type I objects in N drawn without replacement from population.

hypergeom.pmf(k, M, n, N) = choose(n,k)*choose(M-n,N-k)/choose(M,N) for N - (M-n) <= k <= min(m,N)
"""

def _get_prob(k, M, n):
    #get k right out of M object with n good ones
    topsum = 0.0
    bottomsum = 0.0
    for draws in range(k, M + 1):
        prob_dens = hypergeom.pmf(k, M, n, draws) * k/draws
        topsum += prob_dens * draws
        bottomsum += prob_dens    
    
    print topsum/bottomsum
    return topsum/bottomsum

def get_reference_array():
    return [0, 338.06249999996948, 676.12500000002126, 1014.1875000000172, 1352.2499999999927, 1690.3124999999836, 2028.3749999999909, 2366.4374999999732, 2704.4999999999764, 3042.5624999999914, 3380.6249999999941, 3718.6874999999995, 4056.7499999999932, 4394.8124999999964, 4732.8750000000027, 5070.9374999999982]

def _get_prob_bailey(JC, CC, JJ):
    a = (CC + 1)*(CC + CC*CC - 2*CC*JC + JC*(JC - 1)) * comb(CC, JC - 1) / ((JC + 1)*(JC + 2)*comb(CC, JJ))
    
    print a
    
if __name__ == "__main__":
    _get_prob_bailey(1, 3, 1)
    _get_prob(1, 3, 1)
    """
    a = []
    for i in range(0, MAX_GOOD):
        print i+1
        a.append(_get_prob(i+1, MAX_CAND, 15))
    
    print a
    """