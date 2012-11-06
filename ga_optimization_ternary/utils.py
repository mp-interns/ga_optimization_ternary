#!/usr/bin/env python

from __future__ import division
'''
Created on Jul 11, 2012
'''
from ga_optimization_ternary.database import MAX_GOOD_LS, NUM_CANDS

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jul 11, 2012"

import random
import os
import pickle

"""
Hypergeometric distribution

Models drawing objects from a bin. M is total number of objects, n is total number of Type I objects. RV counts number of Type I objects in N drawn without replacement from population.

hypergeom.pmf(k, M, n, N) = choose(n,k)*choose(M-n,N-k)/choose(M,N) for N - (M-n) <= k <= min(m,N)
"""

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
"""

def get_reference_array():
    return [0.0, 901.380952381, 1802.76190476, 2704.14285714, 3605.52380952, 4506.9047619, 5408.28571429, 6309.66666667, 7211.04761905, 8112.42857143, 9013.80952381, 9915.19047619, 10816.5714286, 11717.952381, 12619.3333333, 13520.7142857, 14422.0952381, 15323.4761905, 16224.8571429, 17126.2380952, 18027.6190476]

"""
def _get_prob_bailey(JC, CC, JJ):
    a = (CC + 1)*(CC + CC*CC - 2*CC*JC + JC*(JC - 1)) * comb(CC, JC - 1) / ((JC + 1)*(JC + 2)*comb(CC, JJ))
    
    print a
"""


def _get_prob(k, M, n):
    return n*(M+1)/(k+1)


def get_random_genomes():
    
    filename = "random_initializations.p"
    if os.path.exists(filename):
        with open(filename) as f:
            return pickle.load(f)
    
    gen_length = 2000  # works up to 2000 popsize
    num_gens = 50  # works for fifty iterations
    
    r_init = [0] * num_gens
    
    for gen in range(num_gens):
        for cand in range(gen_length):
            if cand == 0:
                r_init[gen] = []  
            A = random.randint(0, 51)
            B = random.randint(0, 51)
            X = random.randint(0, 6)
            vector = [A, X, B]
            while vector in r_init[gen]:
                A = random.randint(0, 51)
                B = random.randint(0, 51)
                X = random.randint(0, 6)
                vector = [A, X, B]
            r_init[gen].append(vector)
    
    # dump the data
    with open(filename, "wb") as f:
        pickle.dump(r_init, f)
        
if __name__ == "__main__":
    print get_reference_array()[20]
