#!/usr/bin/env python

from __future__ import division

'''
Created on Mar 14, 2012
'''

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Mar 14, 2012"

import pymongo
from collections import defaultdict
from pymatgen.core.periodic_table import Element
from scipy.interpolate import interp1d
import numpy as np

MAX_GOOD = 15
MAX_CAND = 5408

class M_Database():
    
    def __init__(self):
        connection = pymongo.Connection('localhost', 12345)
        self._db = connection.unc.data_process
        self._all_data = defaultdict(lambda: defaultdict(dict))
        
        for data in self._db.find():
            val = (data['gllbsc_dir-gap'], data['gllbsc_ind-gap'], data['heat_of_formation_all'], data['VB_dir'], data['CB_dir'], data['VB_ind'], data['CB_ind'])
            self._all_data[Element(data['A']).Z][Element(data['B']).Z][data['contains_anion']] = val
        
    def get_data(self, A, B, contains_N):
        return self._all_data[A][B][contains_N]


class Stats_Database():
    def __init__(self, clear=False):
        connection = pymongo.Connection('localhost', 12345)
        self._stats_raw = connection.unc.stats_raw
        self._stats_process = connection.unc.stats_process
        if clear:
            self._stats_raw.remove()
            self._stats_process.remove()
    
    def add_stats_raw(self, ps, stats):
        doc = {}
        self._stats_raw.remove({"unique_key":ps.unique_key()})
        doc['unique_key'] = ps.unique_key()
        doc['parameters'] = ps.to_dict()
        doc['iterations'] = []
        for stat in stats:
            it_data = []
            for gen in range(len(stat.generation_ngood)):
                it_data.append({"n_cand":stat.generation_ncandidates[gen], "n_good":stat.generation_ngood[gen]})
            doc['iterations'].append(it_data)
                
        self._stats_raw.insert(doc)
    
    def process_stats(self):
        
        for expt in self._stats_raw.find():
            
            # remove any previous document
            self._stats_process.remove({"unique_key": expt['unique_key']})
            
            # initalize doc
            doc = {}
            doc['parameters'] = expt['parameters']
            doc['unique_key'] = expt['unique_key']
            
            # initialize data
            num_iterations = len(expt['iterations'])
            it_ng_nc = np.zeros((num_iterations, MAX_GOOD+1))  # [iteration, numgood] = (# candidates needed)
            
            for it, dat in enumerate(expt['iterations']):
                ng = 0  # number good that we're trying to initialize
                for gen in dat:
                    while gen['n_good'] >= ng:
                        it_ng_nc[it][ng] = gen['n_cand']
                        ng += 1
                        
            ng_it_nc = np.transpose(it_ng_nc)  # [numgood, iteration] = (# candidates needed)
            
            ng_avg =  [np.average(ng_it_nc[idx]) for idx in range(len(ng_it_nc))]  # [numgood] = (avg # of candidates needed)
            ng_stdev = [np.std(ng_it_nc[idx]) for idx in range(len(ng_it_nc))] # [numgood] = (stdev # of candidates needed)
            ng_min = [np.min(ng_it_nc[idx]) for idx in range(len(ng_it_nc))] 
            ng_max = [np.max(ng_it_nc[idx]) for idx in range(len(ng_it_nc))]
            ng_range = [ng_max[idx] - ng_min[idx] for idx in range(len(ng_max))]

            doc['ng_it_nc'] = ng_it_nc.tolist()
            doc['ng_avg'] = ng_avg
            doc['ng_stdev'] = ng_stdev
            doc['ng_min'] = ng_min
            doc['ng_max'] = ng_max
            doc['ng_range'] = ng_range
            doc['all'] = ng_avg[15]  # shorthand, avg number of candidates needed to get all good cands
            doc['ten'] = ng_avg[10]  # shorthand, avg number of candidates needed to get ten good cands
            
            print doc
            self._stats_process.insert(doc)
                
        
if __name__ == "__main__":
    s_db = Stats_Database()
    s_db.process_stats()
    """
    m_db = M_Database()
    
    i = 0
    for A in m_db._all_data:
        for B in m_db._all_data[A]:
            data = m_db._all_data[A][B][False]
            i = i+1
            if fitness_evaluators.eval_fitness_simple(data[0],data[1],data[2], data[3], data[4], data[5], data[6]) >= 30:
                print A, B, Element.from_Z(A), Element.from_Z(B)
    
    print i
    """