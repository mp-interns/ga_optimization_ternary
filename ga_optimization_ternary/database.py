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
import os
import pickle

GOOD_CANDS_LS = [(3,23,0), (11,51,0), (12,73,1), (20,32,0), (20,50,0), (20,73,1), (38,32,0), (38,50,0), (38,73,1), (39,73,2), (47,41,0), (50,22,0), (55,41,0), (56,31,4), (56,49,4), (56,50,0), (56,73,1), (57,22,1), (57,73,2), (82,31,4)]
MAX_GOOD_LS = len(GOOD_CANDS_LS)
NUM_CANDS = 18928


class M_Database():
    
    def __init__(self):
         
        filename = "data_process.p"
        if os.path.exists(filename):
            with open(filename) as f:
                self._all_data = pickle.load(f)
        
        else:
            connection = pymongo.Connection('localhost', 12345)
            self._db = connection.unc.data_process
            self._all_data = defaultdict(lambda: defaultdict(dict))
            
            for data in self._db.find():
                val = (data['gllbsc_dir-gap'], data['gllbsc_ind-gap'], data['heat_of_formation_all'], data['VB_dir'], data['CB_dir'], data['VB_ind'], data['CB_ind'])
                self._all_data[Element(data['A']).Z][Element(data['B']).Z][data['anion_idx']] = val
            
            # dump the data
            with open(filename, "wb") as f:
                pickle.dump(dict(self._all_data), f)
        
    def get_data(self, A, B, anion_idx):
        return self._all_data[A][B][anion_idx]
    

class Stats_Database():
    def __init__(self, clear=False, extension=""):
        connection = pymongo.Connection('localhost', 12345)
        self._stats_raw = connection.unc["stats_raw" + extension]
        self._stats_process = connection.unc["stats_process" + extension]
        if clear:
            self._stats_raw.remove()
            self._stats_process.remove()
    
    def add_stat_raw(self, ps, stat, iteration):
        doc = {}
        doc['unique_key'] = ps.unique_key()
        doc['iteration'] = iteration
        if iteration == 0:
            doc['parameters'] = ps.to_dict()
        
        it_data = []
        for gen in range(len(stat.generation_ngood)):
            it_data.append({"n_cand":stat.generation_ncandidates[gen], "n_good":stat.generation_ngood[gen]})
        
        doc['data'] = it_data
        doc['num_breakouts'] = stat.num_breakouts
                
        self._stats_raw.insert(doc)
    
    def process_stats_new(self):
        num_iterations = 15
        for key in self._stats_raw.distinct("unique_key"):
            
            #TODO: fix me
            if self._stats_raw.find({"unique_key": key}).count() >= num_iterations:
                # we have a good param set ... go through the iterations
                doc = {}
                it_ng_nc = np.zeros((num_iterations, MAX_GOOD_LS + 1))  # [iteration, numgood] = (# candidates needed)
                
                for it in range(num_iterations):
                    expt = self._stats_raw.find_one({"unique_key":key, "iteration":it})
                    if it == 0:
                        # remove any previous document
                        self._stats_process.remove({"unique_key": expt['unique_key']})
                        doc['parameters'] = expt['parameters']
                        doc['unique_key'] = expt['unique_key']
                        print doc
                       
                    ng = 0  # number good that we're trying to initialize
                    for gen in expt['data']:
                        while gen['n_good'] >= ng:
                            it_ng_nc[it][ng] = gen['n_cand']
                            ng += 1
                    
                ng_it_nc = np.transpose(it_ng_nc)  # [numgood, iteration] = (# candidates needed)
                ng_avg = [np.average(ng_it_nc[idx]) for idx in range(len(ng_it_nc))]  # [numgood] = (avg # of candidates needed)
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
                doc['all'] = ng_avg[MAX_GOOD_LS]  # shorthand, avg number of candidates needed to get all good cands
                doc['ten'] = ng_avg[10]  # shorthand, avg number of candidates needed to get ten good cands
                doc['fifteen'] = ng_avg[15]  # shorthand, avg number of candidates needed to get 15 good cands (about 2/3)
                print doc
                self._stats_process.insert(doc)
        
    def process_stats(self):
        raise ValueError("This must be updated...")
        for expt in self._stats_raw.find():
            
            # remove any previous document
            self._stats_process.remove({"unique_key": expt['unique_key']})
            
            # initalize doc
            doc = {}
            doc['parameters'] = expt['parameters']
            doc['unique_key'] = expt['unique_key']
            
            # initialize data
            num_iterations = len(expt['iterations'])
            it_ng_nc = np.zeros((num_iterations, MAX_GOOD_LS + 1))  # [iteration, numgood] = (# candidates needed)
            
            for it, dat in enumerate(expt['iterations']):
                ng = 0  # number good that we're trying to initialize
                for gen in dat:
                    while gen['n_good'] >= ng:
                        it_ng_nc[it][ng] = gen['n_cand']
                        ng += 1
                        
            ng_it_nc = np.transpose(it_ng_nc)  # [numgood, iteration] = (# candidates needed)
            
            ng_avg = [np.average(ng_it_nc[idx]) for idx in range(len(ng_it_nc))]  # [numgood] = (avg # of candidates needed)
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
            doc['all'] = ng_avg[MAX_GOOD_LS]  # shorthand, avg number of candidates needed to get all good cands
            doc['ten'] = ng_avg[10]  # shorthand, avg number of candidates needed to get ten good cands
            doc['fifteen'] = ng_avg[15]  # shorthand, avg number of candidates needed to get 15 good cands (about 2/3)
            
            print doc
            self._stats_process.insert(doc)

                
class InitializationDB():
    def __init__(self):
        filename = "random_initializations.p"
        with open(filename) as f:
            self.random_inits = pickle.load(f)
    
    def get_initial_list(self, iteration):
        return self.random_inits[iteration]

    
if __name__ == "__main__":
    #s_db = Stats_Database()
    #s_db.process_stats_new()
    
    m_db = M_Database()
    from ga_optimization_ternary.fitness_evaluators import eval_fitness_simple, eval_fitness_complex
    hits = 0
    for A in m_db._all_data:
        for B in m_db._all_data[A]:
            for anion in m_db._all_data[A][B]:
                data = m_db._all_data[A][B][anion]
                if eval_fitness_complex(data[0],data[1],data[2], data[3], data[4], data[5], data[6]) >= 30:
                    #print '('+Element.from_Z(A).symbol+","+Element.from_Z(B).symbol+","+str(anion)+"),",
                    print '[{}, {}, {}],'.format(A, anion, B),
                    hits +=1
    
    print hits
