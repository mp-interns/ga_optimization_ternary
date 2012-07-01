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

import fitness_evaluators

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
        # THIS CODE IS ALL UGLY AND NON PYTHONIC. IT's LATE.
        
        for expt in self._stats_raw.find():
            doc = {}
            # remove any previous document
            self._stats_process.remove({"unique_key": expt['unique_key']})
        
            my_big_array = []  # each idx is one iteration. The val is an ngood array
            my_final_array = []
            
            doc['parameters'] = expt['parameters']
            doc['unique_key'] = expt['unique_key']
            
            for it in expt['iterations']:
                ngood = []  # index is n_candidates tried, val is n_candidates found
                
                x = [0]
                y = [0]
                for gen in it:
                    if gen['n_good'] > y[-1]:
                        x.append(gen['n_cand'])
                        y.append(gen['n_good'])
                
                if not x[-1] == MAX_CAND: 
                    x.append(MAX_CAND)
                    y.append(MAX_GOOD)
                f = interp1d(x, y)
                for idx in range(MAX_CAND + 1):
                    ngood.append(float(f(idx)))
                
                my_big_array.append(ngood)
            
            # my_big_array is initialized
            # get the final_array
            sum_array = []
            num_iterations = len(my_big_array)
            #not pythonic, who's looking
            for cand_no, t in enumerate(my_big_array[0]):
                m_sum = 0
                for it_no, t in enumerate(my_big_array):
                    m_sum += my_big_array[it_no][cand_no]
                    #print m_sum, it_no
                #print m_sum
                sum_array.append(m_sum)
                
            final_array = [val / num_iterations for val in sum_array]
            
            doc['averaged_ngood'] = final_array
            
            tries_needed = []
            for idx in range(MAX_CAND + 1):
                for tries, val in enumerate(final_array):
                    if val >= idx:
                        tries_needed.append(tries)
                        break
            
            print tries_needed
            doc['tries_needed'] = tries_needed
            doc['all'] = tries_needed[MAX_GOOD]
            doc['ten'] = tries_needed[10]
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