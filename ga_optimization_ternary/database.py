#!/usr/bin/env python

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

import fitness_evaluators


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
 
   
if __name__ == "__main__":
    m_db = M_Database()
    
    i = 0
    for A in m_db._all_data:
        for B in m_db._all_data[A]:
            data = m_db._all_data[A][B][False]
            i = i+1
            if fitness_evaluators.eval_fitness_simple(data[0],data[1],data[2], data[3], data[4], data[5], data[6]) >= 30:
                print A, B, Element.from_Z(A), Element.from_Z(B)
    
    print i