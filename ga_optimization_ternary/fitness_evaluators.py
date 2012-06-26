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

from pymatgen.core.periodic_table import Element
from database import M_Database
import math

class FitnessEvaluatorZ():
    
    def __init__(self, fitness=None):
        
        all_indices = range(52)
        all_Z = [3, 4, 5, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83]
        self._Z_dict = dict(zip(all_indices, all_Z))
        self._reverse_dict = dict(zip(all_Z, all_indices))
        self._db = M_Database()
        self._fitness = fitness if fitness else eval_fitness_simple

    def array_to_score(self, genome):
        genome = self.convert_raw_to_Z(genome)
        gap_dir, gap_ind, heat, vb_dir, cb_dir, vb_ind, cb_ind = self._db.get_data(genome[0], genome[1], genome[2])
        return self._fitness(gap_dir, gap_ind, heat, vb_dir, cb_dir, vb_ind, cb_ind)
                                          
    def convert_Z_to_raw(self, array):
        # note that we are representing as A, X, B to get symmetry between A & B
        A = self._reverse_dict[array[0]]
        B = self._reverse_dict[array[2]]
        X = 0 if not array[1] else len(self._Z_dict) - 1
        return (A, B, X)
    
    def convert_raw_to_Z(self, array):
        # note that we are representing as A, X, B to get symmetry between A & B
        A = self._Z_dict[array[0]]
        B = self._Z_dict[array[2]]
        cutoff_idx = len(self._Z_dict) / 2
        X = True if array[1] >= cutoff_idx else False
        return (A, B, X)

'''
class FitnessEvaluatorElectronegativity():
    
    def __init__(self, fitness=None):
        all_indices = range(52)
        all_Z = [3, 4, 5, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83]
        all_Z.sort(key=lambda i: Element.from_Z(i).X)
        self._Z_dict = dict(zip(all_indices, all_Z))
        self._reverse_dict = dict(zip(all_Z, all_indices))
        self._db = M_Database()
        self._fitness = fitness if fitness else eval_fitness_simple

    def array_to_score(self, genome):
        el_A = self._Z_dict[genome[0]]
        el_B = self._Z_dict[genome[1]]
        containsN = False
        gap_dir, gap_ind, heat = self._db.get_data(el_A, el_B, containsN)
        return self._fitness(gap_dir, gap_ind, heat)
                                  
    def convert_Z_to_raw(self, array):
        return tuple([self._reverse_dict[i] for i in array])
    
    def convert_raw_to_Z(self, array):
        return tuple([self._Z_dict[i] for i in array])

'''
def eval_fitness_complex(gap_dir, gap_ind, heat_of_formation, vb_dir, cb_dir, vb_ind, cb_ind):
        raise ValueError("Not updated for handle band edges")
        score_1 = 0
        
        if (gap_dir >= 1) or (gap_ind >= 1):
            score_1 = (max(gaussian_pdf(gap_dir, 2.25), gaussian_pdf(gap_ind, 2.25)) * 2.5)
        
        score_2 = 1-1/(1+math.exp(-heat_of_formation))
        
        #print gap_dir, gap_ind, heat_of_formation, (max(gaussian_pdf(gap_dir, 2.25), gaussian_pdf(gap_ind, 2.25)) * 2.5), 1-1/(1+math.exp(-heat_of_formation)), score
        return score_1 * score_2

def eval_fitness_simple(gap_dir, gap_ind, heat_of_formation, vb_dir, cb_dir, vb_ind, cb_ind):
        score_1 = 0
        
        if (gap_dir >= 1.5 and gap_dir <= 3) or (gap_ind >= 1.5 and gap_ind <= 3):
            score_1 = 10
        
        score_2 = 10 + 0.2 - heat_of_formation
        
        score_3 = 0
        # print gap_dir, gap_ind, heat_of_formation, vb_dir, cb_dir, vb_ind, cb_ind
        '''
        if (vb_dir > 5.73 or vb_ind > 5.73):
            score_3 += 2.5
        
        if (cb_dir < 4.5 or cb_ind < 4.5:
            score_3 += 2.5
        '''
        if (vb_dir > 5.73 and cb_dir < 4.5) or (vb_ind > 5.73 and cb_ind < 4.5):
            score_3 += 10
        
        #print gap_dir, gap_ind, heat_of_formation, (max(gaussian_pdf(gap_dir, 2.25), gaussian_pdf(gap_ind, 2.25)) * 2.5), 1-1/(1+math.exp(-heat_of_formation)), score
        return score_1 + score_2 + score_3

def gaussian_pdf(x, mean=0):
    return (1/math.sqrt(2*math.pi))*math.exp(-0.5*(x-mean)*(x-mean))
    

if __name__ == "__main__":
    fe = FitnessEvaluatorZ()
    #print fe.array_to_score((49,41,20))
    print [i for i in fe.convert_Z_to_raw((12, 73, True))]
    print [i for i in fe.convert_Z_to_raw((20, 73, True))]
    print [i for i in fe.convert_Z_to_raw((38, 73, True))]
    print [i for i in fe.convert_Z_to_raw((56, 73, True))]
    print [i for i in fe.convert_Z_to_raw((12, 73, True))]
    print [i for i in fe.convert_Z_to_raw((12, 73, True))]
    print [i for i in fe.convert_Z_to_raw((12, 73, True))]
    #print gaussian_pdf(1)
    '''
12 73
20 73
38 73
56 73
57 22
82 73
    '''
