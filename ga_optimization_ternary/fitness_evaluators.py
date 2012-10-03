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


def gaussian_pdf(x, mean=0):
    return (1/math.sqrt(2*math.pi))*math.exp(-0.5*(x-mean)*(x-mean))


class FitnessEvaluator():
    
    def __init__(self, fitness, temperature):
        
        all_indices = range(52)
        all_Z = [3, 4, 5, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83]
        self._Z_dict = dict(zip(all_indices, all_Z))
        self._reverse_dict = dict(zip(all_Z, all_indices))
        self._db = M_Database()
        self._fitness = fitness
        self._temp = temperature

    def array_to_score(self, genome):
        genome = self.convert_raw_to_Z(genome)
        gap_dir, gap_ind, heat, vb_dir, cb_dir, vb_ind, cb_ind = self._db.get_data(genome[0], genome[1], genome[2])
        raw_fit = self._fitness(gap_dir, gap_ind, heat, vb_dir, cb_dir, vb_ind, cb_ind)
        # scaled_fit = math.exp(raw_fit/self._temp)
        return raw_fit
                                          
    def convert_Z_to_raw(self, array):
        # note that we are representing as A, X, B to get symmetry between A & B
        A = self._reverse_dict[array[0]]
        B = self._reverse_dict[array[2]]
        X = array[1]
        return (A, B, X)
    
    def convert_raw_to_Z(self, array):
        # note that we are representing as A, X, B to get symmetry between A & B
        A = self._Z_dict[array[0]]
        B = self._Z_dict[array[2]]
        X = array[1]
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
        stab_score = 0
        gap_dir_score = 0
        gap_ind_score = 0
        
        if (gap_dir >= 1.5 and gap_dir <= 3):
            gap_dir_score += 10
        else:
            gap_dir_score += gaussian_pdf(gap_dir, 2.25)
            
        if (gap_ind >= 1.5 and gap_ind <= 3):
            gap_ind_score += 10
        else:
            gap_ind_score += gaussian_pdf(gap_ind, 2.25)
            
        if heat_of_formation <= 0.2:
            stab_score = 10
        else:
            stab_score = 20 * (1-1/(1+math.exp(((-heat_of_formation) + 0.2) * 3.5)))
        
        if vb_dir >= 5.73:
            gap_dir_score += 5
        else:
            gap_dir_score += min(0, 5.73 - vb_dir)
         
        if vb_ind >= 5.73:
            gap_ind_score += 5
        else:
            gap_ind_score += min(0, 5.73 - vb_dir)
             
        if cb_dir <= 4.5:
            gap_dir_score += 5
        else:
            gap_dir_score += min(0, cb_dir - 4.5)
        
        if cb_ind <= 4.5:
            gap_ind_score += 5
        else:
            gap_ind_score += min(0, cb_ind - 4.5)
        
        return max(gap_ind_score, gap_dir_score) + stab_score


def eval_fitness_simple(gap_dir, gap_ind, heat_of_formation, vb_dir, cb_dir, vb_ind, cb_ind):
        stab_score = 0
        gap_dir_score = 0
        gap_ind_score = 0
        
        if (gap_dir >= 1.5 and gap_dir <= 3):
            gap_dir_score += 10
            
        if (gap_ind >= 1.5 and gap_ind <= 3):
            gap_ind_score += 10

        if heat_of_formation <= 0.5:
            stab_score += 5
                    
        if heat_of_formation <= 0.2:
            stab_score += 5
        
        if (vb_dir >= 5.73 and cb_dir <= 4.5):
            gap_dir_score += 10
        
        if (vb_ind >= 5.73 and cb_ind <= 4.5):
            gap_ind_score += 10
        
        return max(gap_ind_score, gap_dir_score) + stab_score


if __name__ == "__main__":
    print eval_fitness_complex(3.01, 3.01, 0.5, 1, 1, 2, 2)
    # print gaussian_pdf(gap_ind, 2.25) * 33
    #fe = FitnessEvaluator()
    #print fe.array_to_score((49,41,20))
    #print [i for i in fe.convert_Z_to_raw((12, 1, 23))]
    #print gaussian_pdf(1)
