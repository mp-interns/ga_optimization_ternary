#!/usr/bin/env python

'''
Created on Mar 14, 2012
'''
from ga_optimization_ternary.ranked_list_optimization import get_excluded_list

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Mar 14, 2012"

from pymatgen.core.periodic_table import Element
from database import M_Database
import math


def gaussian_pdf(x, mean=0, width = 0.5):
    return (1/math.sqrt(2*math.pi))*math.exp(-width*(x-mean)*(x-mean))


class FitnessEvaluator():
    
    def __init__(self, fitness, temperature):
        
        all_indices = range(52)
        all_Z = [3, 4, 5, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83]
        self._Z_dict = dict(zip(all_indices, all_Z))
        self._reverse_dict = dict(zip(all_Z, all_indices))
        self._db = M_Database()
        self._fitness = fitness
        self._temp = temperature
        self._exclusions = get_excluded_list()

    def array_to_score(self, genome):
        genome = self.convert_raw_to_Z(genome)
        gap_dir, gap_ind, heat, vb_dir, cb_dir, vb_ind, cb_ind = self._db.get_data(genome[0], genome[1], genome[2])
        if self._fitness.__name__ == "eval_fitness_simple_exclusion":
            if genome in self._exclusions:
                return 0
            return eval_fitness_simple(gap_dir, gap_ind, heat, vb_dir, cb_dir, vb_ind, cb_ind)
        
        if self._fitness.__name__ == "eval_fitness_complex_exclusion":
            if genome in self._exclusions:
                return 0
            return eval_fitness_complex(gap_dir, gap_ind, heat, vb_dir, cb_dir, vb_ind, cb_ind)     

        if self._fitness.__name__ == "eval_fitness_complex_product_exclusion":
            if genome in self._exclusions:
                return 0
            return eval_fitness_complex_product(gap_dir, gap_ind, heat, vb_dir, cb_dir, vb_ind, cb_ind)   
                    
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
        elif gap_dir == 0:
            gap_dir_score += 0
        else:
            gap_dir_score += 33 * gaussian_pdf(gap_dir, 2.25)
            
        if (gap_ind >= 1.5 and gap_ind <= 3):
            gap_ind_score += 10
        elif gap_ind == 0:
            gap_ind_score += 0
        else:
            gap_ind_score += 33 * gaussian_pdf(gap_ind, 2.25)
            
        if heat_of_formation <= 0.2:
            stab_score = 10
        else:
            stab_score = 20 * (1-1/(1+math.exp(((-heat_of_formation) + 0.2) * 3.5)))
        
        if vb_dir >= 5.73:
            gap_dir_score += 5
        else:
            distance = (5.73 - vb_dir) * 5
            gap_dir_score += 10 * (1-1/(1+math.exp(-distance)))
         
        if vb_ind >= 5.73:
            gap_ind_score += 5
        else:
            distance = (5.73 - vb_ind) * 5
            gap_ind_score += 10 * (1-1/(1+math.exp(-distance)))
             
        if cb_dir <= 4.5:
            gap_dir_score += 5
        else:
            distance = (cb_dir - 4.5) * 5
            gap_dir_score += 10 * (1-1/(1+math.exp(-distance)))
        
        if cb_ind <= 4.5:
            gap_ind_score += 5
        else:
            distance = (cb_ind - 4.5) * 5
            gap_ind_score += 10 * (1-1/(1+math.exp(-distance)))
        
        return max(gap_ind_score, gap_dir_score) + stab_score


def eval_fitness_complex_product(gap_dir, gap_ind, heat_of_formation, vb_dir, cb_dir, vb_ind, cb_ind):
        stab_score = 0
        gap_dir_score = 0
        gap_ind_score = 0
        
        if (gap_dir >= 1.5 and gap_dir <= 3):
            gap_dir_score += 10
        elif gap_dir == 0:
            gap_dir_score += 0
        else:
            gap_dir_score += 33 * gaussian_pdf(gap_dir, 2.25)
            
        if (gap_ind >= 1.5 and gap_ind <= 3):
            gap_ind_score += 10
        elif gap_ind == 0:
            gap_ind_score += 0
        else:
            gap_ind_score += 33 * gaussian_pdf(gap_ind, 2.25)
            
        if heat_of_formation <= 0.2:
            stab_score = 10
        else:
            stab_score = 20 * (1-1/(1+math.exp(((-heat_of_formation) + 0.2) * 3.5)))
        
        if vb_dir >= 5.73:
            gap_dir_score += 5
        else:
            distance = (5.73 - vb_dir) * 5
            gap_dir_score += 10 * (1-1/(1+math.exp(-distance)))
         
        if vb_ind >= 5.73:
            gap_ind_score += 5
        else:
            distance = (5.73 - vb_ind) * 5
            gap_ind_score += 10 * (1-1/(1+math.exp(-distance)))
             
        if cb_dir <= 4.5:
            gap_dir_score += 5
        else:
            distance = (cb_dir - 4.5) * 5
            gap_dir_score += 10 * (1-1/(1+math.exp(-distance)))
        
        if cb_ind <= 4.5:
            gap_ind_score += 5
        else:
            distance = (cb_ind - 4.5) * 5
            gap_ind_score += 10 * (1-1/(1+math.exp(-distance)))
            
        return max(gap_ind_score, gap_dir_score) * stab_score * 0.15
    

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
        
        if (vb_dir >= 5.73):
            gap_dir_score += 5
        
        if (cb_dir <= 4.5):
            gap_dir_score += 5
            
        if (vb_ind >= 5.73):
            gap_ind_score += 5
            
        if (cb_ind <= 4.5):
            gap_ind_score += 5
            
        return max(gap_ind_score, gap_dir_score) + stab_score

def eval_fitness_simple_exclusion():
    pass

def eval_fitness_complex_exclusion():
    pass

def eval_fitness_complex_product_exclusion():
    pass


def eval_fitness_complex_oxide_shield(gap_dir, gap_ind, heat_of_formation, vb_dir, cb_dir, vb_ind, cb_ind):
        stab_score = 0
        gap_score = 0
        gap_dir_score = 0
        gap_ind_score = 0
        
        if gap_dir >= 3:
            gap_score += 10
        elif gap_dir == 0:
            gap_score = 0
        else:
            distance = (3 - gap_dir) * 5
            gap_score += 20 * (1-1/(1+math.exp(-distance)))
        
        if heat_of_formation <= 0.2:
            stab_score = 10
        else:
            stab_score = 20 * (1-1/(1+math.exp(((-heat_of_formation) + 0.2) * 3.5)))
        
        if (vb_dir >= 6.2 and vb_dir <= 7):
            gap_dir_score += 10
        else:
            gap_dir_score += 55 * gaussian_pdf(vb_dir, 6.6, 5)
        
        if (vb_ind >= 6.2 and vb_ind <= 7):
            gap_ind_score += 10
        else:
            gap_ind_score += 55 * gaussian_pdf(vb_ind, 6.6, 5)
        
        return gap_score + stab_score + max(gap_ind_score, gap_dir_score)


def eval_fitness_complex_product_oxide_shield(gap_dir, gap_ind, heat_of_formation, vb_dir, cb_dir, vb_ind, cb_ind):
        stab_score = 0
        gap_score = 0
        gap_dir_score = 0
        gap_ind_score = 0
        
        if gap_dir >= 3:
            gap_score += 10
        elif gap_dir == 0:
            gap_score = 0
        else:
            distance = (3 - gap_dir) * 5
            gap_score += 20 * (1-1/(1+math.exp(-distance)))
        
        if heat_of_formation <= 0.2:
            stab_score = 10
        else:
            stab_score = 20 * (1-1/(1+math.exp(((-heat_of_formation) + 0.2) * 3.5)))
        
        if (vb_dir >= 6.2 and vb_dir <= 7):
            gap_dir_score += 10
        else:
            gap_dir_score += 55 * gaussian_pdf(vb_dir, 6.6, 5)
        
        if (vb_ind >= 6.2 and vb_ind <= 7):
            gap_ind_score += 10
        else:
            gap_ind_score += 55 * gaussian_pdf(vb_ind, 6.6, 5)
        
        return gap_score/10 * stab_score/10 * max(gap_ind_score, gap_dir_score)/10 * 30


def eval_fitness_simple_oxide_shield(gap_dir, gap_ind, heat_of_formation, vb_dir, cb_dir, vb_ind, cb_ind):
        stab_score = 0
        gap_score = 0
        gap_dir_score = 0
        gap_ind_score = 0
        
        if gap_dir >= 3:
            gap_score += 10
            
        if heat_of_formation <= 0.5:
            stab_score += 5
                    
        if heat_of_formation <= 0.2:
            stab_score += 5
        
        if (vb_dir >= 6.2 and vb_dir <= 7):
            gap_dir_score += 10
        
        if (vb_ind >= 6.2 and vb_ind <= 7):
            gap_ind_score += 10
        
        return gap_score + stab_score + max(gap_ind_score, gap_dir_score)


def score (cb_dir):
        if cb_dir <= 4.5:
            return 5
        else:
            return max(0, 5 - (cb_dir - 4.5))


if __name__ == "__main__":
    #print score(4.8)
    print eval_fitness_complex_oxide_shield(2, 7.00001, 0.5, 6.199999, 1, 2, 2)
    # print gaussian_pdf(gap_ind, 2.25) * 33
    #fe = FitnessEvaluator()
    #print fe.array_to_score((49,41,20))
    #print [i for i in fe.convert_Z_to_raw((12, 1, 23))]
    #print gaussian_pdf(1)
