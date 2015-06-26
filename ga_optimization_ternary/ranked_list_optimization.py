#!/usr/bin/env python
from __future__ import division
'''
Created on Jul 23, 2012
'''

from pymatgen.core.periodic_table import Element
import math
import pickle
import os
from ga_optimization_ternary.database import GOOD_CANDS_LS, GOOD_CANDS_OS

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jul 23, 2012"

"""
def get_ranked_list_goldschmidt():
    #TODO: get the oxidation state right!!!
    
    filename = "goldschmidt_rank.p"
    if os.path.exists(filename):
        with open(filename) as f:
            return pickle.load(f)
    
    all_AB = FitnessEvaluator(eval_fitness_simple, 10)._reverse_dict.keys()
    print 'generating goldschmidt ranks...'
    cand_score = {}  # dictionary of cand_tuple:score. a high score is BAD
    for a in all_AB:
        for b in all_AB:
            for x in range(7):
                r_a = Element.from_Z(a).average_ionic_radius  # TODO: get the correct oxidation state!
                r_b = Element.from_Z(b).average_ionic_radius
                r_x = None
                if x == 0:
                    r_x = Element("O").ionic_radii[-2]
                elif x == 1:
                    r_x = Element("O").ionic_radii[-2] * 2/3 + Element("N").ionic_radii[-3] * 1/3
                elif x == 2:
                    r_x = Element("O").ionic_radii[-2] * 1/3 + Element("N").ionic_radii[-3] * 2/3
                elif x == 3:
                    r_x = Element("N").ionic_radii[-3]
                elif x == 4:
                    r_x = Element("O").ionic_radii[-2] * 2/3 + Element("F").ionic_radii[-1] * 1/3
                elif x == 5:
                    r_x = Element("O").ionic_radii[-2] * 1/3 + Element("F").ionic_radii[-1] * 1/3 + Element("N").ionic_radii[-3] * 1/3
                elif x == 6:
                    r_x = Element("O").ionic_radii[-2] * 2/3 + Element("S").ionic_radii[-2] * 1/3
                
                goldschmidt = (r_a + r_x)/(math.sqrt(2) *(r_b+r_x))
                score = abs(goldschmidt - 1)  # a high score is bad, like golf
                cand_score[(a, b, x)] = score
            
    results = sorted(cand_score, key=cand_score.get)
    with open(filename, "wb") as f:
        pickle.dump(results, f)
    return results


def get_ranked_list_goldschmidt_smart():
    
    # Despite this trying to be smart, it actually performs worse...
    
    filename = "goldschmidt_rank_smart.p"
    if os.path.exists(filename):
        with open(filename) as f:
            return pickle.load(f)
    
    all_AB = FitnessEvaluator(eval_fitness_simple)._reverse_dict.keys()
    print 'generating goldschmidt ranks...'
    cand_score = {}  # dictionary of cand_tuple:score. a high score is BAD
    for a in all_AB:
        for b in all_AB:
            for x in (True, False):
                el_a = Element.from_Z(a)
                el_b = Element.from_Z(b)
                scores = []
                for a_oxi in el_a.oxidation_states:
                    for b_oxi in el_b.oxidation_states:
                        score = -1
                            
                        r_a = el_a.ionic_radii[a_oxi] if a_oxi in el_a.ionic_radii else el_a.average_ionic_radius
                        r_b = el_b.ionic_radii[b_oxi] if b_oxi in el_b.ionic_radii else el_b.average_ionic_radius

                        r_x = None
                        if x:
                            r_x = Element("O").ionic_radii[-2] * 2/3 + Element("N").ionic_radii[-3] * 1/3
                        else:
                            r_x = Element("O").ionic_radii[-2]
            
                        goldschmidt = (r_a + r_x)/(math.sqrt(2) *(r_b+r_x))
                        score = abs(goldschmidt - 1)  # a high score is bad, like golf
                        scores.append(score)
                
                cand_score[(a, b, x)] = min(scores)
    
    results = sorted(cand_score, key=cand_score.get)
    with open(filename, "wb") as f:
        pickle.dump(results, f)
    return results
"""

def get_excluded_list():
    
    filename = "excluded_compounds.p"
    if os.path.exists(filename):
        with open(filename) as f:
            return pickle.load(f)
    
    from ga_optimization_ternary.fitness_evaluators import FitnessEvaluator, eval_fitness_simple
    all_AB = FitnessEvaluator(eval_fitness_simple, 10)._reverse_dict.keys()
    print 'generating exclusion ranks...'
    exclusions = []
    
    for a in all_AB:
        for b in all_AB:
            for x in range(7):
                
                #nelectrons must be even
                ne_a = Element.from_Z(a).Z
                ne_b = Element.from_Z(b).Z
                ne_x = None
                if x == 0:
                    ne_x = Element("O").Z * 3
                elif x == 1:
                    ne_x = Element("O").Z * 2 + Element("N").Z
                elif x == 2:
                    ne_x = Element("O").Z + Element("N").Z * 2
                elif x == 3:
                    ne_x = Element("N").Z
                elif x == 4:
                    ne_x = Element("O").Z * 2 + Element("F").Z
                elif x == 5:
                    ne_x = Element("O").Z + Element("F").Z + Element("N").Z
                elif x == 6:
                    ne_x = Element("O").Z * 2 + Element("S").Z
                
                #modify the score based on charge-balance
                even_found = False
                el_a = Element.from_Z(a)
                el_b = Element.from_Z(b)
                val_x = 0
                if x == 0:
                    val_x = -2 * 3
                elif x == 1:
                    val_x = (-2 * 2) + (-3 * 1)
                elif x == 2:
                    val_x = (-2 * 1) + (-3 * 2)
                elif x == 3:
                    val_x = -3 * 3
                elif x == 4:
                    val_x = (-2 * 2) + (-1 * 1)
                elif x == 5:
                    val_x = (-2 * 1) + (-1 * 1) + (-3 * 1)
                elif x == 6:
                    val_x = (-2 * 2) + (-2 * 1)
                
                for a_oxi in el_a.oxidation_states:
                    for b_oxi in el_b.oxidation_states:
                        if (ne_a + ne_b + ne_x) % 2 == 0 and (a_oxi + b_oxi + val_x) == 0:
                            even_found = True
                
                if not even_found:
                    exclusions.append((a, b, x))
                
    with open(filename, "wb") as f:
        pickle.dump(exclusions, f)
    return exclusions

def get_ranked_list_goldschmidt_halffill():
    #TODO: get the oxidation state right!!!
    
    filename = "goldschmidt_rank_halffill.p"
    if os.path.exists(filename):
        with open(filename) as f:
            return pickle.load(f)
    
    from ga_optimization_ternary.fitness_evaluators import FitnessEvaluator, eval_fitness_simple
    all_AB = FitnessEvaluator(eval_fitness_simple, 10)._reverse_dict.keys()
    print 'generating goldschmidt ranks...'
    cand_score = {}  # dictionary of cand_tuple:score. a high score is BAD
    
    for a in all_AB:
        for b in all_AB:
            for x in range(7):
                r_a = Element.from_Z(a).average_ionic_radius  # TODO: get the correct oxidation state!
                r_b = Element.from_Z(b).average_ionic_radius
                r_x = None
                if x == 0:
                    r_x = Element("O").ionic_radii[-2]
                elif x == 1:
                    r_x = Element("O").ionic_radii[-2] * 2/3 + Element("N").ionic_radii[-3] * 1/3
                elif x == 2:
                    r_x = Element("O").ionic_radii[-2] * 1/3 + Element("N").ionic_radii[-3] * 2/3
                elif x == 3:
                    r_x = Element("N").ionic_radii[-3]
                elif x == 4:
                    r_x = Element("O").ionic_radii[-2] * 2/3 + Element("F").ionic_radii[-1] * 1/3
                elif x == 5:
                    r_x = Element("O").ionic_radii[-2] * 1/3 + Element("F").ionic_radii[-1] * 1/3 + Element("N").ionic_radii[-3] * 1/3
                elif x == 6:
                    r_x = Element("O").ionic_radii[-2] * 2/3 + Element("S").ionic_radii[-2] * 1/3
                
                goldschmidt = (r_a + r_x)/(math.sqrt(2) *(r_b+r_x))
                score = abs(goldschmidt - 1)  # a high score is bad, like golf
                
                #nelectrons must be even
                ne_a = Element.from_Z(a).Z
                ne_b = Element.from_Z(b).Z
                ne_x = None
                if x == 0:
                    ne_x = Element("O").Z * 3
                elif x == 1:
                    ne_x = Element("O").Z * 2 + Element("N").Z
                elif x == 2:
                    ne_x = Element("O").Z + Element("N").Z * 2
                elif x == 3:
                    ne_x = Element("N").Z
                elif x == 4:
                    ne_x = Element("O").Z * 2 + Element("F").Z
                elif x == 5:
                    ne_x = Element("O").Z + Element("F").Z + Element("N").Z
                elif x == 6:
                    ne_x = Element("O").Z * 2 + Element("S").Z
                
                #modify the score based on charge-balance
                even_found = False
                el_a = Element.from_Z(a)
                el_b = Element.from_Z(b)
                val_x = 0
                if x == 0:
                    val_x = -2 * 3
                elif x == 1:
                    val_x = (-2 * 2) + (-3 * 1)
                elif x == 2:
                    val_x = (-2 * 1) + (-3 * 2)
                elif x == 3:
                    val_x = -3 * 3
                elif x == 4:
                    val_x = (-2 * 2) + (-1 * 1)
                elif x == 5:
                    val_x = (-2 * 1) + (-1 * 1) + (-3 * 1)
                elif x == 6:
                    val_x = (-2 * 2) + (-2 * 1)
                
                for a_oxi in el_a.oxidation_states:
                    for b_oxi in el_b.oxidation_states:
                        if (ne_a + ne_b + ne_x) % 2 == 0 and (a_oxi + b_oxi + val_x) == 0:
                            even_found = True
                
                if not even_found:
                    score = score + 100
                cand_score[(a, b, x)] = score
            
    results = sorted(cand_score, key=cand_score.get)
    with open(filename, "wb") as f:
        pickle.dump(results, f)
    return results
    
        
def get_stats(ranked_list):
    # Input is a ranked list
    # each element of the list is a vector of (A, X, B)
    num_good = [0]
    num_cands = [0]
    good_cands = GOOD_CANDS_LS
    
    counter = 0
    for cand in ranked_list:
        counter += 1
        if cand in good_cands:
            num_good.append(num_good[-1] + 1)
            num_cands.append(counter)
    
    return (num_good, num_cands)


def get_stats_OS(ranked_list):
    # Input is a ranked list
    # each element of the list is a vector of (A, X, B)
    num_good = [0]
    num_cands = [0]
    good_cands = GOOD_CANDS_OS
    
    counter = 0
    for cand in ranked_list:
        counter += 1
        if cand in good_cands:
            num_good.append(num_good[-1] + 1)
            num_cands.append(counter)
    
    return (num_good, num_cands)


if __name__ == "__main__":
    a = get_stats_OS(get_ranked_list_goldschmidt_halffill())
    print a