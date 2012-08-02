#!/usr/bin/env python
from __future__ import division
'''
Created on Jul 23, 2012
'''
from ga_optimization_ternary.fitness_evaluators import FitnessEvaluator,\
    eval_fitness_simple
from pymatgen.core.periodic_table import Element
import math
from ga_optimization_ternary.simple_ga import StatTrack
import pickle
import os

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jul 23, 2012"


def get_ranked_list_goldschmidt():
    #TODO: get the oxidation state right!!!
    
    filename = "goldschmidt_rank.p"
    if os.path.exists(filename):
        with open(filename) as f:
            return pickle.load(f)
    
    all_AB = FitnessEvaluator(eval_fitness_simple)._reverse_dict.keys()
    print 'generating goldschmidt ranks...'
    cand_score = {}  # dictionary of cand_tuple:score. a high score is BAD
    for a in all_AB:
        for b in all_AB:
            for x in (True, False):
                r_a = Element.from_Z(a).average_ionic_radius  # TODO: get the correct oxidation state!
                r_b = Element.from_Z(b).average_ionic_radius
                r_x = None
                if x:
                    r_x = Element("O").ionic_radii[-2] * 2/3 + Element("N").ionic_radii[-3] * 1/3
                else:
                    r_x = Element("O").ionic_radii[-2]
                
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

    
        
def get_stats(ranked_list):
    # Input is a ranked list
    # each element of the list is a vector of (A, X, B)
    num_good = [0]
    num_cands = [0]
    good_cands = StatTrack(None)._good_candidates
    
    counter = 0
    for cand in ranked_list:
        counter += 1
        if cand in good_cands:
            num_good.append(num_good[-1] + 1)
            num_cands.append(counter)
    
    return (num_good, num_cands)
        
if __name__ == "__main__":
    
    a, b = get_stats(get_ranked_list_goldschmidt())
    print a
    print b