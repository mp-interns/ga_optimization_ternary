#!/usr/bin/env python

from __future__ import division
'''
Created on Mar 14, 2012
'''
from pyevolve.Mutators import G1DListMutatorIntegerRange

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Mar 14, 2012"

from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Statistics, Crossovers
from pyevolve.DBAdapters import DBFileCSV
from fitness_evaluators import FitnessEvaluatorZ

from scipy.interpolate import interp1d
import numpy as np

# This function is the evaluation function, we want
# to give high score to more zero'ed chromosomes
ALL_FOUND = False

class StatTrack():
    
    def __init__(self, fe):
        self._candidates_tried = set()
        self._candidates_good = set()
        #this is in Z representation
        self._good_candidates = [(3, 23, False), (11, 51, False), (20, 32, False), (20, 50, False), (38, 32, False), (38, 50, False), (47, 41, False), (50, 22, False), (55, 41, False), (56, 50, False)]
        self._good_candidates.extend([(12, 73, True), (20, 73, True), (38, 73, True), (56, 73, True), (57, 22, True)])
        self.generation_ncandidates = [0]
        self.generation_ngood = [0]
        self._fitness_evaluator = fe
    
    def updateStats(self, generation_num, population):
        for i in population:
            cand = self._fitness_evaluator.convert_raw_to_Z((i[0], i[1], i[2]))
            self._candidates_tried.add(cand)
            
            if cand in self._good_candidates:
                self._candidates_good.add(cand)
            
        self.generation_ncandidates.append(len(self._candidates_tried))
        self.generation_ngood.append(len(self._candidates_good))
        if len(self._candidates_good) == len(self._good_candidates):
            AllFound.ALL_FOUND = True
            
        return self.generation_ncandidates[-1] - self.generation_ncandidates[-2]
        
    def evolve_callback(self, ga):
        cands_added = self.updateStats(ga.currentGeneration, ga.getPopulation().internalPop)
        if cands_added < 10:
            ga.setMutationRate(0.75)
        else:
            ga.setMutationRate(0.02)
    
    def get_interpolation_function(self):
        return interp1d(self.generation_ncandidates, self.generation_ngood)
    
def run_main(fe, popsize, gens):
    
    #Fitness function
    st = StatTrack(fe)
   
    # Genome instance, 1D List of 50 elements
    genome = G1DList.G1DList(3)

    # Sets the range max and min of the 1D List
    genome.setParams(rangemin=0, rangemax=51)
    genome.crossover.set(Crossovers.G1DListCrossoverUniform)
    genome.mutator.set(G1DListMutatorIntegerRange)

    # The evaluator function (evaluation function)
    genome.evaluator.set(fe.array_to_score)
   
    # Genetic Algorithm Instance
    ga = GSimpleGA.GSimpleGA(genome)

    # Set the Roulette Wheel selector method, the number of generations and
    # the termination criteria
    ga.selector.set(Selectors.GTournamentSelector)
    ga.terminationCriteria.set(AllFoundCriteria)
    ga.stepCallback.set(st.evolve_callback)
    ga.setPopulationSize(popsize)
    ga.setGenerations(gens)
    ga.evolve()
    
    return st
        
def iteration_test(fe, iterations = 20, popsize = 500, gens = 500):
    statslist = []
    
    min_candidates = 10000000
    
    for i in range(iterations):
        AllFound.ALL_FOUND=False
        stat=run_main(fe, popsize, gens)
        
        my_hit_rate = stat.generation_ngood[-1]/float(stat.generation_ncandidates[-1])
        usual_hit_rate = len(stat._good_candidates)/5408  # 5408 is total number of calcs in DB
        print stat.generation_ncandidates[-1], stat.generation_ngood[-1], (my_hit_rate)/(usual_hit_rate)
        statslist.append(stat.get_interpolation_function())
        if stat.generation_ncandidates[-1] < min_candidates:
            min_candidates = stat.generation_ncandidates[-1]
        
    print '----'
    
    for i in range(min_candidates):
        values = [stats(i) for stats in statslist]
        print i,'\t', float(sum(values))/len(values)


class AllFound():
    ALL_FOUND = False
                
def AllFoundCriteria(ga_engine):
    """ Terminate the evolution based on the raw stats
    """
    return AllFound.ALL_FOUND
    
    
if __name__ == "__main__":
    iteration_test(FitnessEvaluatorZ())
    
