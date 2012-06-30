#!/usr/bin/env python

from __future__ import division
'''
Created on Mar 14, 2012
'''
from ga_optimization_ternary.fitness_evaluators import eval_fitness_simple,\
    eval_fitness_complex, eval_fitness_partial
from test.test_descrtut import defaultdict
from collections import OrderedDict

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Mar 14, 2012"

from pyevolve import G1DList, Mutators, Initializators
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Statistics, Crossovers
from fitness_evaluators import FitnessEvaluatorZ
from database import Stats_Database

from scipy.interpolate import interp1d

class AllFound():
    ALL_FOUND = False
                

def AllFoundCriteria(ga_engine):
    """ Terminate the evolution based on the raw stats
    """
    return AllFound.ALL_FOUND

class ParameterSet():
    
    def __init__(self, crossover_fnc, fitness_fnc, selection_fnc, mutation_fnc, initialization_fnc, popsize, elitism_num, niching_bool):
        self.crossover_fnc = crossover_fnc
        self.fitness_fnc = fitness_fnc
        self.popsize = popsize
        self.selection_fnc = selection_fnc
        self.mutation_fnc = mutation_fnc
        self.initialization_fnc = initialization_fnc
        self.elitism_num = elitism_num
        self.niching_bool = niching_bool
    
    def to_dict(self):
        d = OrderedDict()
        d['crossover_fnc'] = self.crossover_fnc.__name__
        d['fitness_fnc'] = self.fitness_fnc.__name__
        d['selection_fnc'] = self.selection_fnc.__name__
        d['mutation_fnc'] = self.mutation_fnc.__name__
        d['initialization_fnc'] = self.initialization_fnc.__name__
        d['popsize'] = self.popsize
        d['elitism_num'] = self.elitism_num
        d['niching'] = self.niching_bool
        return d
    
    def unique_key(self):
        return '-'.join([str(val) for val in self.to_dict().values()])


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


def run_simulation(pset, max_generations):
    
    # Genome instance, 1D List of 50 elements
    genome = G1DList.G1DList(3)
    # Sets the range max and min of the 1D List
    genome.setParams(rangemin=0, rangemax=51)
    
    #Fitness function
    fe = FitnessEvaluatorZ(pset.fitness_fnc)
    st = StatTrack(fe)
    
    genome.crossover.set(pset.crossover_fnc)
    genome.mutator.set(pset.mutation_fnc)
    genome.evaluator.set(fe.array_to_score)
    genome.initializator.set(pset.initialization_fnc)

    ga = GSimpleGA.GSimpleGA(genome)

    ga.selector.set(pset.selection_fnc)
    ga.terminationCriteria.set(AllFoundCriteria)
    ga.stepCallback.set(st.evolve_callback)
    ga.setPopulationSize(pset.popsize)
    ga.setGenerations(max_generations)
    if pset.elitism_num > 0:
        ga.setElitism(True)
        ga.setElitismReplacement(pset.elitism_num)
    else:
        ga.setElitism(False)
    
    # TODO: figure out niching
    
    ga.evolve()
    
    return st
    

def main_loop():
    max_generations = 10000  # should always work...(hopefully)
    num_iterations = 50
    db = Stats_Database(clear=True)
    popsizes = [125, 250, 500, 1000]
    fitness_fncs = [eval_fitness_simple, eval_fitness_partial]
    crossover_fncs = [Crossovers.G1DListCrossoverUniform, Crossovers.G1DListCrossoverSinglePoint, Crossovers.G1DListCrossoverTwoPoint]  # TODO: is Single point different than 2 point?
    selection_fncs = [Selectors.GTournamentSelector, Selectors.GRouletteWheel, Selectors.GUniformSelector]  # Rank selector is SLOW...
    mutator_fncs = [Mutators.G1DListMutatorIntegerRange]
    elitisms = [0, 1, 2, 5]
    nichings = [False]  # TODO: implement True
    initialization_fncs = [Initializators.G1DListInitializatorInteger]  # TODO: add data-mined initializors
    
    for popsize in popsizes:
        for fitness_fnc in fitness_fncs:
            for crossover_fnc in crossover_fncs:
                for selection_fnc in selection_fncs:
                    for mutator_fnc in mutator_fncs:
                        for elitism in elitisms:
                            for niching in nichings:
                                for initialization_fnc in initialization_fncs:
                                    stats = []
                                    for iteration in range(num_iterations):
                                        AllFound.ALL_FOUND = False  # reset the simulation
                                        ps = ParameterSet(crossover_fnc, fitness_fnc, selection_fnc, mutator_fnc, initialization_fnc, popsize, elitism, niching)  # set up the parameters
                                        stat = run_simulation(ps, max_generations)
                                        print stat.generation_ncandidates[-1], ps.to_dict(), ps.unique_key()
                                        stats.append(stat)
                                    db.add_stats_raw(ps, stats)
                                    
                                    
    
    
if __name__ == "__main__":
    main_loop()
