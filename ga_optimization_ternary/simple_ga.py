#!/usr/bin/env python

from __future__ import division
'''
Created on Mar 14, 2012
'''
from ga_optimization_ternary.fitness_evaluators import eval_fitness_simple, eval_fitness_complex
from test.test_descrtut import defaultdict
from collections import OrderedDict
from ga_optimization_ternary.database import MAX_GOOD_LS, GOOD_CANDS_LS

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Mar 14, 2012"

from pyevolve import G1DList, Mutators, Initializators, GAllele
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Statistics, Crossovers
from fitness_evaluators import FitnessEvaluator
from database import Stats_Database

from scipy.interpolate import interp1d

import multiprocessing
import math

class AllFound():
    ALL_FOUND = False
                

def AllFoundCriteria(ga_engine):
    """ Terminate the evolution based on the raw stats
    """
    return AllFound.ALL_FOUND

class ParameterSet():
    
    def __init__(self, crossover_fnc, fitness_fnc, fitness_temp, selection_fnc, tournament_rate, mutation_fnc, mutation_rate, initialization_fnc, popsize, elitism_num, niching_bool):
        self.crossover_fnc = crossover_fnc
        self.fitness_fnc = fitness_fnc
        self.fitness_temp = fitness_temp
        self.popsize = popsize
        self.selection_fnc = selection_fnc
        self.mutation_fnc = mutation_fnc
        self.mutation_rate = mutation_rate
        self.tournament_rate = tournament_rate
        self.initialization_fnc = initialization_fnc
        self.elitism_num = elitism_num
        self.niching_bool = niching_bool
    
    def to_dict(self):
        d = OrderedDict()
        d['crossover_fnc'] = self.crossover_fnc.__name__
        d['fitness_fnc'] = self.fitness_fnc.__name__
        d['fitness_temp'] = self.fitness_temp
        d['selection_fnc'] = self.selection_fnc.__name__
        d['mutation_fnc'] = self.mutation_fnc.__name__
        d['mutation_rate'] = self.mutation_rate
        d['tournament_rate'] = self.tournament_rate
        d['initialization_fnc'] = self.initialization_fnc.__name__
        d['popsize'] = self.popsize
        d['elitism_num'] = self.elitism_num
        d['niching'] = self.niching_bool
        return d
    
    def unique_key(self):
        return '-'.join([str(val) for val in self.to_dict().values()])


class StatTrack():
    
    def __init__(self, fe, mutation_rate, tournament_rate):
        self._candidates_tried = set()
        self._candidates_good = set()
        self.generation_ncandidates = [0]
        self.generation_ngood = [0]
        self._fitness_evaluator = fe
        self.num_breakouts = 0
        self.mutation_rate = mutation_rate
        self.tournament_rate = tournament_rate
    
    def updateStats(self, generation_num, population):
        for i in population:
            cand = self._fitness_evaluator.convert_raw_to_Z((i[0], i[1], i[2]))
            self._candidates_tried.add(cand)
            
            if cand in GOOD_CANDS_LS:
                self._candidates_good.add(cand)
            
        self.generation_ncandidates.append(len(self._candidates_tried))
        self.generation_ngood.append(len(self._candidates_good))
        if len(self._candidates_good) == MAX_GOOD_LS:
            AllFound.ALL_FOUND = True
            
        return self.generation_ncandidates[-1] - self.generation_ncandidates[-2]
        
    def evolve_callback(self, ga):
        ga.getPopulation().setParams(tournamentPool=(int)(math.ceil(self.tournament_rate * len(ga.getPopulation()))))
        cands_added = self.updateStats(ga.currentGeneration, ga.getPopulation().internalPop)
        
        breakout_cutoff = (int)(math.ceil(0.2 * len(ga.getPopulation().internalPop)))
        #breakout_cutoff = 10
        if cands_added < breakout_cutoff:
            ga.setMutationRate(0.75)
        else:
            ga.setMutationRate(self.mutation_rate)
            self.num_breakouts += 1
            
        #if ga.currentGeneration == 0:
            '''
            print 'YAY'
            pop = ga.getPopulation()
            for p in pop:
                print p.__class__
            # CHANGE THE ENTIRE POPULATION "pop"
            #pop[i].genomeList = [random.random(), 0 for i in xrange(chromosome/2)]
            #newPop.evaluate()
            #newPop.sort()
            '''
        return False
        
    
    def get_interpolation_function(self):
        return interp1d(self.generation_ncandidates, self.generation_ngood)


def run_simulation(pset, max_generations):
    # Genome instance
    setOfAlleles = GAllele.GAlleles()
    setOfAlleles.add(GAllele.GAlleleRange(0, 51))
    setOfAlleles.add(GAllele.GAlleleRange(0, 6))
    setOfAlleles.add(GAllele.GAlleleRange(0, 51))

    # Genome instance, 1D List of 50 elements
    genome = G1DList.G1DList(3)
    genome.setParams(allele=setOfAlleles)
    
    #Fitness function
    
    fe = FitnessEvaluator(pset.fitness_fnc, pset.fitness_temp)
    
    st = StatTrack(fe, pset.mutation_rate, pset.tournament_rate)
    
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
        ga.setElitismReplacement((int)(math.ceil(pset.elitism_num * pset.popsize)))
    else:
        ga.setElitism(False)
        
    ga.setMutationRate(pset.mutation_rate)
    
    # TODO: figure out niching
    ga.evolve()
    
    return st
    

def main_loop():
    ncores = 7
    clear = False
    # clear the Stats DB
    db = Stats_Database(clear=clear)
    popsizes = [16, 100, 500, 1000]  #TODO: add 250 and 2000 later, if it is needed based on the data
    fitness_fncs = [eval_fitness_simple, eval_fitness_complex]
    fitness_temps = [5, 10, 30]
    crossover_fncs = [Crossovers.G1DListCrossoverUniform, Crossovers.G1DListCrossoverSinglePoint, Crossovers.G1DListCrossoverTwoPoint]
    selection_fncs = [Selectors.GTournamentSelectorAlternative, Selectors.GRankSelector, Selectors.GRouletteWheel, Selectors.GUniformSelector]
    mutator_fncs = [Mutators.G1DListMutatorAllele]
    tournament_rates = [0.05, 0.1, 0.5]
    mutation_rates = [0.01, 0.1, 0.25]
    elitisms = [0, 0.1, 0.5]
    nichings = [False]  # TODO: implement True
    initialization_fncs = [Initializators.G1DListInitializatorAllele]  # TODO: add data-mined initializors
    
    all_ps = []
    for popsize in popsizes:
        for fitness_fnc in fitness_fncs:
            for crossover_fnc in crossover_fncs:
                for selection_fnc in selection_fncs:
                    for mutator_fnc in mutator_fncs:
                        for elitism in elitisms:
                            for niching in nichings:
                                for initialization_fnc in initialization_fncs:
                                    for m_rate in mutation_rates:
                                        for fitness_temp in fitness_temps:
                                            for tournament_rate in tournament_rates:
                                                if (tournament_rate == tournament_rates[0] or selection_fnc == selection_fncs[0] or selection_fnc == selection_fncs[1]):
                                                    all_ps.append(ParameterSet(crossover_fnc, fitness_fnc, fitness_temp, selection_fnc, tournament_rate, mutator_fnc, m_rate, initialization_fnc, popsize, elitism, niching))  # set up the parameters

    print 'the number of parameter sets is', len(all_ps)
    process_parallel(all_ps, ncores)
    # process_serial(all_ps)

    
def process_parameterset(ps):
    num_iterations = 15
    production = True
    max_generations = 20000  # should always work...(hopefully)
    
    if production:
        db = Stats_Database(clear=False)
        if not db._stats_raw.find({"unique_key": ps.unique_key()}).count() >= num_iterations:
            db._stats_raw.remove({"unique_key": ps.unique_key()}, safe=True)
            for iteration in range(num_iterations):
                AllFound.ALL_FOUND = False  # reset the simulation
                stat = run_simulation(ps, max_generations)
                print stat.generation_ncandidates[-1], len(stat.generation_ngood), stat.num_breakouts, ps.to_dict(), ps.unique_key()
                if production:
                    db.add_stat_raw(ps, stat, iteration)
        else:
            print 'SKIPPING KEY' + ps.unique_key()
        # stats.append(stat)
    
    #if production:
    #    db.add_stats_raw(ps, stats)
    
    return 0


def process_parallel(all_ps, ncores):
    pool = multiprocessing.Pool(ncores)
    states = pool.map(process_parameterset, all_ps)
    state = 0 if all([s == 0 for s in states]) else -1
    print "FINISHED with state", state


def process_serial(all_ps):
    for ps in all_ps:
        process_parameterset(ps)

if __name__ == "__main__":
    main_loop()
