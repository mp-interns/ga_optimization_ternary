#!/usr/bin/env python
from __future__ import division

'''
Created on Jul 6, 2012
'''
from ga_optimization_ternary.database import M_Database, Stats_Database, MAX_GOOD_LS,\
    NUM_CANDS
from ga_optimization_ternary.utils import get_reference_array
from ga_optimization_ternary import ranked_list_optimization
from ga_optimization_ternary.ranked_list_optimization import get_ranked_list_goldschmidt

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jul 6, 2012"

import matplotlib.pyplot as plt
import pymongo
import numpy as np

def get_pretty_name(ugly_name):
    d = {}
    d['crossover_fnc'] = "Crossover Function"
    d['selection_fnc'] = "Selection Function"
    d['fitness_fnc'] = "Fitness Function"
    d['G1DListCrossoverSinglePoint'] = "Single Point"
    d['G1DListCrossoverTwoPoint'] = "Two Point"
    d['G1DListCrossoverUniform'] = "Uniform"
    d['GRouletteWheel'] = "Roulette"
    d['GTournamentSelector'] = "Tournament"
    d['GUniformSelector'] = "Uniform"
    d['eval_fitness_partial'] = "Partial"
    d['eval_fitness_simple'] = "All-or-Nothing"
    d['initialization_fnc'] = "Initialization Function"
    d['G1DListInitializatorAllele'] = "Random"
    d['popsize'] = "Population size"
    d['elitism_num'] = "Number of Elite"
    d['eval_fitness_complex'] = "Partial"
    return d.get(ugly_name, ugly_name)

class PerformancePlot():

    def __init__(self, format=None):
        plt.figure(1, figsize=(8,6))
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        num_exps = self.stats_process.count()
        
        self.lw = 4
        self.fontsize = 14
        self.fontname = "Trebuchet MS"
        
        self.get_reference_data()  # reference
        self.get_goldschmidt_data()  # goldschmidt reference
        self.get_data(0, "best GA", "blue")  # best
        self.get_data(num_exps-1, "worst GA", "red", "right")  # ~worst
        
        plt.xlabel("Average number of calculations", fontname=self.fontname, fontsize=self.fontsize)
        plt.ylabel("Potential solar light splitting materials", fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_xticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_yticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.ylim((0, MAX_GOOD_LS + 0.5))
        plt.xlim((0, NUM_CANDS))
        
        if format:
            plt.savefig("performance_plot."+format)
        else:
            plt.show()
        
    def get_data(self, idx, label, color, pos="left", crit="fifteen"):
        x = []
        y = []
        xerr = []
        data = self.stats_process.find({}, sort = [(crit, pymongo.ASCENDING)])[(int)(idx)]
        for idx, val in enumerate(data['ng_avg']):
            y.append(idx)
            x.append(val)
            xerr.append(data['ng_stdev'][idx])
        ha, va, xytext = "right", "bottom", (-5, 5)
        
        if pos == "right":
            ha, va, xytext = "left", "bottom", (15, 5)
            
        plt.errorbar(x, y, xerr=xerr, lw = self.lw, markersize=9, elinewidth=1, ecolor="black", marker="o", capsize=3, color=color, barsabove=True)
        plt.annotate(label, xy = (x[(int)(MAX_GOOD_LS*3/4)], y[(int)(MAX_GOOD_LS*3/4)]), xytext = xytext, color=color, textcoords = 'offset points', ha = ha, va = va, fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    
    def get_goldschmidt_data(self):
        color = [.996, .415, 0]
        y, x = ranked_list_optimization.get_stats(get_ranked_list_goldschmidt())
        plt.errorbar(x, y, lw=self.lw, color=color)
        plt.annotate("goldschmidt", xy = (x[15], y[15]), xytext = (5, -15), color=color, textcoords = 'offset points', ha = 'left', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    
    def get_reference_data(self):
        x = [0, get_reference_array()[MAX_GOOD_LS]]
        y = [0, MAX_GOOD_LS]
        plt.errorbar(x, y, lw=self.lw, color="black")
        plt.annotate("random", xy = (x[1] * 0.75, y[1] * 0.75), xytext = (-5, 0), color="black", textcoords = 'offset points', ha = 'right', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
        

class ComparisonPlot():

    def __init__(self, format=None):
        plt.figure(2)
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        self.num_exps = self.stats_process.count()
        
        self.lw = 4
        self.fontsize = 16
        self.fontname = "Trebuchet MS"
        
        plt.xlabel("GA parameter set rank", fontname=self.fontname, fontsize=self.fontsize)
        plt.ylabel("Efficiency vs. random", fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_xticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_yticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.xlim((0, self.num_exps))
        self.get_reference()
        self.get_data() 
        if format:
            plt.savefig("comparison_plot."+format)
        else:
            plt.show()
    
    def get_data(self):
        x = []
        y = []
        idx = 1
        for data in self.stats_process.find({}, sort = [("ten", pymongo.ASCENDING)]):
            x.append(idx)
            idx = idx + 1
            y.append(get_reference_array()[10]/data["ten"])
        plt.errorbar(x, y, lw=self.lw, color="red")
        plt.annotate("10 materials", xy = (x[50], y[50]), xytext = (10, 0), color="red", textcoords = 'offset points', ha = 'left', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
        
        color = [.996, .415, 0]
        x = []
        y = []
        idx = 1
        for data in self.stats_process.find({}, sort = [("fifteen", pymongo.ASCENDING)]):
            x.append(idx)
            idx = idx + 1
            y.append(get_reference_array()[15]/data["fifteen"])
        plt.errorbar(x, y, lw=self.lw, color=color)
        plt.annotate("15 materials", xy = (x[75], y[75]), xytext = (50, -10), color=color, textcoords = 'offset points', ha = 'left', va = 'top', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
        
        x = []
        y = []
        idx = 1
        for data in self.stats_process.find({}, sort = [("all", pymongo.ASCENDING)]):
            x.append(idx)
            idx = idx + 1
            y.append(get_reference_array()[MAX_GOOD_LS]/data["all"])
        plt.errorbar(x, y, lw=self.lw, color="blue")
        plt.annotate("all 22 materials", xy = (x[100], y[100]), xytext = (25, 5), color="blue", textcoords = 'offset points', ha = 'left', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    
    def get_reference(self):
        x = [0, self.num_exps]
        y = [1, 1]
        plt.plot(x, y, "k--", lw=3)


class ParametersPlot():

    def __init__(self, format=None):
        plt.figure(1, figsize=(16,8))
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        self.num_exps = self.stats_process.count()
        
        self.fontname = "Trebuchet MS"
        
        params = ["crossover_fnc", "popsize", "selection_fnc", "tournament_rate", "mutation_rate", "elitism_num", "fitness_fnc", "fitness_temp", "initialization_fnc"]
        
        plt.subplot(3, 3, 1)
        for idx, param in enumerate(params):
            plt.subplot(3, 3, idx+1)
            self.get_data(param)
        
        if format:
            plt.savefig("parameters_plot."+format)
        else:
            plt.show()
    
    def get_data(self, parameter):
        plt.title(get_pretty_name(parameter), fontname=self.fontname)
        
        data = []
        labels = self.stats_process.distinct("parameters."+parameter)
        labels.sort()
        for label in labels:
            # get the best
            best = self.stats_process.find({"parameters."+parameter:label}, sort=[("all",pymongo.ASCENDING)])[0]
            data.append(get_reference_array()[MAX_GOOD_LS]/best["all"])
        
        pretty_labels = [get_pretty_name(n) for n in labels]
        xlocations = np.array(range(len(data))) + 0.5
        width = 0.75
        plt.bar(xlocations, data, width=width)
        plt.xticks(xlocations + width/2, pretty_labels, fontname=self.fontname)
        plt.yticks(fontname=self.fontname)
        plt.xlim(0, xlocations[-1] + width * 2)
    
class TenAllPlot():

    def __init__(self):
        plt.figure(3)
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        num_exps = self.stats_process.count()
        
        self.get_data()
        
    
    def get_data(self):
        x = []
        y = []
        for data in self.stats_process.find({}, sort = [("ten", pymongo.ASCENDING)]):
            x.append(data['ten'])
            y.append(data['all'])
        plt.scatter(x, y)
    
class DataTable():
    
    def __init__(self):
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        
        for it in self.stats_process.find():
            p = it['parameters']
            
            print ("{}\t{}\t{}\t{}\t{}\t{}\t{}").format(p['popsize'], get_pretty_name(p['selection_fnc']), get_pretty_name(p['fitness_fnc']),get_pretty_name(p['crossover_fnc']), p['elitism_num'], it['ten'], it['all'])
    
if __name__ == "__main__":
    #PerformancePlot()
    #ComparisonPlot()
    ParametersPlot()
    #plt.show()
    
    # DataTable()
    
"""
# example data
x = [1, 2, 3, 4, 5]
y = [1, 2, 3, 4, 5]
x2 = [2, 4, 6, 8, 10]
xerr = [1, 2, 3, 4, 5]

# First illustrate basic pyplot interface, using defaults where possible.
plt.xlim((0,6))
plt.ylim((0,6))
plt.errorbar(x, y, xerr=xerr)
plt.errorbar(x2, y, xerr=xerr)

print 'yay'
plt.show()
"""