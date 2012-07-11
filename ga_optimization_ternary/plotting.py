#!/usr/bin/env python
from __future__ import division

'''
Created on Jul 6, 2012
'''
from ga_optimization_ternary.database import M_Database, Stats_Database, MAX_CAND, MAX_GOOD

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jul 6, 2012"

import matplotlib.pyplot as plt
import pymongo
from database import MAX_GOOD
import numpy as np

class PerformancePlot():

    def __init__(self):
        plt.figure(1, figsize=(8,6))
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        num_exps = self.stats_process.count()
        
        self.lw = 4
        self.fontsize = 16
        self.fontname = "Gill Sans"
        
        self.get_reference_data()  # reference
        self.get_data(0, "best GA", "blue")  # best
        self.get_data(num_exps/2, "median GA", "red")  # ~medium
        
        plt.xlabel("Candidates calculated", fontname=self.fontname, fontsize=self.fontsize)
        plt.ylabel("Matching materials found", fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_xticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_yticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.ylim((0, MAX_GOOD + 0.5))
        
    def get_data(self, idx, label, color, crit="ten"):
        x = []
        y = []
        xerr = []
        data = self.stats_process.find({}, sort = [(crit, pymongo.ASCENDING)])[(int)(idx)]
        for idx, val in enumerate(data['ng_avg']):
            y.append(idx)
            x.append(val)
            xerr.append(data['ng_stdev'][idx])
        
        plt.errorbar(x, y, xerr=xerr, lw = self.lw, markersize=9, elinewidth=1, ecolor="black", marker="o", capsize=3, color=color, barsabove=True)
        plt.annotate(label, xy = (x[(int)(MAX_GOOD*3/4)], y[(int)(MAX_GOOD*3/4)]), xytext = (-5, 5), color=color, textcoords = 'offset points', ha = 'right', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    
    def get_reference_data(self):
        x = [0, MAX_CAND]
        y = [0, MAX_GOOD]
        plt.errorbar(x, y, lw=self.lw, color="black")
        plt.annotate("random", xy = (x[1] * 0.66, y[1] * 0.66), xytext = (-15, -30), color="black", textcoords = 'offset points', ha = 'left', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
        

class ComparisonPlot():

    def __init__(self):
        plt.figure(2)
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        self.num_exps = self.stats_process.count()
        
        self.lw = 4
        self.fontsize = 16
        self.fontname = "Gill Sans"
        
        plt.xlabel("GA parameter set rank", fontname=self.fontname, fontsize=self.fontsize)
        plt.ylabel("Efficiency vs. random", fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_xticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_yticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.xlim((0, self.num_exps))
        self.get_reference()
        self.get_data()
        
    
    def get_data(self):
        x = []
        y = []
        idx = 1
        for data in self.stats_process.find({}, sort = [("ten", pymongo.ASCENDING)]):
            x.append(idx)
            idx = idx + 1
            y.append((10/MAX_GOOD * MAX_CAND)/data["ten"])
        plt.errorbar(x, y, lw=self.lw, color="red")
        plt.annotate("10 materials", xy = (x[50], y[50]), xytext = (10, 0), color="red", textcoords = 'offset points', ha = 'left', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
        
        x = []
        y = []
        idx = 1
        for data in self.stats_process.find({}, sort = [("all", pymongo.ASCENDING)]):
            x.append(idx)
            idx = idx + 1
            y.append(MAX_CAND/data["all"])
        plt.errorbar(x, y, lw=self.lw, color="blue")
        plt.annotate("all 15 materials", xy = (x[100], y[100]), xytext = (-25, 10), color="blue", textcoords = 'offset points', ha = 'right', va = 'top', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    
    def get_reference(self):
        x = [0, self.num_exps]
        y = [1, 1]
        plt.plot(x, y, "k--", lw=3)

class ParametersPlot():

    def __init__(self):
        plt.figure(2)
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        self.num_exps = self.stats_process.count()
        
        self.get_data("crossover_fnc")
        
    
    def get_data(self, param):
        plt.title(self.get_pretty_name(param))
        parameter = "crossover_fnc"
        
        labels = ["Baseline", "System"]
        data = [3.75, 4.75]
        
        
        xlocations = np.array(range(len(data))) + 0.5
        width = 0.75
        plt.bar(xlocations, data, width=width)
        plt.xticks(xlocations + width/2, labels)
        plt.xlim(0, xlocations[-1] + width * 2)
    
    def get_pretty_name(self, ugly_name):
        d = {}
        d['crossover_fnc'] = "Crossover Function"
        
        return d[ugly_name]
    
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
        
if __name__ == "__main__":
    ParametersPlot()
    plt.show()

    
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