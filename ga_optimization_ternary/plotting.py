#!/usr/bin/env python
from __future__ import division

'''
Created on Jul 6, 2012
'''
from ga_optimization_ternary.database import M_Database, Stats_Database, NUM_CANDS,\
    GOOD_CANDS_LS, GOOD_CANDS_OS
from ga_optimization_ternary.utils import get_reference_array,\
    get_reference_array_OS
from ga_optimization_ternary import ranked_list_optimization
from ga_optimization_ternary.ranked_list_optimization import get_ranked_list_goldschmidt_halffill
from twisted.test.test_text import lineWidth

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jul 6, 2012"

import matplotlib as mpl
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
    d['GTournamentSelector'] = "T"
    d['GRankSelector'] = "Rank"
    d['GTournamentSelectorAlternative'] = "Tourn."
    d['GUniformSelector'] = "Uniform"
    d['eval_fitness_partial'] = "Smooth"
    d['eval_fitness_simple'] = "Discontinuous"
    d['initialization_fnc'] = "Initialization Function"
    d['G1DListInitializatorAllele'] = "Random"
    d['popsize'] = "Population size"
    d['mutation_rate'] = "Mutation Rate"
    d['elitism_num'] = "Number of Elite"
    d['eval_fitness_complex'] = "Smooth"
    d["GTournamentSelectorAlternative-0.25-1.25"] = "T-0.25"
    d["GTournamentSelectorAlternative-0.05-1.25"] = "T-0.05"
    d["GTournamentSelectorAlternative-0.1-1.25"] = "T-0.10"
    d["GTournamentSelectorAlternative-2-1.25"] = "T-Bin"
    d["GRouletteWheel-0.05-1.25"] = "R-1.25"
    d["GRouletteWheel-0.05-2.5"] = "R-2.5"
    d["GRouletteWheel-0.05-5"] = "R-5"
    d["GRouletteWheel-0.05-10"] = "R-10"
    d["GRouletteWheel-0.05-20"] = "R-20"
    d["GUniformSelector-0.05-1.25"] = "U-0"
    return d.get(ugly_name, ugly_name)


def get_short_name(ugly_name):
    d = {}
    d['crossover_fnc'] = "X"
    d['mutation_rate'] = "M"
    d['elitism_num'] = "E"
    d['fitness_fnc'] = "F"
    d['selection_overall'] = "S"
    d['popsize'] = "P"
    
    d['G1DListCrossoverSinglePoint'] = "1P"
    d['G1DListCrossoverTwoPoint'] = "2P"
    d['G1DListCrossoverUniform'] = "U"
    d['eval_fitness_partial'] = "S"
    d['eval_fitness_simple'] = "D"
    d['eval_fitness_complex'] = "S"
    d["GTournamentSelectorAlternative-0.25-1.25"] = "T-0.25"
    d["GTournamentSelectorAlternative-0.05-1.25"] = "T-0.05"
    d["GTournamentSelectorAlternative-0.1-1.25"] = "T-0.10"
    d["GTournamentSelectorAlternative-2-1.25"] = "T-Bin"
    d["GRouletteWheel-0.05-1.25"] = "R-1.25"
    d["GRouletteWheel-0.05-2.5"] = "R-2.5"
    d["GRouletteWheel-0.05-5"] = "R-5"
    d["GRouletteWheel-0.05-10"] = "R-10"
    d["GRouletteWheel-0.05-20"] = "R-20"
    d["GUniformSelector-0.05-1.25"] = "U-0"
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
        # self.get_goldschmidt_data()  # goldschmidt reference
        self.get_data(0, "best GA", "dodgerblue", pos="right", crit="mixed")  # best
        # self.get_data(0, "best GA (ten)", "green", pos="right", crit="ten")  # best
        self.get_data(num_exps-1, "worst GA", "tomato", "right")  # ~worst
        
        plt.xlabel("Average number of calculations", fontname=self.fontname, fontsize=self.fontsize)
        plt.ylabel("Potential solar light splitting materials", fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_xticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_yticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.ylim((0, len(GOOD_CANDS_LS) + 0.5))
        plt.xlim((0, NUM_CANDS))
        
        if format:
            plt.savefig("performance_plot."+format)
        else:
            plt.show()
        
    def get_data(self, idx, label, color, pos="left", crit="all"):
        x = []
        y = []
        xerr = []
        if crit == "mixed":
            data = self.stats_process.find({"ten":{"$lte":1060}, "all":{"$lte":4507}}, sort = [(crit, pymongo.ASCENDING)])[(int)(idx)]
        else:
            data = self.stats_process.find({}, sort = [(crit, pymongo.ASCENDING)])[(int)(idx)]
        for idx, val in enumerate(data['ng_avg']):
            y.append(idx)
            x.append(val)
            xerr.append(data['ng_stdev'][idx])
        ha, va, xytext = "right", "bottom", (-5, 5)
        
        if pos == "right":
            ha, va, xytext = "left", "bottom", (25, 5)
            
        plt.errorbar(x, y, xerr=xerr, lw = self.lw, markersize=9, elinewidth=1, ecolor="black", marker="o", capsize=3, color=color, barsabove=True)
        plt.annotate(label, xy = (x[(int)(len(GOOD_CANDS_LS)*3/4)], y[(int)(len(GOOD_CANDS_LS)*3/4)]), xytext = xytext, color=color, textcoords = 'offset points', ha = ha, va = va, fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    
    def get_goldschmidt_data(self):
        color = [.996, .415, 0]
        y, x = ranked_list_optimization.get_stats(get_ranked_list_goldschmidt_halffill())
        plt.errorbar(x, y, lw=self.lw, color=color)
        plt.annotate("chemical\nrules", xy = (x[19], y[19]), xytext = (-15, 15), color=color, textcoords = 'offset points', ha = 'right', va = 'top', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    
    def get_reference_data(self):
        x = [0, get_reference_array()[len(GOOD_CANDS_LS)]]
        y = [0, len(GOOD_CANDS_LS)]
        plt.errorbar(x, y, lw=self.lw, color="black")
        plt.annotate("random", xy = (x[1] * 0.95, y[1] * 0.95), xytext = (-5, 5), color="black", textcoords = 'offset points', ha = 'right', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)

class PerformancePlotExclusion():

    def __init__(self, format=None):
        plt.figure(1, figsize=(8,6))
        self.db = Stats_Database()
        self.db2 = Stats_Database(extension="_exclusion")
        self.stats_process = self.db._stats_process
        self.stats_process2 = self.db2._stats_process
        
        self.lw = 4
        self.fontsize = 14
        self.fontname = "Trebuchet MS"
        
        self.get_goldschmidt_data()  # goldschmidt reference
        self.get_data(0, "best GA", "dodgerblue", pos="right", crit="all")  # best
        self.get_data_exclusion(0, "best GA + chemical", "#00B31B", pos="left", crit="all")  # best
        
        plt.xlabel("Average number of calculations", fontname=self.fontname, fontsize=self.fontsize)
        plt.ylabel("Potential solar light splitting materials", fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_xticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_yticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.ylim((0, len(GOOD_CANDS_LS) + 0.5))
        plt.xlim((0, 5500))
        
        if format:
            plt.savefig("performance_plot_exclusion."+format)
        else:
            plt.show()
        
    def get_data(self, idx, label, color, pos="left", crit="all"):
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
            ha, va, xytext = "left", "bottom", (25, 5)
            
        plt.errorbar(x, y, xerr=xerr, lw = self.lw, markersize=9, elinewidth=1, ecolor="black", marker="o", capsize=3, color=color, barsabove=True)
        plt.annotate(label, xy = (x[(int)(len(GOOD_CANDS_LS)*3/4)], y[(int)(len(GOOD_CANDS_LS)*3/4)]), xytext = xytext, color=color, textcoords = 'offset points', ha = ha, va = va, fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    
    def get_data_exclusion(self, idx, label, color, pos="left", crit="all"):
        x = []
        y = []
        xerr = []
        data = self.stats_process2.find({}, sort = [(crit, pymongo.ASCENDING)])[(int)(idx)]
        for idx, val in enumerate(data['ng_avg']):
            y.append(idx)
            x.append(val)
            xerr.append(data['ng_stdev'][idx])
        ha, va, xytext = "right", "bottom", (-10, 5)
        
        if pos == "right":
            ha, va, xytext = "left", "bottom", (25, 5)
            
        plt.errorbar(x, y, xerr=xerr, lw = self.lw, markersize=9, elinewidth=1, ecolor="black", marker="o", capsize=3, color=color, barsabove=True)
        plt.annotate(label, xy = (x[(int)(len(GOOD_CANDS_LS)*19/20)], y[(int)(len(GOOD_CANDS_LS)*19/20)]), xytext = xytext, color=color, textcoords = 'offset points', ha = ha, va = va, fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
        
    def get_goldschmidt_data(self):
        color = "tomato"
        y, x = ranked_list_optimization.get_stats(get_ranked_list_goldschmidt_halffill())
        plt.errorbar(x, y, lw=self.lw, color=color)
        plt.annotate("chemical", xy = (x[19], y[19]), xytext = (10, 15), color=color, textcoords = 'offset points', ha = 'right', va = 'top', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    
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
        plt.errorbar(x, y, linestyle="--", lw=self.lw, color="blue")
        plt.annotate("10 materials", xy = (x[50], y[50]), xytext = (10, 0), color="blue", textcoords = 'offset points', ha = 'left', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
        
        x = []
        y = []
        idx = 1
        for data in self.stats_process.find({}, sort = [("all", pymongo.ASCENDING)]):
            x.append(idx)
            idx = idx + 1
            y.append(get_reference_array()[len(GOOD_CANDS_LS)]/data["all"])
        plt.errorbar(x, y, lw=self.lw, color="blue")
        plt.annotate("all 20 materials", xy = (x[100], y[100]), xytext = (25, 5), color="blue", textcoords = 'offset points', ha = 'left', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    
    def get_reference(self):
        x = [0, self.num_exps]
        y = [1, 1]
        plt.plot(x, y, "k-", lw=3)
        plt.annotate("no gain", xy = (x[0], y[0]), xytext = (25, 5), color="black", textcoords = 'offset points', ha = 'left', va = 'bottom', fontname=self.fontname, fontsize=self.fontsize, arrowprops = None)
    

class ParametersPlot():

    def __init__(self, format=None, criteria="all"):
        plt.figure(1, figsize=(16,8))
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        self.num_exps = self.stats_process.count()
        self.criteria = criteria
        self.fontname = "Trebuchet MS"
        
        params = ["crossover_fnc", "popsize", "selection_overall", "mutation_rate", "elitism_num", "fitness_fnc"]
        
        plt.subplot(3, 2, 1)
        for idx, param in enumerate(params):
            plt.subplot(3, 2, idx+1)
            self.get_data(param)
        
        if format:
            mpl.rcParams['figure.figsize'] = (16, 8)
            print mpl.rcParams['figure.figsize']
            
            plt.savefig("parameters_plot."+format)
        else:
            plt.show()
    
    def get_data(self, parameter):
        plt.title(get_pretty_name(parameter), fontname=self.fontname)
        
        data = []
        labels = self.stats_process.distinct("parameters."+parameter)
        if "-" not in str(labels[0]):
            labels.sort()
        else:
            labels.sort(key=lambda label: ord(label.split("-")[0][1])*10 + float(label.split("-")[1])+float(label.split("-")[2]))
        for label in labels:
            # get the best
            best = self.stats_process.find({"parameters."+parameter:label}, sort=[(self.criteria, pymongo.ASCENDING)])[0]
            if self.criteria == "all":
                data.append(get_reference_array()[len(GOOD_CANDS_LS)]/best["all"])
            elif self.criteria == "ten":
                data.append(get_reference_array()[10]/best["ten"])
            '''
            # get the average
            m_sum = 0.0
            npoints = 0.0
            num_to_average = 1  # average the num_to_average BEST runs
            for item in self.stats_process.find({"parameters." + parameter: label}, {"all": 1}, sort=[("all",pymongo.ASCENDING)]).limit(num_to_average):
                m_sum = m_sum + get_reference_array()[len(GOOD_CANDS_LS)] / item["all"]
                npoints = npoints + 1
            data.append(m_sum / npoints)
            '''
        pretty_labels = [get_pretty_name(n) for n in labels]
        xlocations = np.array(range(len(data))) + 0.5
        width = 0.75
        plt.bar(xlocations, data, width=width)
        plt.xticks(xlocations + width/2, pretty_labels, fontname=self.fontname)
        plt.yticks(fontname=self.fontname)
        plt.xlim(0, xlocations[-1] + width * 2)



class HeatMapPlot():

    def __init__(self, format=None, criteria="all"):
        plt.figure(1, figsize=(12,8))
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        self.num_exps = self.stats_process.count()
        self.criteria = criteria
        self.fontname = "Trebuchet MS"
        self.fontsize = 14
        
        params = ["crossover_fnc", "popsize", "selection_overall", "elitism_num", "fitness_fnc", "mutation_rate"]
        param_vals = []
        
        
        for p1 in params:
            p1_labels = self.stats_process.distinct("parameters."+p1)
            if "-" not in str(p1_labels[0]):
                p1_labels.sort()
            else:
                p1_labels.sort(key=lambda label: ord(label.split("-")[0][1])*10 + float(label.split("-")[1])+float(label.split("-")[2]))
            
            for p1x in p1_labels:
                param_vals.append((p1,p1x))
        
        values = []
        
        for x1 in param_vals:
            for y1 in param_vals:
                if x1[0] == y1[0] and x1[1] != y1[1]:
                    values.append(0.5)
                elif x1[0] == y1[0]:
                    constraint = {"parameters."+x1[0]:x1[1]}
                    values.append(self.get_score(constraint))
                else:
                    constraint = {"parameters."+x1[0]:x1[1], "parameters."+y1[0]:y1[1]}
                    values.append(self.get_score(constraint))
                
        values = np.array(values)
        values.shape = (len(param_vals),len(param_vals))
              
        plt.hot()
        plt.pcolormesh(values)
        x = [str(get_short_name(x[1]))+"="+get_short_name(x[0]) for x in param_vals]
        ticklok = [z+0.5 for z in range(0, len(param_vals))]
        plt.xticks(ticklok, x, ha='center', rotation=90, fontsize=14, fontname=self.fontname)
        plt.yticks(ticklok, x, ha='right', fontsize=14, fontname=self.fontname)
        
        last_name = param_vals[0][0]
        
        for idx, label in enumerate(param_vals):
            if label[0] != last_name:
                plt.plot([idx,idx], [0, len(param_vals)], color="darkslategray", linewidth=3)
                plt.plot([0, len(param_vals)], [idx,idx],  color="darkslategray", linewidth=3)
                last_name = label[0]
        
        plt.xlim(0, len(param_vals))
        plt.ylim(0, len(param_vals))
        plt.colorbar() 
            
        if format:
            plt.savefig("heatmap_plot."+format)
        else:
            plt.show()
    
    def get_score(self, constraint):
        m_sum = 0.0
        npoints = 0.0
        num_to_average = 50000  # average the num_to_average BEST runs
        for item in self.stats_process.find(constraint, {"all": 1}, sort=[("all",pymongo.ASCENDING)]).limit(num_to_average):
            m_sum = m_sum + get_reference_array()[len(GOOD_CANDS_LS)] / item["all"]
            npoints = npoints + 1
        return (m_sum / npoints)
    


class BreakoutPlot():

    def __init__(self, format=None):
        plt.figure(3)
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        self.fontname = "Trebuchet MS"
        self.fontsize = 14
        
        self.get_data()
        
        plt.xlabel("Average number of global mutations", fontname=self.fontname, fontsize=self.fontsize)
        plt.ylabel("Efficiency (all 20 materials)", fontname=self.fontname, fontsize=self.fontsize)
        
        if format:
            plt.savefig("breakout_plot."+format)
        else:
            plt.show()
    
    def get_data(self):
        
        colors = {100: "tomato", 500: "dodgerblue", 1000: "springgreen"}
        markers = {100: "s", 500: "o", 1000: "D"}

        for popsize in colors:
            x = []
            y = []
            for data in self.stats_process.find({"parameters.popsize":popsize}):
                x.append(data['breakouts_avg'])
                y.append(get_reference_array()[len(GOOD_CANDS_LS)]/data['all'])
            plt.scatter(x, y, c=[colors[popsize]]*len(x), marker=markers[popsize], s=30, label="popsize="+str(popsize))
        
        plt.legend(loc=0)
        yxis = plt.ylim()
        plt.ylim(0.75, yxis[1])
        plt.yticks(np.arange(1,plt.ylim()[1], 1))
        
        xxis = plt.xlim()
        plt.xlim(0, xxis[1])


class TenAllPlot():

    def __init__(self, format=None):
        plt.figure(3)
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        self.fontname = "Trebuchet MS"
        self.fontsize = 14
        
        self.get_data()
        
        plt.xlabel("Efficiency (10 materials)", fontname=self.fontname, fontsize=self.fontsize)
        plt.ylabel("Efficiency (all 20 materials)", fontname=self.fontname, fontsize=self.fontsize)
        
        if format:
            plt.savefig("tenall_plot."+format)
        else:
            plt.show()
    
    def get_data(self):
        
        colors = {100: "tomato", 500: "dodgerblue", 1000: "springgreen"}
        markers = {100: "s", 500: "o", 1000: "D"}

        for popsize in colors:
            x = []
            y = []
            for data in self.stats_process.find({"parameters.popsize":popsize}):
                x.append(get_reference_array()[10]/data['ten'])
                y.append(get_reference_array()[len(GOOD_CANDS_LS)]/data['all'])
            plt.scatter(x, y, c=[colors[popsize]]*len(x), marker=markers[popsize], s=30, label="popsize="+str(popsize))
        
        plt.legend(loc=4)
        yxis = plt.ylim()
        plt.ylim(0.75, yxis[1])
        xxis = plt.xlim()
        plt.xlim(0.75, xxis[1])
        plt.xticks(np.arange(1,plt.xlim()[1], 1))
        plt.yticks(np.arange(1,plt.ylim()[1], 1))
        
        plt.annotate("best GA", xy= (9.5, 4.2), xytext=(9, 4.3), arrowprops={"arrowstyle":"->", "connectionstyle":"arc"})

class LSOSPlot():

    def __init__(self, format=None):
        plt.figure(1, figsize=(8,6))
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        
        self.db_os = Stats_Database(extension="_OS")
        self.stats_process_os = self.db_os._stats_process
        
        print self.stats_process_os.count()
        
        self.fontname = "Trebuchet MS"
        self.fontsize = 14
        
        self.get_data()
        
        plt.xlabel("Efficiency (one-photon light splitters)", fontname=self.fontname, fontsize=self.fontsize)
        plt.ylabel("Efficiency (photoanode oxide shields)", fontname=self.fontname, fontsize=self.fontsize)
        
        if format:
            plt.savefig("lsos_plot."+format)
        else:
            plt.show()
    
    def get_data(self):
        x = []
        y = []
        for data in self.stats_process.find({}):
            unique_key = data['unique_key']
            unique_key = unique_key.replace("eval_fitness_complex", "eval_fitness_complex_oxide_shield")
            unique_key = unique_key.replace("eval_fitness_simple", "eval_fitness_simple_oxide_shield")
            data2 = self.stats_process_os.find_one({"unique_key":unique_key})
            
            x.append(get_reference_array()[len(GOOD_CANDS_LS)]/data['all'])
            y.append(get_reference_array_OS()[len(GOOD_CANDS_OS)]/data2['all'])
        plt.scatter(x, y, c="dodgerblue", s=30)
        yxis = plt.ylim()
        plt.ylim(0.75, yxis[1])
        xxis = plt.xlim()
        plt.xlim(0.75, xxis[1])
        plt.xticks(np.arange(1,plt.xlim()[1], 1))
        plt.yticks(np.arange(1,plt.ylim()[1], 1))


        
class DataTable():
    
    def __init__(self):
        self.db = Stats_Database()
        self.stats_process = self.db._stats_process
        
        for it in self.stats_process.find():
            p = it['parameters']
            
            print ("{}\t{}\t{}\t{}\t{}\t{}\t{}").format(p['popsize'], get_pretty_name(p['selection_fnc']), get_pretty_name(p['fitness_fnc']),get_pretty_name(p['crossover_fnc']), p['elitism_num'], it['ten'], it['all'])
    
if __name__ == "__main__":
    format = None

    #print get_reference_array()[10]/9.5
    #print get_reference_array()[20]/3671
    
    if format:
        mpl.rcParams['savefig.dpi'] = 160
        #LSOSPlot(format=format)
        HeatMapPlot(format=format)
        #BreakoutPlot(format=format)
        #PerformancePlot(format=format)
        #ComparisonPlot(format=format)
        #ParametersPlot(format=format)
        #PerformancePlotExclusion(format=format)
        #TenAllPlot(format=format)
    
    else:
        #LSOSPlot()
        #TenAllPlot()
        #PerformancePlot()
        #PerformancePlotExclusion()
        #ComparisonPlot()
        #ParametersPlot()
        #BreakoutPlot()
        HeatMapPlot()
    
    # DataTable()
