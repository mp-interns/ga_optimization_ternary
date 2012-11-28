#!/usr/bin/env python

'''
Created on Oct 28, 2012
'''
from __future__ import division

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Oct 28, 2012"

import matplotlib.pyplot as plt
import numpy as np
import math

def gaussian_pdf(x, mean=0):
    return (1/math.sqrt(2*math.pi))*math.exp(-0.5*(x-mean)*(x-mean))

class FitnessPlot():

    def __init__(self, format=None):
        plt.figure(1, figsize=(16, 8))
        
        self.fontname = "Trebuchet MS"
        self.lw = 4
        self.fontsize = 14
        
        params = [("Formation E (eV)", "Discontinuous"), ("Formation E (eV)", "Smooth"), ("Band Gap (eV)", "Discontinuous"), ("Band Gap (eV)", "Smooth"), ("Band edge distance (eV)", "Discontinuous"), ("Band edge distance (eV)", "Smooth")]
        plt.subplot(3, 2, 1)
        for idx, param in enumerate(params):
            plt.subplot(3, 2, idx+1)
            self.get_data(param)
        
        plt.subplots_adjust(left=0.125, bottom=None, right=None, top=None, wspace=0.25, hspace=0.5)
        
        if format:
            plt.savefig("fitness_plot."+format)
        else:
            plt.show()
    
    def get_data(self, parameter):
        
        x = [1, 2, 3, 4, 5]
        y = [1, 2, 3, 4, 5]
        
        if parameter[0] == "Formation E (eV)":
            x = get_interval(-1, 2, 0.01)
            if parameter[1] == "Smooth":
                y = [heat_of_formation_complex(val) for val in x]
            else:
                y = [heat_of_formation_simple(val) for val in x]
        
        if parameter[0] == "Band Gap (eV)":
            x = get_interval(0, 5, 0.01)
            if parameter[1] == "Smooth":
                y = [gap_complex(val) for val in x]
            else:
                y = [gap_simple(val) for val in x]
        
        if parameter[0] == "Band Gap (eV)":
            x = get_interval(0, 5, 0.01)
            if parameter[1] == "Smooth":
                y = [gap_complex(val) for val in x]
            else:
                y = [gap_simple(val) for val in x]

        if parameter[0] == "Band edge distance (eV)":
            x = get_interval(-1, 2, 0.01)
            if parameter[1] == "Smooth":
                y = [edge_complex(val) for val in x]
            else:
                y = [edge_simple(val) for val in x]
                                
        plt.xlabel(parameter[0], fontname=self.fontname, fontsize=self.fontsize)
        plt.ylabel("Fitness "+parameter[1], fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_xticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.setp(plt.gca().get_yticklabels(), fontname=self.fontname, fontsize=self.fontsize)
        plt.errorbar(x, y, lw=self.lw, color="blue")
        plt.ylim((0, max(y)+1))
        #plt.xlim((0, NUM_CANDS))
    

def heat_of_formation_complex (heat_of_formation):
    if heat_of_formation <= 0.2:
        return 10
    else:
        return 20 * (1-1/(1+math.exp(((-heat_of_formation) + 0.2) * 3.5)))


def heat_of_formation_simple(heat_of_formation):
    if heat_of_formation <= 0.2:
        return 10
    
    if heat_of_formation <= 0.5:
        return 5
    
    return 0


def gap_complex(gap_dir):
    if gap_dir == 0:
        return 0
    if (gap_dir >= 1.5 and gap_dir <= 3):
        return 10
    else:
        return 33 * gaussian_pdf(gap_dir, 2.25)


def gap_simple(gap_dir):
    if (gap_dir >= 1.5 and gap_dir <= 3):
        return 10
    return 0


def edge_complex(cb_dir):
    if cb_dir <= 0:
        return 5
    else:
        distance = cb_dir * 5
        return 10 * (1-1/(1+math.exp(-distance)))


def edge_simple(cb_dir):
    if cb_dir <= 0:
        return 5
    return 0

                          
def get_interval(min, max, interval):
    counter = min
    m_a = []
    while counter <= max:
        m_a.append(counter)
        counter += interval
    return m_a

if __name__ == "__main__":
    format = "png"
    
    if format:
        FitnessPlot(format=format)
    
    else:
        FitnessPlot()
    
    
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