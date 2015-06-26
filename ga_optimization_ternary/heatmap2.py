#!/usr/bin/env python

'''
Created on Oct 31, 2012
'''

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Oct 31, 2012"

import matplotlib.pyplot as plt
import numpy as np
import random

R = xrange(10000)
data = [random.choice(R) for i in range(25)]
data = np.array(data)
data.shape = (5,5)

fig = plt.figure()
x = ['a','b','c','d','e']
plt.xticks([0.5, 1.5, 2.5, 3.5, 4.5], x, ha='center')
plt.yticks([0.5, 1.5, 2.5, 3.5, 4.5], x, ha='center')
plt.hot()
plt.pcolormesh(data)
plt.colorbar() 
plt.savefig('example.png')