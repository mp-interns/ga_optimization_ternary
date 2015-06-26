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

import pylab

chromnames = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
values = pylab.zeros((24, 24))

values[2][6] = 3
values[6][2] = 3
values[23][1] = 13
values[1][23] = 13
values[5][2] = 5
values[2][5] = 5
values[12][14] = 11
values[14][12] = 11
values[19][9] = 6
values[9][19] = 6


f = pylab.Figure()

#pylab.pcolor(values)
#imagename = "pcolor.png"

# pylab.imshow(values)
# imagename = "imshow.png"

pylab.matshow(values)
imagename = "matshow.png"


pylab.xticks(range(len(chromnames)), chromnames, rotation=90)
pylab.yticks(range(len(chromnames)), chromnames)

pylab.savefig(imagename)
