'''
8/12/2015
this program just plots the flat and normalised flat alond the line
'''

#from __future__ import division, print_function

__author__ = "<isimion@ast.cam.ac.uk>"

# Standard library
import os, sys
import logging
import glob
import pyfits,numpy,pylab
from scipy import interpolate

# Third-party
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import leastsq

plot_path = "/Users/iuliasimion/Observing2014/hac_spectra/allyson/data_reduction_example/plots/"
def Openfile(i):
   return pyfits.open(i)[0].data.copy()

arc = Openfile('/Users/iuliasimion/Observing2014/hac_spectra/allyson/data_reduction_example/data/n2.0021.fit')

#take the median of columns 
arc = numpy.median(arc[:,100:200],1)
x = numpy.arange(arc.size)

fig = plt.figure(figsize=(18.5, 2.3))
plt.subplots_adjust(left=0.17, bottom=0.3, right=0.90, top=0.90, wspace=0.2)
ax = fig.add_subplot(111)
plt.plot(x, arc, lw = 2)
plt.xlim(0, arc.shape[0])
#plt.title('CCD row', fontsize = 15)
plt.xlabel('X pixel (dispersion direction)', fontsize = 20)
plt.ylabel('Intensity', fontsize = 20)
plt.show()
fig.savefig(plot_path+'arc.eps')

