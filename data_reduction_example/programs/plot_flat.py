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

flat = Openfile('/Users/iuliasimion/Observing2014/hac_spectra/allyson/data_reduction_example/data/Flat.fits')
nFlat = Openfile('/Users/iuliasimion/Observing2014/hac_spectra/allyson/data_reduction_example/data/nFlat.fits')

#take the median of columns 
flatnorm = numpy.median(flat,1)
x = numpy.arange(flatnorm.size)
model = interpolate.splrep(x,flatnorm,t=x[25:-25:50])
omodel = interpolate.splev(x,model)

flatnorm = (flat.T/omodel).T

#plot flatnorm along a line after you take the median over all the columns
flatnorm = numpy.median(flatnorm,1)
x = numpy.arange(flatnorm.size)

#normalised flat created with pyraf
flatnorm_pyraf = numpy.median(nFlat,1)
x_pyraf = numpy.arange(flatnorm_pyraf.size)

fig = plt.figure(figsize=(5.5, 3.5))
plt.subplots_adjust(left=0.12, bottom=0.2, right=0.90, top=0.90, wspace=0.2)
ax = fig.add_subplot(111)
plt.plot(x, flatnorm, lw = 2)
#plt.plot(x_pyraf, flatnorm_pyraf, lw = 2, color = 'red') --looks slightly different 
#plt.legend()
plt.xlim(0, flat.shape[0])
plt.title('normalised flat', fontsize = 15)
plt.xlabel('X', fontsize = 15)
plt.ylim([0.94,1.06])
plt.show()
fig.savefig(plot_path+'flatnorm.eps')


fig = plt.figure(figsize=(5.5, 3.5))
plt.subplots_adjust(left=0.12, bottom=0.2, right=0.90, top=0.90, wspace=0.2)
ax = fig.add_subplot(111)
flatnorm = numpy.median(flat,1)
plt.plot(x, flatnorm, lw = 2)
plt.plot(x, omodel, lw = 1)
#plt.legend()
plt.xlim(0, flat.shape[0])
plt.title('flat', fontsize = 15)
plt.xlabel('X', fontsize = 15)
plt.show()
fig.savefig(plot_path+'flat.eps')
