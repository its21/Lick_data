# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 12:21:58 2015

@author: isimion

Collection of pyraf scripts primarily for reducing MDM spectra:
This is an example on night 2 of observations
for the other nights, the reduction steps were the same

Load pyraf (I installed UREKA - http://ssb.stsci.edu/ureka/ -, so I open an UREKA terminal and type $pyraf )
$ pyraf
$ cd /Users/iuliasimion/Observing2014/hac_spectra/allyson/data_reduction_example/programs
$ pyexecute data_reduction.py

If you modify it:
$ import data_reduction
$ reload(data_reduction)

#ATTENTION: for each new task you need to create a .par file (e.g. '/Users/iuliasimion/iraf/pars/normalise_flats.par')
which contains the line:
mode,s,h,'al'


###run these functions in order, one could automatise this
pyraf$: bias_correction()
pyraf$: combine_flats()
pyraf$: normalise_flats()
pyraf$: flat_correction()
pyraf$: cr_removal() #(this takes 30 mins and sigclip value is not really perfect) try to avoid running this by mistake
pyraf$: image_shift()
pyraf$: image_combine
#####(in the data_path_reduced):
pyraf$: run_apall() or run_apall_1exp() - I prefer 1exp so I use that in the template fitting
pyraf$: run_identify() or run_identify_1exp()
pyraf$: run_rvcorrect() or run_rvcorrect_1exp()
"""

from __future__ import division

import os
from pyraf import iraf
import shutil
import numpy as np
import glob
import uuid
from astropy.table import Table
from astropy.io import fits
#import pyfits, pylab
import matplotlib.pyplot as plt
import collections
from scipy import optimize

homedir = os.path.expanduser("~")
night = 'data-test/'
data_path = os.path.expanduser("/Users/iuliasimion/work/2018/observing/Lick/Lick_data/")
data_path_reduced = os.path.expanduser("/Users/iuliasimion/work/2018/observing/Lick/Lick_data/")
   
def Openfile(i):
   return fits.open(i)[0].data.copy()

files = glob.glob(os.path.join(data_path, night+"????.fits"))# glob.glob('/data/n2.????.fit')
files.sort()
for i, f in enumerate(files):
   h = fits.open(f)[0].header
   print(i, h['OBJECT'],h['EXPTIME'])

#160, 16     1, 162
obj = Openfile(files[160])# + Openfile(files[161]) + Openfile(files[162])
offslit = obj[2140:2170,220]
onslit = obj[2140:2170,206]
sn = np.mean(onslit)/np.std(offslit)

plt.imshow(obj, vmin = 990, vmax =1100)
plt.show()

sum_down = np.sum(obj, axis=0) #sum along columns  = len 295
idx = int(np.argmax(sum_down)) #central column of the spectrum 135

spectrum = np.sum(obj[:, 203:215], axis=1)

sumc = sum_down[203:215]
idx = int(np.argmax(sumc)) #central column of the spectrum 135


yv = sumc[idx-5:idx+5]
xv = np.arange(idx-5, idx+5,1)
#xn = np.sum(xv*yv)/np.sum(yv)
#width = np.sqrt(abs(np.sum((xv-xn)**2*yv)/np.sum(yv)))
