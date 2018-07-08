from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt

import os
import shutil
import numpy as np
import glob
import uuid
from astropy.table import Table
from astropy.io import fits
import pyfits, pylab
import matplotlib.pyplot as plt
import collections


'''
run in /Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/n2/reduced
plots the sky lines OI at 5577 and 6300
'''
homedir = os.path.expanduser("~")
night = 'n2'
data_path = os.path.expanduser("/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/"+night+'/')
data_path_reduced = os.path.expanduser("/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/"+night+'/'+'reduced')

#read rv
#rvfile = np.genfromtxt(os.path.join(data_path_reduced+'/vhelio', 'HAC177_rvcorrect.stdout')) 

def Openfile(i):
   return pyfits.open(i)[0].data.copy()

def check_files():
    #check what you have in the folder
    files = glob.glob(os.path.join(data_path, night+".????.fit"))# glob.glob('/data/n2.????.fit')
    files.sort()
    for i, f in enumerate(files):
       h = pyfits.open(f)[0].header
       print i, h['OBJECT'],h['EXPTIME']

flat_frames = []
bias_frames = []
object_files = []
arc_frames = []
targets = []
for filename in glob.glob(os.path.join(data_path, night+"*.fit")):
    hdr = fits.getheader(filename, 0)    
    obje = hdr['OBJECT']
    try:
        if obje == "FLAT":
            flat_frames.append(filename)
        elif obje == "ZERO":
            bias_frames.append(filename)
        elif obje == "COMP":
            arc_frames.append(filename)
        else:
            object_files.append(night+filename.lstrip(data_path))
            targets.append(obje)
    except:
        print "finished reading"


def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

def plot_spec(targets, all = 'no', continuum = 'yes'):
	target_names = f7(targets)

	if all == 'no':
		for i in range(0,len(target_names)):
			f = target_names[i] +'_'+str(targets.count(target_names[i]))+'.cal.fits'
			#read data
			spec = pyfits.open(f)[0].data.copy()
			hdulist=fits.open(f)
			h = pyfits.open(f)[0].header
			dw=h['CRVAL1']+np.arange(0, spec.shape[2],1)*h['CD1_1']
			df=spec[0,0,:]

			# Plot template and data
			fig = plt.figure(figsize=(10.5, 7))
			plt.title(target_names[i]  +' '+str(targets.count(target_names[i])) + 'exp.')
			plt.plot(dw, df, 'r.-')
			plt.show()
			fig.savefig(target_names[i]+'.eps')

	if all == 'yes':
		fig = plt.figure(figsize=(10.5, 7))
		for i in range(0,len(target_names)):
			f = target_names[i] +'_'+str(targets.count(target_names[i]))+'.cal.fits'
			#read data
			spec = pyfits.open(f)[0].data.copy()
			hdulist=fits.open(f)
			h = pyfits.open(f)[0].header
			dw=h['CRVAL1']+np.arange(0, spec.shape[2],1)*h['CD1_1']
			df=spec[0,0,:]
			print spec.shape


			#data
			plt.plot(dw, df, '-', label = target_names[i])
		plt.legend()
		plt.show()
		#fig.savefig('night_sky.eps')


	if all == 'yes':
		fig = plt.figure(figsize=(5.5, 3))
		ax = fig.add_subplot(1,1,1)
		plt.subplots_adjust(left=0.17, bottom=0.2, right=0.90, top=0.90, wspace=0.2)
		for i in range(1,len(target_names)-10):
			f = 'vh_'+target_names[i] +'_'+str(targets.count(target_names[i]))+'.cal.fits'
			#read data
			spec = pyfits.open(f)[0].data.copy()
			hdulist=fits.open(f)
			h = pyfits.open(f)[0].header
			dw=h['CRVAL1']+np.arange(0, spec.shape[2],1)*h['CD1_1']
			df=spec[2,0,:]


			plt.plot(dw, df/max(df), 'o-', label = target_names[i])
			if (i == 25) | (i ==26):
				plt.plot(dw, df/max(df), 'o-', label = target_names[i], lw = 10)
		#data
		plt.plot([5889.951, 5889.951], [0,0.8], '--', color = 'black')
		plt.plot([5895.924, 5895.924], [0,0.8], '--', color = 'black')
		plt.plot([6300.304, 6300.304], [0,0.8], '--', color = 'black')
		plt.plot([5577.338, 5577.338], [0,0.8], '--', color = 'black')
		plt.plot([4861.33, 4861.33], [0,0.8], '--', color = 'black')
		plt.plot([4861.33-4, 4861.33-4], [0,0.8], '--', color = 'red')
		plt.plot([4861.33+4, 4861.33+4], [0,0.8], '--', color = 'red')#error 250km/s

		plt.plot([6562.80, 6562.80], [0,0.8], '--', color = 'black')
		plt.plot([6562.80-5.5, 6562.80-5.5], [0,0.8], '--', color = 'red')
		plt.plot([6562.80+5.5, 6562.80+5.5], [0,0.8], '--', color = 'red')

		plt.plot([4340.47, 4340.47], [0,0.8], '--', color = 'black')
		ax.set_xlim(5570, 5585)
		ax.set_xlabel('$\lambda$ (A)')
		ax.set_ylabel('Relative intensity')

		ax.set_title('[OI] - 5577 A')
		ax.set_ylim(0, 1.05)

		plt.legend(fontsize = 10)
		plt.show()
		fig.savefig('OI5577.eps')

		fig = plt.figure(figsize=(5.5, 3))
		ax = fig.add_subplot(1,1,1)
		plt.subplots_adjust(left=0.17, bottom=0.2, right=0.90, top=0.90, wspace=0.2)
		for i in range(1,len(target_names)-10):
			f = 'vh_'+target_names[i] +'_'+str(targets.count(target_names[i]))+'.cal.fits'
			#read data
			spec = pyfits.open(f)[0].data.copy()
			hdulist=fits.open(f)
			h = pyfits.open(f)[0].header
			dw=h['CRVAL1']+np.arange(0, spec.shape[2],1)*h['CD1_1']
			df=spec[2,0,:]


			plt.plot(dw, df/max(df), 'o-', label = target_names[i])
			if (i == 25) | (i ==26):
				plt.plot(dw, df/max(df), 'o-', label = target_names[i], lw = 10)
		#data
		plt.plot([5889.951, 5889.951], [0,0.8], '--', color = 'black')
		plt.plot([5895.924, 5895.924], [0,0.8], '--', color = 'black')
		plt.plot([6300.304, 6300.304], [0,0.8], '--', color = 'black')
		plt.plot([5577.338, 5577.338], [0,0.8], '--', color = 'black')
		plt.plot([4861.33, 4861.33], [0,0.8], '--', color = 'black')
		plt.plot([4861.33-4, 4861.33-4], [0,0.8], '--', color = 'red')
		plt.plot([4861.33+4, 4861.33+4], [0,0.8], '--', color = 'red')#error 250km/s

		plt.plot([6562.80, 6562.80], [0,0.8], '--', color = 'black')
		plt.plot([6562.80-5.5, 6562.80-5.5], [0,0.8], '--', color = 'red')
		plt.plot([6562.80+5.5, 6562.80+5.5], [0,0.8], '--', color = 'red')

		plt.plot([4340.47, 4340.47], [0,0.5], '--', color = 'black')
		ax.set_xlim(6293, 6307)
		ax.set_ylim(0, 0.6)

		ax.set_xlabel('$\lambda$ (A)')
		ax.set_title('[OI] - 6300 A')
		ax.set_ylabel('Relative intensity')

		plt.legend(fontsize = 10)
		plt.show()
		fig.savefig('OI6300.eps')


	if continuum == 'yes':
		fig = plt.figure(figsize=(10.5, 7))
		ax = fig.add_subplot(1,1,1)
		for i in range(0,len(target_names)-1):
			f = 'c_'+'vh_'+target_names[i] +'_'+str(targets.count(target_names[i]))+'.cal.fits'
			#read data
			spec = pyfits.open(f)[0].data.copy()
			hdulist=fits.open(f)
			h = pyfits.open(f)[0].header
			dw=h['CRVAL1']+np.arange(0, spec.shape[2],1)*h['CD1_1']
			df=spec[0,0,:]


			plt.plot(dw, df, '-', label = target_names[i])
			if (i == 25) | (i ==26):
				plt.plot(dw, df, '-', label = target_names[i], lw = 10)
		#data
		plt.plot([5889.951, 5889.951], [0,0.8], '--', color = 'black', label='NaI')
		plt.plot([5895.924, 5895.924], [0,0.8], '--', color = 'black', label='NaI')
		plt.plot([6300.304, 6300.304], [0,0.8], '--', color = 'black', label='OI')
		plt.plot([5577.338, 5577.338], [0,0.8], '--', color = 'black', label='OI')
		plt.plot([4861.33, 4861.33], [0,0.8], '--', color = 'black', label='Hbeta')
		plt.plot([4861.33-4, 4861.33-4], [0,0.8], '--', color = 'red')
		plt.plot([4861.33+4, 4861.33+4], [0,0.8], '--', color = 'red')#error 250km/s

		plt.plot([6562.80, 6562.80], [0,0.8], '--', color = 'black', label='Halpha')
		plt.plot([6562.80-5.5, 6562.80-5.5], [0,0.8], '--', color = 'red', label='Halpha')
		plt.plot([6562.80+5.5, 6562.80+5.5], [0,0.8], '--', color = 'red', label='Halpha')

		plt.plot([4340.47, 4340.47], [0,0.8], '--', color = 'black', label='Hgamma')
		ax.set_xlim(4200, 6700)
		ax.set_ylim(0.2, 1.2)

		plt.legend()
		plt.show()
		#fig.savefig('night_sky_scaled.eps')
	return None

plot_spec(targets, all = 'yes', continuum = 'no')