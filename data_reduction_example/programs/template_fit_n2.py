import numpy as np
import glob
import os
import scipy.signal as signal
from astropy.io import fits
import pyfits, pylab
import matplotlib.pyplot as plt
import collections
from scipy import optimize
import numdifftools as nd
import scipy
import numdifftools
import cPickle
import scipy as sc

'''
- calculates the velocity and the errors
'''

#choose the path for your main folder:
folder_path = "/Users/iuliasimion/Observing2014/hac_spectra/allyson/data_reduction_example/"

data_path = os.path.expanduser(folder_path+"data/")
data_path_reduced = os.path.expanduser(folder_path+"data_reduced/")
results_path = os.path.expanduser(folder_path+"results/") #for Allyson
#templates_path = os.path.expanduser(folder_path+"T3500-7250/") #for Allyson -- you need to download the library
templates_path = os.path.expanduser("/Users/iuliasimion/Observing2014/T3500-7250/")
plot_path = os.path.expanduser(folder_path+"plots/")
#I save here: results_path = os.path.expanduser("/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/results/") 


def files():
    flat_frames = []
    bias_frames = []
    object_files = []
    arc_frames = []
    targets = []
    arc_targets = []
    for filename in glob.glob(os.path.join(data_path, night+"*.fit")):
        hdr = fits.getheader(filename, 0)    
        obje = hdr['OBJECT']
        try:
            if (obje == "FLAT")  | (obje == "flat"):
                flat_frames.append(filename)
            elif (obje == "ZERO"):
                bias_frames.append(filename)
            elif obje == "COMP":
                arc_frames.append(filename)
                arc_targets.append(obje)
            else:
                object_files.append(night+filename.lstrip(data_path))
                targets.append(obje)
        except:
            print "finished reading"

    return flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets

def read_spec(spec_file, var='no'):
	spec = pyfits.open(spec_file)[0].data.copy()
	hdulist=fits.open(spec_file)	
	h = pyfits.open(spec_file)[0].header
	tw=h['CRVAL1']+np.arange(0, spec.shape[2],1)*h['CD1_1']
	if var =='yes':
		tf=spec[3,0,:]
	else:
		tf=spec[0,0,:]
	return tw, tf

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

def fit_template(w, f, dv, target, N = 15):
	'''
	w : wavelengths data, 1d array
	N : degree of the polynomial
	f : spectrum
	dv: sigma spectrum
	target : target name
	'''

	#feature scaling (-0.5<tw<0.5) : important because w ~ 1e3 and you have exponentials in the poly
	ww = w[((w<(np.log10(6300-15))) | (w>np.log10((6300+15)))) & ((w<np.log10((5577-9))) | (w>np.log10((5577+9)))) & ((w < np.log10((6868-9))) | (w > np.log10((6868+12)))) &  ((w < np.log10((6363-5))) | (w > np.log10((6363+5)))) ]
	ff = f[((w<(np.log10(6300-15))) | (w>np.log10((6300+15)))) & ((w<np.log10((5577-9))) | (w>np.log10((5577+9)))) & ((w < np.log10((6868-9))) | (w > np.log10((6868+12)))) &  ((w < np.log10((6363-5))) | (w > np.log10((6363+5))))  ]
	dvv = dv[((w<(np.log10(6300-15))) | (w>np.log10((6300+15)))) & ((w<np.log10((5577-9))) | (w>np.log10((5577+9)))) & ((w < np.log10((6868-9))) | (w > np.log10((6868+12)))) &  ((w < np.log10((6363-5))) | (w > np.log10((6363+5)))) ]
	dw_s = (ww-ww[len(ww)/2])/[ww[-1]-ww[0]]

	#matrix of powers for a poly of degree N
	c0 = np.arange(0,N+1,1)
	Powers = np.repeat(c0, len(ww), axis=0)
	Powers = Powers.reshape(N+1, len(ww))

	#Features matrix: n = N+1, m=len(w) (columns=m, row=n) -needs to be transposed
	A = np.power(dw_s, Powers)

	#velocity ranges
	rvmin, rvmax = -600, 400
	drv = 1. 
	drvs = np.arange(rvmin, rvmax, drv)

	#choose your templates grid
	temperature = np.arange(6000, 8250, 250)
	logg = np.arange(20,40,5)
	met = np.arange(5,25,5)
	best_vel = []

	#clean previous saved results
	#*****_templatefit_log_correct.txt
	savefile = target+'_templatefit_log.txt'
   	if os.path.exists( os.path.join(data_path, savefile) ):
		os.remove( os.path.join( data_path, savefile) )
		print 'Removing file ' + os.path.join( data_path, savefile)

	#Read the template
	counter = 0
	for T in temperature:
		for G in logg:
			for M in met:
				if M < 10:
					filename = os.path.join(templates_path+'T_0'+str(T)+"/", "T0"+str(T)+'G'+str(G)+"M0"+str(M)+"V000K2ANWNVD01F.ASC")
				else:
					filename = os.path.join(templates_path+'T_0'+str(T)+"/", "T0"+str(T)+'G'+str(G)+"M"+str(M)+"V000K2ANWNVD01F.ASC")
				#print filename
				templ_file = np.loadtxt(filename)
				fixpix = 0.00011
				logw_regular = np.arange(min(np.log10(tw)), max(np.log10(tw)), fixpix)
				#get your template on a regular grid
				templ = sc.interpolate.interp1d(np.log10(tw), templ_file)
				win = signal.gaussian(len(logw_regular), std = mean_sigma/fixpix)
				tf = signal.convolve(templ(logw_regular), win, mode='same') / sum(win)

				redchi = np.zeros(len(drvs))
				chi = np.zeros(len(drvs))

				restrict = ((ww>(np.log10(6562-50))) & (ww<np.log10((6562+50)))) | ((ww>np.log10((4861-50))) & (ww<np.log10((4861+50)))) 
				for i, rv in enumerate(drvs):
					#fi = sc.interpolate.interp1d(np.log10(tw)+np.log10(1.0 + rv/c), tf)
					fi = sc.interpolate.interp1d(logw_regular+np.log10(1.0 + rv/c), tf)
					AT = np.transpose(A*fi(ww)/dvv)
					coeff_poly, residuals = scipy.linalg.lstsq(AT, ff/dvv)[0] , scipy.linalg.lstsq(AT, ff/dvv)[1]

					rank  = scipy.linalg.lstsq(AT, ff/dvv)[2]
					singular  = scipy.linalg.lstsq(AT, ff/dvv)[3]
					redchi[i] = residuals/(len(ww[restrict])-N)
					chi[i] = residuals

				velindx = np.argmin(chi)
				vel = drvs[velindx]
				#templ = sc.interpolate.interp1d(logw_regular+np.log10(1.0 + vel/c), tf)
				#save values to a file
				myfile = open(os.path.join(results_path, savefile), 'a') 
				myfile.write(str(T) + '\t' +str(G)+'\t'+str(M)+'\t'+str(vel)  +'\t'+str(min(chi))+'\t'+str(min(redchi)) + '\n') #add also the error
				myfile.close()

				counter = counter + 1
				print 'fitting template', counter,'/'+str(len(temperature)*len(logg)*len(met))+' for object', target
	return None

#read wavelengths template (features)
tw = np.loadtxt(os.path.join(templates_path, "LAMBDA_D1A.dat"))

def fit_templates():
	'''
	fit the templates grid for each target
	'''
	
	flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()
	target_names = f7(targets)
	print 'targets:', target_names
	print 'object_files:', object_files

	for j in range(1,len(target_names)):
		fd_new = []
		dvd_new = []
		for filename in object_files[1:-10]: #[:-3] for n4
			filen = 'vh_s_cr_f_'+filename.rstrip('.fit')+'.1d.cal.fits'   
			target_to_cc =os.path.join(data_path_reduced,filen)     
			hdr = fits.getheader(target_to_cc, 0) 
			obje = hdr['OBJECT']
			if target_names[j] == obje:
				print target_names[j], filen
				wd, fd = read_spec(target_to_cc)
				dw, dvd = read_spec(target_to_cc, var = 'yes') #variance
				
				#print 'appending', filen, target_names[j]
				fd_new.append(fd)
				dvd_new.append(dvd)
				wd = np.log10(wd)
				dw = np.log10(dw)
			else:
				continue #print target_names[j],'not found'
		dvd_new = np.array(dvd_new)
		fd_new = np.array(fd_new)
		fd_f, dvd_f = 0, 0
		for ss in range(len(fd_new)):
			fd_f += fd_new[ss]#+fd_new[1] +fd_new[2] 
			dvd_f += dvd_new[ss]**2.#+dvd_new[1]**2.+dvd_new[2]**2.
		restrict = ((wd<(np.log10(6562-50))) | (wd>np.log10((6562+50)))) & ((wd<np.log10((4861-50))) | (wd>np.log10((4861+50)))) 
		dvd_f[restrict] = 100*dvd_f[restrict]
		dvd_f[dvd_f==0]=1
		print 'fitting night', night, 'with sigma and mean', sigma_vel_cal, mean_sigma
		fit_template(wd, fd_f, np.sqrt(dvd_f), target_names[j], N = 15)
	return None

c =  299792.458

sigmas = [12.723465477714441, 8.3742481968970246 , 9.7718080194882226, 14.172317339249838, 9.0600211435116105, 3.8823774213063063]
means = [-9.4653852149442552e-05,  -9.2634790369026325e-05, -8.8542384544618228e-05, -3.6228700237005092e-06, -3.8056416407812281e-06, -9.914964525885501e-05]

#choose the night you want to reduce 'n1', 'n2', 'n3' , 'n4' , 'n5', 'n6'
n=2
night = 'n'+str(n)

sigma_vel_cal, mean_sigma = sigmas[n-1], means[n-1]
night = 'n2'

fit_templates()
