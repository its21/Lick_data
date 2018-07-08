import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
import glob
import os
from PyAstronomy import pyasl
from astropy.io import fits
import pyfits
import scipy as sc

from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
from scipy.interpolate import interp1d


#choose the night
night = 'n2'
data_path_reduced = os.path.expanduser("/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/"+night+'/'+'reduced')
data_path = os.path.expanduser("/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/"+night+'/')

c =  299792.458

def Openfile(i):
   return pyfits.open(i)[0].data.copy()

def check_files():
    #check what you have in the folder
    files = glob.glob(os.path.join(data_path, night+".????.fit"))# glob.glob('/data/n2.????.fit')
    files.sort()
    for i, f in enumerate(files):
       h = pyfits.open(f)[0].header
       #print i, h['OBJECT'],h['EXPTIME']
    return None

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
            if obje == "FLAT":
                flat_frames.append(filename)
            elif obje == "ZERO":
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

def read_spec(spec_file, var='no', sky='no'):
	spec = pyfits.open(spec_file)[0].data.copy()
	hdulist=fits.open(spec_file)	
	h = pyfits.open(spec_file)[0].header
	tw=h['CRVAL1']+np.arange(0, spec.shape[2],1)*h['CD1_1']
	if var =='yes':
		tf=spec[3,0,:]
	if sky =='yes':
		tf=spec[2,0,:]
	else:
		tf=spec[0,0,:]
	return tw, tf

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]


def model(p, x):
    A = p['A'].value
    mu = p['mu'].value
    sigma = p['sigma'].value
    #b = p['b'].value
    return A*np.exp(-0.5 * (x - mu)**2 / (2 * sigma**2)) 

   
def func(p, x, y):
    return y - model(p, x)
 

def calc_cal_shift():
    '''
    fit [OI] sky lines with Gaussians and find the centre and sigma of each line
    '''

    #OI lines used for finding the calibration errors
    lines = [5577, 6300]

    lines_mu1, lines_mu2 = [], []
    lines_sigma = []
    velocity_err = []

    flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()
    target_names = f7(targets)
    letters = set('C')
    targs = []
    for word in target_names:
        if letters & set(word):
            targs.append(word)
    for line in lines:
        for name in targs:
            #choose your target
            target_to_cc =os.path.join(data_path_reduced,'vh_'+name +'_'+str(targets.count(name))+'_med.cal.fits')
            w, f = read_spec(target_to_cc)
            w, var = read_spec(target_to_cc, var = 'yes') #variance
            w, sky = read_spec(target_to_cc, sky = 'yes') #variance

            #isolate the wavelentgh of the sky line
            w_skyline = w[(w > line-7) & (w < line+7)].copy()
            idx = int(np.argmax(sky[(w > (line-7)) & (w < (line+7))])) 
            w_max = w_skyline[idx]         #find the peak and the width
            sky_line = sky[(w > (line-7)) & (w < (line+7))]
            wn = np.sum(w_skyline*sky_line)/np.sum(sky_line)
            width = np.sqrt(abs(np.sum((w_skyline-wn)**2*sky_line)/np.sum(sky_line)))

            #fit to find parameters of the line
            params = Parameters()
            params.add('A', value = 0.5)
            params.add('mu', value = np.log10(line))
            params.add('sigma', value = 0.0002)
            #params.add('b', value= 50)

            result = minimize(func, params, args=(np.log10(w_skyline),sky_line), maxfev = 8000)
            fit = func(params, np.log10(w_skyline), sky_line)

            #reduced chi-square
            chi = result.redchi
            fwhm = 2.355*params['sigma'].value
            velocity = (10**fwhm - 1)*c #width of the line in km/s
            sigma = params['sigma'].value 

            if line  == lines[0]:
                lines_mu1.append(params['mu'].value)
            elif line  == lines[1]:
                lines_mu2.append(params['mu'].value)

            lines_sigma.append(params['sigma'].value)
            velocity_err.append(velocity)

    return lines_mu1, lines_mu2, lines_sigma 

lines_mu1, lines_mu2, lines_sigma  = calc_cal_shift()

def calibration_error():
    #find average velocity error introduced by the calibration and error on the average
    lines = [5577, 6300]
    lines_mu1, lines_mu2, lines_sigma = calc_cal_shift()

    lines_mu1 = np.array(lines_mu1)
    lines_mu2 = np.array(lines_mu2)
    lines_sigma = np.array(lines_sigma) #width of the line in log10 space
    lines_fwhm = 2.355*abs(lines_sigma)

    deltalogw =  list(lines_mu1-np.log10(lines[0])) + list((lines_mu2) - np.log10(lines[1]))
    deltalogw = np.array(deltalogw)
    vel = (10**deltalogw - 1)*c #same as 'DeltaLambda/Lambda *c '
    print 'delta v = ', vel

    ###calculate the calibration error
    nobins = 16
    h, edge = np.histogram(abs(np.array(vel)), nobins, range = [0,60])
    x = edge[:-1]+(edge[1]-edge[0])/2.
    params = Parameters()
    params.add('A', value= 100)
    params.add('mu', value= 20)
    params.add('sigma', value= 10)
    params.add('b', value= 0)

    result = minimize(func, params, args=(x,h), maxfev = 8000)
    #calibration error:
    sigma_vel_cal = params['sigma'].value #for n2 this is more correct (because one target, HAC210 that I eliminate 
        #in the results section because of phase is bad)
    #sigma_vel_cal = np.std(vel) #for other nights I also use this
    #velocity shift:
    mu_vel_cal = params['mu'].value
    print 'systematic shift = ', mu_vel_cal


    ###find the sigma of the LSF used to degrade the Munari+2005 templates
    #lines_sigma is the width of each sky line
    #taking the mean of the sigmas, will give you the average sigma you need to degrade the templates with a Gaussian convolution
    h, edge = np.histogram(np.array(lines_sigma), nobins)
    x = edge[:-1]+(edge[1]-edge[0])/2.
    params = Parameters()
    params.add('A', value= 100)
    params.add('mu', value= 0.0004)
    params.add('sigma', value= 0.001)
    params.add('b', value= 0)

    result = minimize(func, params, args=(x,h), maxfev = 8000)
    sigma_sigma = params['sigma'].value
    mean_sigma = params['mu'].value
    #mean_sigma = np.mean(lines_sigma)


    return sigma_vel_cal, mean_sigma


sigma_vel_cal, mean_sigma = calibration_error()