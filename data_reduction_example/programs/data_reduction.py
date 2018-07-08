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
import pyfits, pylab
import matplotlib.pyplot as plt
import collections
from scipy import optimize

homedir = os.path.expanduser("~")
night = 'n2'
data_path = os.path.expanduser("/Users/iuliasimion/Observing2014/hac_spectra/allyson/data_reduction_example/data/")
data_path_reduced = os.path.expanduser("/Users/iuliasimion/Observing2014/hac_spectra/allyson/data_reduction_example/data_reduced")

   
def Openfile(i):
   return pyfits.open(i)[0].data.copy()

def check_files():
    #check what you have in the folder
    files = glob.glob(os.path.join(data_path, night+".????.fit"))# glob.glob('/data/n2.????.fit')
    files.sort()
    for i, f in enumerate(files):
       h = pyfits.open(f)[0].header
       print i, h['OBJECT'],h['EXPTIME']

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

parfile = iraf.osfn( os.path.join(homedir,'/Users/iuliasimion/iraf/pars/files.par') )
t = iraf.IrafTaskFactory(taskname="files", value=parfile, function=files)

# fit a Gaussian to summed 1D profile to find width
def model(p, x):
    A, sigma, mu, b = p
    return A*np.exp(-0.5 * (x - mu)**2 / (2 * sigma**2)) + b
    
def func(p, x, y):
    return y - model(p, x)

def shift_files():
    '''
    calculates the shift in trace in between exposures and plots the fit
    input: it needs the targets and object_files calculated above
    output:
     - input_shift.csv: the list of flat and bias corrected objects (all exposures)
     - output_shift.csv: the list of shifted objects
     - shift.csv: the shift for each object
    cr = 'no' : you shift the images with cosmic rays (no lacos applied)
    cr = 'yes': you shift the cr cleaned images
    '''
    flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()
    target_names = f7(targets)

    if os.path.exists(os.path.join(data_path, "input_shift.csv")):
        os.remove(os.path.join(data_path, "input_shift.csv"))
        print "removing input_shift.csv"
    if os.path.exists(os.path.join(data_path_reduced, "output_shift.csv")):
        os.remove(os.path.join(data_path_reduced, "output_shift.csv"))
        print "removing output_shift.csv"
    if os.path.exists(os.path.join(data_path, "shift.csv")):
        os.remove(os.path.join(data_path, "shift.csv"))
        print "removing shift.csv"


    for item, count in collections.Counter(targets).items():
        if count > 1:
            print (item,count)
            center = [] 
            count_exposure = 0

            plt.clf()

            for filename in glob.glob(os.path.join(data_path, night+"*.fit")):
                hdr = fits.getheader(filename, 0)    
                obje = hdr['OBJECT']
                
                if obje == item:
                    #    filen = 'cr_'+'f_'+night+filename.lstrip(data_path)+'s'
                    #if targets.count(item) > 2:
                    #    filen = 'f_'+night+filename.lstrip(data_path)+'s'
                    filen = 'cr_'+'f_'+night+filename.lstrip(data_path)+'s'
                    if os.path.exists(os.path.join(data_path, 's_'+filen)):
                        os.remove(os.path.join(data_path, 's_'+filen))
                        print "removing ",os.path.join(data_path,'s_'+filen)
                    count_exposure = count_exposure+1
                    print count_exposure
                    #if count_exposure > 1:
                    print 'writing exposure',count_exposure, filen
                    myfile = open(os.path.join(data_path, 'input_shift.csv'), 'a') 
                    myfile.write(os.path.join(data_path, filen)+ '\n') 
                    myfile.close()
                    myfile = open(os.path.join(data_path_reduced, 'output_shift.csv'),'a') 
                    myfile.write(os.path.join(data_path_reduced, 's_'+filen) + '\n') 
                    myfile.close()

                    obj = Openfile(os.path.join(data_path, filen))
                    sum_down = np.sum(obj, axis=0) #sum along columns  = len 295
                    sumc = sum_down[130:160]
                    idx = int(np.argmax(sumc)) #central column of the spectrum 135
                    
                    yv = sumc[idx-5:idx+5]
                    xv = np.arange(idx-5, idx+5,1)
                    print 'lengths', idx, len(sumc), len(yv), len(xv)
                    #xn = np.sum(xv*yv)/np.sum(yv)
                    #width = np.sqrt(abs(np.sum((xv-xn)**2*yv)/np.sum(yv)))

                    #keep this
                    p_opt, ier = optimize.leastsq(func,x0=[100000,2,idx,1000],args=(xv,yv))
                    A, sigma, mu, b = p_opt
                    center.append(mu)
                    #print p_opt

                    #####print and plot
                    #print center
                    #print item, idx, mu, item, filen
                    fit = lambda t : A*np.exp(-0.5*(t-mu)**2/(2*sigma**2))+b
                    
                    plt.plot(xv,yv, lw = 3, color = 'red', label = 'data')
                    plt.plot(np.arange(idx-7, idx+7,0.1),fit(np.arange(idx-7, idx+7,0.1)), ls = '-.',label = filen)
            center = np.array(center)
            deltay = center[1:]-center[0]
            #you don't want a shift between frames larger than 2pixels
            print 'READ. if this is < 2, it is okay! DELTAY=',deltay
            for i in range(-1,len(deltay)):
                if i == -1:
                    myfile = open(os.path.join(data_path, 'shift.csv'),'a') 
                    myfile.write(str(0) + ' 0.0' + '\n') #doesnt shift the first exposure
                    myfile.close()
                else:
                    myfile = open(os.path.join(data_path, 'shift.csv'),'a') 
                    myfile.write(str(-deltay[i]) + ' 0.0'+ '\n') 
                    myfile.close()

            #plt.legend()
            #plt.show()
    return None

parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/shift_files.par') )
t = iraf.IrafTaskFactory(taskname="shift_files", value=parfile, function=shift_files)


def all_files_list():
    #if you re running this it means you want to delete the previous files
    #run it in the main folder
    if os.path.exists('input.csv'):
        os.remove('input.csv')
    if os.path.exists('output.csv'):
        os.remove('output.csv')

    for filename in glob.glob(os.path.join(data_path, night+"*.fit")):
        #original files
        myfile = open('input.csv','a') #write it in the directory where u run it
        myfile.write(filename + '\n') 
        myfile.close()
        #bias corrected and trimmed files
        myfile = open('output.csv','a') 
        myfile.write(filename +'s'+ '\n') 
        myfile.close()
    return



def bias(bias_frames):  
    '''
    *improve: run imstat on the bias frames and discard the ones with st dev too big or max too big 
    $ imstat bias.fit 
    the std dev should be = Read noise/ gain / sqrt(no of bias frames)
    our value is bigger but it's okay because we do not clean the CR at this point and the effect is not big 
    we get:
    --> imstat bias.fit #(in the data folder)
    #               IMAGE      NPIX      MEAN    STDDEV       MIN       MAX
             bias.fit    604136     4.733     2.209    -3.293     12.38

    '''
    bias = []
    for i in bias_frames:
        bias.append(Openfile(i))
    bias = np.median(bias,0)
    pyfits.PrimaryHDU(bias).writeto(os.path.join(data_path, "bias.fit"),clobber=True)
    return

def bias_correction():
    '''
    trim images 
    correct all images for bias using the master bias and the overscan region in each image
    $ bias_correction() 
    $ epar ccdproc
    '''
    #iraf packages
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)

    #check if the master bias exists, and if it does, remove it
    if os.path.exists(os.path.join(data_path, "bias.fit")):
        os.remove(os.path.join(data_path, "bias.fit"))

    #check if the bias corrected files exist already, and if they do, remove them
    for filename in glob.glob(os.path.join(data_path, night+"*.fits")):
        if os.path.exists(filename):
            os.remove(filename)
            
    #generate lists with input and output file names
    #if (os.path.isfile('input.csv') == False):
    all_files_list()

    #create the master bias
    bias_frames = []
    for filename in glob.glob(os.path.join(data_path, "n*")):
        hdr = fits.getheader(filename, 0)    
        obje = hdr['OBJECT']
        if obje == "ZERO":
            bias_frames.append(filename)
    bias(bias_frames)

    #settings
    iraf.ccdproc.setParam('images','@input.csv')
    iraf.ccdproc.setParam('output','@output.csv')
    iraf.ccdproc.setParam('overscan','yes')
    iraf.ccdproc.setParam('trim','yes')
    iraf.ccdproc.setParam('zerocor','yes')
    iraf.ccdproc.setParam('flatcor','no')
    iraf.ccdproc.setParam('readaxis','line')
    iraf.ccdproc.setParam('biassec',"["+"305:360,5:2045"+"]")
    iraf.ccdproc.setParam('trimsec',"["+"5:300,5:2045"+"]")
    iraf.ccdproc.setParam('zero',os.path.join(data_path, "bias.fit"))
    iraf.ccdproc.setParam('flat','')
    iraf.ccdproc.setParam('interactive','no')
    
    iraf.ccdproc()
    return None

    

parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/bias_correction.par') )
t = iraf.IrafTaskFactory(taskname="bias_correction", value=parfile, function=bias_correction)


def combine_flats():

    """
    Co-add flat-field images. Check no variation night to night and then add all together.

    Once combined should use imarith to divide each by the average and check see only
    noise and no residuals

    First should lcpixmap flat fields

    """

    #iraf packages
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)

    #print 'In directory ' + flatdir
    #print 'Combining flats...'

    if os.path.exists('Flat.fits'):
        os.remove('Flat.fits')

    if os.path.exists('flats.csv'):
        os.remove('flats.csv')

    for filename in glob.glob(os.path.join(data_path, night+"*.fits")):
        hdr = fits.getheader(filename, 0)    
        obje = hdr['OBJECT']
        if obje == "FLAT":
            myfile = open('flats.csv','a') 
            myfile.write(filename + '\n') 
            myfile.close()
            
#    # Give write permission to files
#    for name in names:
#        os.chmod(name.replace('[1]',''),0644)

    iraf.flatcombine.setParam('input', '@flats.csv' )
    iraf.flatcombine.setParam('rdnoise',7.9)
    iraf.flatcombine.setParam('gain',2.7)
    iraf.flatcombine.setParam('combine','average')
    iraf.flatcombine.setParam('reject','crreject')
    iraf.flatcombine.setParam('ccdtype','')
    iraf.flatcombine.setParam('scale','mode')
    iraf.flatcombine.setParam('lsigma',3.0)
    iraf.flatcombine.setParam('hsigma',3.0)
    iraf.flatcombine.setParam('output',os.path.join(data_path, "Flat.fits"))

    iraf.flatcombine()

    return None

parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/combine_flats.par' ) )
t = iraf.IrafTaskFactory(taskname="combine_flats", value=parfile, function=combine_flats)

def normalise_flats():

    """
    Normalise flats

    Dispersion axis (1=along lines, 2=along columns, 3=along z) (1:3) (1): 2
    """

    #print 'In directory ' + flatdir
    #print 'Normalising combinined flats...'

    if os.path.exists(os.path.join(data_path, "nFlat.fits")):
        os.remove(os.path.join(data_path, "nFlat.fits"))
        #print 'Removing file ' + os.path.join(flatdir,'nFlat.fits')

    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)

    iraf.response.setParam('calibration', os.path.join(data_path, "Flat.fits"))
    iraf.response.setParam('normalization', os.path.join(data_path, "Flat.fits"))
    iraf.response.setParam('response', os.path.join(data_path, "nFlat.fits"))
    iraf.response.setParam('low_reject', 3.)
    iraf.response.setParam('high_reject', 3.)
    iraf.response.setParam('function', 'spline3')
    iraf.response.setParam('order',40)

    iraf.response()

    return None

parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/normalise_flats.par') )
t = iraf.IrafTaskFactory(taskname="normalise_flats", value=parfile, function=normalise_flats)

def flat_correction():

    """
    Flat field correction
    """

    if os.path.exists('input_flat.csv'):
        os.remove('input_flat.csv')
    if os.path.exists('output_flat.csv'):
        os.remove('output_flat.csv')

    for filename in glob.glob(os.path.join(data_path, night+"*.fits")):
        #bias corrected files (do all even if you don't need to)
        myfile = open('input_flat.csv','a') 
        myfile.write(filename + '\n') 
        myfile.close()
        #flat and bias corrected and trimmed files
        myfile = open('output_flat.csv','a') 
        myfile.write(os.path.join(data_path,'f_'+night+filename.lstrip(data_path)) + '\n') 
        myfile.close()

 
    #settings
    iraf.ccdproc.setParam('images','@input_flat.csv')
    iraf.ccdproc.setParam('overscan','no')
    iraf.ccdproc.setParam('trim','no')
    iraf.ccdproc.setParam('zerocor','no')
    iraf.ccdproc.setParam('flatcor','yes')
    iraf.ccdproc.setParam('readaxis','line')
    iraf.ccdproc.setParam('biassec',"["+"305:360,5:2045"+"]")
    iraf.ccdproc.setParam('trimsec',"["+"5:300,5:2045"+"]")
    iraf.ccdproc.setParam('zero',os.path.join(data_path, "bias"))
    iraf.ccdproc.setParam('flat',os.path.join(data_path, "nFlat"))
    iraf.ccdproc.setParam('interactive','no')
    iraf.ccdproc.setParam('output', '@output_flat.csv')

    iraf.ccdproc()

    return None

parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/flat_correction.par') )
t = iraf.IrafTaskFactory(taskname="flat_correction", value=parfile, function=flat_correction)


def lacos(file_input, file_output, mask):
    '''
    removes the cosmic rays above a threshold of eg. 10 sigma (sigclip)
    '''
    iraf.noao(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.apextract(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.stsdas(_doprint=0)

    iraf.images(_doprint=0)
    iraf.immatch(_doprint=0)


    iraf.lacos_im.setParam('input',file_input)
    iraf.lacos_im.setParam('output',file_output)
    iraf.lacos_im.setParam('outmask',mask)
    iraf.lacos_im.setParam('readn',7.9)
    iraf.lacos_im.setParam('gain',2.7)
    iraf.lacos_im.setParam('sigclip',14)
    iraf.lacos_im.setParam('niter',4)
    iraf.lacos_im.setParam('verbose','no')

    iraf.lacos_im()

    return None


parfile = iraf.osfn( os.path.join(homedir,'/Users/iuliasimion/iraf/pars/lacos.par') )
t = iraf.IrafTaskFactory(taskname="lacos", value=parfile, function=lacos)

#inputf = 'f_n2.0068.fits'
#outputf = 'cr_f_n2.0068.fits'
#maskf = 'mask_n2.0068'

#lacos(inputf , outputf)

def cr_removal():
    flat_frames = []
    bias_frames = []
    object_files = []
    arc_frames = []
    targets = []


    if os.path.exists(os.path.join(data_path, 'lacos.csv')):
        os.remove(os.path.join(data_path, 'lacos.csv'))
        print "removing ", 'lacos'
    for filename in glob.glob(os.path.join(data_path, "cr_*")):
        if os.path.exists(filename):
            os.remove(filename)
    for filename in glob.glob(os.path.join(data_path, "mask_*")):
        if os.path.exists(filename):
            os.remove(filename)

    #create a list input_cr and output_cr and then iterate on those
    for filename in glob.glob(os.path.join(data_path, 'f_'+night+"*.fits")):
        hdr = fits.getheader(filename, 0)    
        obje = hdr['OBJECT']
        try:
            if obje == "FLAT":
                flat_frames.append(filename)
            elif obje == "ZERO":
                bias_frames.append(filename)
            elif obje == "COMP":
                print 'comp'
                arc_frames.append(filename)
            else:
                print 'here'
                object_files.append(night+filename.lstrip(data_path))
                targets.append(obje)
                input_f = os.path.join(data_path,filename.lstrip(data_path))
                output_f = os.path.join(data_path, 'cr_'+filename.lstrip(data_path))
                mask_f = os.path.join(data_path,'mask_'+filename.lstrip(data_path))
                myfile = open(os.path.join(data_path, 'lacos.csv'), 'a') 
                myfile.write(input_f + '\t' + output_f + '\t' + mask_f + '\n') 
                myfile.close()
                print input_f, output_f, mask_f
                #lacos(input_f, output_f, mask_f)
        except:
            print "finished reading"

    for lines in np.genfromtxt(os.path.join(data_path, 'lacos.csv'), dtype = str):
        print lines
        lacos(lines[0], lines[1], lines[2])
    #clean up after lacos_im

    return None


parfile = iraf.osfn( os.path.join(homedir,'/Users/iuliasimion/iraf/pars/cr_removal.par') )
t = iraf.IrafTaskFactory(taskname="cr_removal", value=parfile, function=cr_removal)


def image_shift():

    """
    uses input_shift.csv, output_shift.csv and shift.csv created with shift_files()
    automatically shift images with multiple exposures to match the first exposure 
    if you want to check the shift: ds9 s_f_n2.0022.fits s_f_n2.0023.fits s_f_n2.0024.fits
    """
    #shift the flat corrected files( they ll be use when combining 3 exposures)
    #shift_files(cr = 'no')
    flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()
    target_names = f7(targets)

    #shift the cr cleaned images (to be used for 1 and 2 exposures)
    shift_files()

    #check if the shifted images exist
    #for filename in glob.glob(os.path.join(data_path, "s_*.fits")):
    #    if os.path.exists(filename):
    #        os.remove(filename)
    #        print 'removing'

    iraf.images(_doprint=0)
    iraf.imgeom(_doprint=0)
    iraf.imshift.setParam('input', '@'+os.path.join(data_path, 'input_shift.csv'))
    iraf.imshift.setParam('output', '@'+os.path.join(data_path_reduced, 'output_shift.csv'))
    iraf.imshift.setParam('shifts_file', os.path.join(data_path,'shift.csv') )

    iraf.imshift()

    return None

parfile = iraf.osfn( os.path.join(homedir,'/Users/iuliasimion/iraf/pars/image_shift.par') )
t = iraf.IrafTaskFactory(taskname="image_shift", value=parfile, function=image_shift)



def image_combine(targets_to_combine, output_image):

    """
    Combine all images. Need to run image_shift before.

    Note: this only gets rid of the cosmic rays if you have 3 exposeures or more
    use lacos_im 
    """
    
    iraf.images(_doprint=0)
    iraf.immatch(_doprint=0)

    iraf.imcombine.setParam('input','@'+targets_to_combine)
    iraf.imcombine.setParam('output',output_image)
    iraf.imcombine.setParam('combine','median')
    iraf.imcombine.setParam('reject','crreject')
    iraf.imcombine.setParam('lsigma',3.0)
    iraf.imcombine.setParam('hsigma',3.0)
    iraf.imcombine.setParam('rdnoise',7.9)
    iraf.imcombine.setParam('gain',2.7)
    iraf.imcombine.setParam('scale','mode')
    iraf.imcombine()
    print 'finished combining', output_image
    return None

#imcombine s_f_n2.0041, s_f_n2.0042 test reject=ccreject rdnoise=7.9, gain=2.7
parfile = iraf.osfn( os.path.join(homedir,'/Users/iuliasimion/iraf/pars/image_combine.par') )
t = iraf.IrafTaskFactory(taskname="image_combine", value=parfile, function=image_combine)

def combine_images():
    flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()
    target_names = f7(targets)

    '''
    clean the junk
    '''
    for j in range(len(target_names)):
        if os.path.exists(os.path.join(data_path_reduced, target_names[j]+'_'+str(targets.count(target_names[j]))+'_med.fits')):
            os.remove(os.path.join(data_path_reduced, target_names[j]+'_'+str(targets.count(target_names[j]))+'_med.fits'))
            print "removing ",os.path.join(data_path_reduced, target_names[j]+'_'+str(targets.count(target_names[j]))+'_med.fits')
 
    '''
    finished cleaning
    '''

    '''
    input files imcombine
    '''
    for j in range(0,len(target_names)):

        if os.path.exists(os.path.join(data_path, "input_combine.csv")):
            os.remove(os.path.join(data_path, "input_combine.csv"))
            print "removing input_combine.csv"

        if (targets.count(target_names[j])) > 1:
            ind = [i for i,x in enumerate(targets) if x == target_names[j]]
            print ind
            for index in ind:
                #use shifted images
                print target_names[j], object_files[index]
                filen = 's_'+'cr_'+'f_'+object_files[index]+'s'
                myfile = open(os.path.join(data_path, 'input_combine.csv'), 'a') 
                myfile.write(os.path.join(data_path, filen) + '\n') 
                myfile.close()
        else:
            #no shift
            index = [i for i,x in enumerate(targets) if x == target_names[j]]
            filen = 'cr_'+'f_'+object_files[np.array(index)]+'s'
            myfile = open(os.path.join(data_path, 'input_combine.csv'), 'a') 
            myfile.write(os.path.join(data_path, filen) + '\n') 
            myfile.close()

        #run imcombine
        image_combine(os.path.join(data_path, 'input_combine.csv'), os.path.join(data_path_reduced, target_names[j]+'_'+str(targets.count(target_names[j]))+'_med.fits'))

    return None


parfile = iraf.osfn( os.path.join(homedir,'/Users/iuliasimion/iraf/pars/combine_images.par') )
t = iraf.IrafTaskFactory(taskname="combine_images", value=parfile, function=combine_images)

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

def extract_spectrum(input_target, output_target, input_arc, output_arc, no_exp):

    """
    Extract spectrum

    Must be in target directory.

    """

    if os.path.exists( os.path.join(data_path_reduced, input_target +'_mo'+'.fits') ):
        os.remove( os.path.join( data_path_reduced, input_target +'_mo'+'.fits') )
        print 'Removing file ' + os.path.join( data_path_reduced, input_target +'_mo'+'.fits')

    #modify readnoise 
    readnoise  = 7.9/np.sqrt(no_exp)
    print 'no of exposures, readnoise =' , no_exp, readnoise

    #modify input frame
    iraf.unlearn('imarith')

    iraf.imarith.setParam('operand1', input_target )
    iraf.imarith.setParam('operand2', no_exp) #modified input target
    iraf.imarith.setParam('op','*')
    iraf.imarith.setParam('result', input_target +'_mo')

    iraf.imarith()

    iraf.noao(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.apextract(_doprint=0)

    #filenames_to_remove = ['targets_to_extract.csv', 'targets_extracted.csv', 'arcs_to_extract.csv', 'arcs_extracted.csv']

    interactive = 'yes'
    trace = 'yes'
    # Doesn't seem to work if I give it absolute path to input!
    iraf.apall.setParam('input',input_target+'_mo') # List of input images #CHANGE
    iraf.apall.setParam('output',output_target) # List of output spectra #CHANGE
    iraf.apall.setParam('apertures', '') # Apertures
    iraf.apall.setParam('format','multispec') # Extracted spectra format
    iraf.apall.setParam('referen','') # List of aperture reference images
    iraf.apall.setParam('profile','') # List of aperture profile images

    iraf.apall.setParam('interactive',interactive) # Run task interactively?
    iraf.apall.setParam('find','yes') # Find apertures?
    iraf.apall.setParam('recenter','yes') # Recenter apertures?
    iraf.apall.setParam('resize','yes') # Resize apertures?
    iraf.apall.setParam('edit','yes') # Edit apertures?
    iraf.apall.setParam('trace','yes') # Trace apertures?
    iraf.apall.setParam('fittrac',interactive) # Fit the traced points interactively?
    iraf.apall.setParam('extract','yes') # Extract spectra?
    iraf.apall.setParam('extras','yes') # Extract sky, sigma, etc.?
    iraf.apall.setParam('review',interactive) # Review extractions?

    iraf.apall.setParam('line','INDEF') # Dispersion line
    iraf.apall.setParam('nsum',20) # Number of dispersion lines to sum or median

                                # DEFAULT APERTURE PARAMETERS

    iraf.apall.setParam('lower',-4.) # Lower aperture limit relative to center
    iraf.apall.setParam('upper',4.) # Upper aperture limit relative to center
    iraf.apall.setParam('apidtab','') # Aperture ID table (optional)

                                # DEFAULT BACKGROUND PARAMETERS

    iraf.apall.setParam('b_funct','chebyshev') # Background function
    iraf.apall.setParam('b_order',1) # Background function order
    iraf.apall.setParam('b_sampl','-80:-10,10:80') # Background sample regions
    iraf.apall.setParam('b_naver',-8) # Background average or median
    iraf.apall.setParam('b_niter',2) # Background rejection iterations
    iraf.apall.setParam('b_low_r',3.) # Background lower rejection sigma
    iraf.apall.setParam('b_high_',3.) # Background upper rejection sigma
    iraf.apall.setParam('b_grow',0.) # Background rejection growing radius

                                # APERTURE CENTERING PARAMETERS

    iraf.apall.setParam('width',8.) # Profile centering width
    iraf.apall.setParam('radius',10.) # Profile centering radius
    iraf.apall.setParam('thresho',0.) # Detection threshold for profile centering

                                # AUTOMATIC FINDING AND ORDERING PARAMETERS

    iraf.apall.setParam('nfind',1) # Number of apertures to be found automatically
    iraf.apall.setParam('minsep',5.) # Minimum separation between spectra
    iraf.apall.setParam('maxsep',1000.) # Maximum separation between spectra
    iraf.apall.setParam('order','increasing') # Order of apertures

                                # RECENTERING PARAMETERS

    iraf.apall.setParam('aprecen','') # Apertures for recentering calculation
    iraf.apall.setParam('npeaks','INDEF') # Select brightest peaks
    iraf.apall.setParam('shift','yes') # Use average shift instead of recentering?

                                # RESIZING PARAMETERS

    iraf.apall.setParam('llimit','INDEF') # Lower aperture limit relative to center
    iraf.apall.setParam('ulimit','INDEF') # Upper aperture limit relative to center
    iraf.apall.setParam('ylevel',0.2) # Fraction of peak or intensity for automatic widt
    iraf.apall.setParam('peak','yes') # Is ylevel a fraction of the peak?
    iraf.apall.setParam('bkg','yes') # Subtract background in automatic width?
    iraf.apall.setParam('r_grow',0.) # Grow limits by this factor
    iraf.apall.setParam('avglimi','no') # Average limits over all apertures?

                                # TRACING PARAMETERS

    iraf.apall.setParam('t_nsum',10) # Number of dispersion lines to sum
    iraf.apall.setParam('t_step', 10) # Tracing step
    iraf.apall.setParam('t_nlost',3) # Number of consecutive times profile is lost befo
    iraf.apall.setParam('t_funct','spline3') # Trace fitting function
    iraf.apall.setParam('t_order',1) # Trace fitting function order
    iraf.apall.setParam('t_sampl','*') # Trace sample regions
    iraf.apall.setParam('t_naver',1) # Trace average or median
    iraf.apall.setParam('t_niter',0) # Trace rejection iterations
    iraf.apall.setParam('t_low_r',3.) # Trace lower rejection sigma
    iraf.apall.setParam('t_high_',3.) # Trace upper rejection sigma
    iraf.apall.setParam('t_grow',0.) # Trace rejection growing radius

                                # EXTRACTION PARAMETERS

    iraf.apall.setParam('backgro','median') # Background to subtract
    iraf.apall.setParam('skybox',1) # Box car smoothing length for sky
    iraf.apall.setParam('weights','variance') # Extraction weights (none|variance)
    iraf.apall.setParam('pfit','fit1d') # Profile fitting type (fit1d|fit2d)
    iraf.apall.setParam('clean','yes') # Detect and replace bad pixels?
    #iraf.apall.setParam('saturat',300000.) # Saturation level
    iraf.apall.setParam('readnoi',readnoise) # Read out noise sigma (photons)
    iraf.apall.setParam('gain',2.7) # Photon gain (photons/data number)
    iraf.apall.setParam('lsigma',4.) # Lower rejection threshold
    iraf.apall.setParam('usigma',4.) # Upper rejection threshold
    iraf.apall.setParam('nsubaps',1) # Number of subapertures per aperture
    iraf.apall.setParam('mode','q') # h = hidden, q = query, l = learn

    iraf.apall()



    # Now extract arc through same aperture for wavelength calibration
    
    print '\n' '\n' '\n'
    print 'Extracting Arc through same aperture...'

    iraf.apall.setParam('input', input_arc)
    iraf.apall.setParam('output', output_arc)
    iraf.apall.setParam('references', input_target )
    iraf.apall.setParam('recenter','no')
    iraf.apall.setParam('trace','no')
    iraf.apall.setParam('background','no')
    iraf.apall.setParam('interactive','no')
    
    iraf.apall()
    
    return None


parfile = iraf.osfn( os.path.join(homedir,'/Users/iuliasimion/iraf/pars/extract_spectrum.par') )
t = iraf.IrafTaskFactory(taskname="extract_spectrum", value=parfile, function=extract_spectrum)

def run_apall():
    '''
    run this in the data_path_reduced
    '''

    flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()
    target_names = f7(targets)
    for i in range(len(target_names)):
        input_targetf = target_names[i] +'_'+str(targets.count(target_names[i]))+'_med'
        output_targetf = target_names[i] +'_'+str(targets.count(target_names[i]))+'_med'+'.1d'
        input_arcf = arc_frames[i] 
        output_arcf = 'arc_' +target_names[i] 
        print 'input_target, output_target, input_arcf, output_arcf', input_targetf, output_targetf, input_arcf, output_arcf

        filenames_to_remove = [output_targetf, output_arcf]
        for item in filenames_to_remove:
            if os.path.exists( os.path.join(data_path_reduced, item+'.fits') ):
                os.remove( os.path.join( data_path_reduced, item+'.fits') )
                print 'Removing file ' + os.path.join( data_path_reduced, item)

        extract_spectrum(input_targetf, output_targetf, input_arcf, output_arcf,targets.count(target_names[i]) )
    return None

parfile = iraf.osfn( os.path.join(homedir,'/Users/iuliasimion/iraf/pars/run_apall.par') )
t = iraf.IrafTaskFactory(taskname="run_apall", value=parfile, function=run_apall)

def run_apall_1exp():
    '''
    run this in the data_path
    here you extract individual exposures and not the combined spectrum
    '''
    print 'here'
    flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()
    target_names = f7(targets)
    print 'here2'

    for filename in object_files:
        filen = 's_cr_f_'+filename.rstrip('.fit')
        print filen
        #hdr = fits.getheader(filename, 0)    
        #obje = hdr['OBJECT']
   
        input_targetf = 's_cr_f_'+filename.rstrip('.fit')
        output_targetf = 's_cr_f_'+filename.rstrip('.fit')+'.1d'

        filenames_to_remove = [output_targetf]
        for item in filenames_to_remove:
            if os.path.exists( os.path.join(data_path_reduced, item+'.fits') ):
                os.remove( os.path.join( data_path_reduced, item+'.fits') )
                print 'Removing file ' + os.path.join( data_path_reduced, item+'.fits')

        print input_targetf, output_targetf
        #dont extract the arc so it doesnt matter inputarcf and outputarcf, comment the arc extraction lines in apall
        extract_spectrum(input_targetf, output_targetf, filen, filen,1)
        

    return None

parfile = iraf.osfn( os.path.join(homedir,'/Users/iuliasimion/iraf/pars/run_apall_1exp.par') )
t = iraf.IrafTaskFactory(taskname="run_apall_1exp", value=parfile, function=run_apall_1exp)

def wavelength_calibration(target, arc):

    """
    Does wavelength calibration.

    Writes every fit to database so make sure it's using the correct one.

    This needs to be run in object directory for database

    """

    print 'Target directory is ' + data_path_reduced
    print 'Doing wavelength calibration...'

    if os.getcwd() != (data_path_reduced):

        print 'Warning: current working directory must be target directory!'

        return None

    iraf.noao(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.apextract(_doprint=0)
    iraf.longslit(_doprint=0)

    iraf.unlearn('identify')

    iraf.identify.setParam('images',arc)
    iraf.identify.setParam('coordli','/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/lines.dat')
    iraf.identify.setParam('niterat',1)
    iraf.identify.setParam('function','spline3')
    iraf.identify.setParam('order',3)
    iraf.identify.setParam('zwidth',200.0) #  Zoom graph width in user units
    iraf.identify.setParam('database','database')

    #iraf.identify(Stdout=os.path.join(data_path_reduced+'/identify', target+"_identify.stdout"))
    iraf.identify()

    # Update fits header: pair the arc with the image.

    print '\n' '\n' '\n'
    print 'Updating fits header...'
    
    iraf.hedit.setParam('images',target)
    iraf.hedit.setParam('fields','REFSPEC1')
    iraf.hedit.setParam('value',arc) # should be wavelength calibrated?
    iraf.hedit.setParam('add','yes')
    iraf.hedit.setParam('verify','yes')
    iraf.hedit.setParam('show','yes')

    iraf.hedit()
    
    return None

parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/wavelength_calibration.par') )
t = iraf.IrafTaskFactory(taskname="wavelength_calibration", value=parfile, function=wavelength_calibration)

def reidentify_lines(arc, reference_arc, target, head = 'no'):
    print 'Target directory is ' + data_path_reduced
    print 'Doing wavelength calibration...'
    iraf.noao(_doprint=0)
    iraf.onedspec(_doprint=0)

    if os.getcwd() != (data_path_reduced):

        print 'Warning: current working directory must be target directory!'

    if head == 'no':
        iraf.hedit.setParam('images',arc)
        iraf.hedit.setParam('fields','REFSPEC1')
        iraf.hedit.setParam('value',arc) # should be wavelength calibrated?
        iraf.hedit.setParam('add','no')
        iraf.hedit.setParam('delete','yes')
        iraf.hedit()
        print 'deleting REFSPEC'
 
    iraf.unlearn('reidentify')

    iraf.reidentify.setParam('reference',reference_arc)
    iraf.reidentify.setParam('images',arc)
    iraf.reidentify.setParam('interactive','YES')
    iraf.reidentify.setParam('newaps','no')
    iraf.reidentify.setParam('override','no')
    iraf.reidentify.setParam('refit','yes') #  Zoom graph width in user units
    iraf.reidentify.setParam('trace','no')
    iraf.reidentify.setParam('coordlist','/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/lines.dat')
    iraf.reidentify.setParam('answer','yes')

    iraf.reidentify.setParam('addfeatures','yes')
    iraf.reidentify.setParam('verbose','yes')

    iraf.reidentify()

    print '\n' '\n' '\n'
    print 'Updating fits header...'

    if head == 'yes':
        iraf.hedit.setParam('images',target)
        iraf.hedit.setParam('fields','REFSPEC1')
        iraf.hedit.setParam('value',arc) # should be wavelength calibrated?
        iraf.hedit.setParam('add','yes')
        iraf.hedit.setParam('verify','yes')
        iraf.hedit.setParam('show','yes')

        iraf.hedit()
    return None


parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/reidentify_lines.par') )
t = iraf.IrafTaskFactory(taskname="reidentify_lines", value=parfile, function=reidentify_lines)


def wavelength_solution(input_target, output_target , islog = 'no'):

    """
    Applies wavelength solution using dispcor. Must be run in target directory.
    """
    print 'Target directory is ' + data_path_reduced
    print 'Doing wavelength calibration...'

    if os.getcwd() != (data_path_reduced):

        print 'Warning: current working directory must be target directory!'

    if os.path.exists(output_target):
        os.remove(output_target)
        print 'removing', output_target

    iraf.unlearn('dispcor')

    iraf.dispcor.setParam('input',input_target)
    iraf.dispcor.setParam('output',output_target)
    iraf.dispcor.setParam('confirm','yes')
    iraf.dispcor.setParam('log',islog)
    iraf.dispcor.setParam('w1',3618.339619574508) # Starting wavelength
    iraf.dispcor.setParam('w2',7704.026691620338) # Ending wavelength
    iraf.dispcor.setParam('nw',2041) # Number of output pixels
    iraf.dispcor.setParam('dw','INDEF') # Wavelength interval per pixel

    iraf.dispcor()

    return None

parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/wavelength_solution.par') )
t = iraf.IrafTaskFactory(taskname="wavelength_solution", value=parfile, function=wavelength_solution)

def run_identify():
    islog = 'no'
    flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()
    run = 2
    target_names = f7(targets)
    #for i in range(len(target_names)):
    if night == 'n2':
        i = 0
    input_arcf = 'arc_' +target_names[i]  #extracted arc with same aperture as target
    input_targetf = target_names[i] +'_'+str(targets.count(target_names[i]))+'_med.1d.fits' #extracted target spectrum
    if islog == 'no':
        output_targetf = target_names[i] +'_'+str(targets.count(target_names[i]))+'_med.cal.fits'
    elif islog == 'yes':
        output_targetf = target_names[i] +'_'+str(targets.count(target_names[i]))+'_med.log.cal.fits'

    print 'input arc, input target:',night, input_arcf, input_targetf
    #for night 1 identify arc 0034
    #to the wavelength calibration on the arcfile and update the header of the target
    wavelength_calibration(input_targetf, input_arcf)
    wavelength_solution(input_targetf, output_targetf)

    for i in range(1,len(target_names)):
        input_arcf = 'arc_' +target_names[i]  #extracted arc with same aperture as target
        input_targetf = target_names[i] +'_'+str(targets.count(target_names[i]))+'_med.1d.fits' #extracted target spectrum
        if islog == 'no':
            output_targetf = target_names[i] +'_'+str(targets.count(target_names[i]))+'_med.cal.fits'
        elif islog == 'yes':
            output_targetf = target_names[i] +'_'+str(targets.count(target_names[i]))+'_med.log.cal.fits'
        print 'input arc, input target:',input_arcf, input_targetf
        #i dont know why this is
        if run == 1:
            reidentify_lines('arc_' +target_names[i], 'arc_' +target_names[0], input_targetf, head = 'yes' )
        else:
            wavelength_calibration(input_targetf, input_arcf)
        wavelength_solution(input_targetf, output_targetf)
    return None

parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/run_identify.par') )
t = iraf.IrafTaskFactory(taskname="run_identify", value=parfile, function=run_identify)



def run_identify_1exp():
    '''
    copy .1d.fits in data_path_reduced and run this there 
    '''

    flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()
    run = 2

    for filename in object_files[1:]:
        input_targetf = 's_cr_f_'+filename.rstrip('.fit')+'.1d.fits'
        hdr = fits.getheader(input_targetf, 0)    
        obje = hdr['OBJECT']

        input_arcf = 'arc_' +obje  #extracted arc with same aperture as target
        output_targetf = 's_cr_f_'+filename.rstrip('.fit')+'.1d.cal.fits'
      
        print 'input arc, input target:',input_arcf, input_targetf
        #i dont know why this is
        if run == 1:
            reidentify_lines(input_arcf, 'arc_' +'HAC156', input_targetf, head = 'yes' )
        else:
            wavelength_calibration(input_targetf, input_arcf)
        wavelength_solution(input_targetf, output_targetf)
    return None

parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/run_identify_1exp.par') )
t = iraf.IrafTaskFactory(taskname="run_identify_1exp", value=parfile, function=run_identify_1exp)



def correct_rv(hdr, target_to_rvcorrect, target):
    '''
    estimate shift in wavelength due to the Earth's rotation
    the correction should be less than 30km/s
    '''

    date = hdr['DATE-OBS']
    year = date[0:4]
    month = date[5:7]
    day = date[8:10]
    time = hdr['TIME-OBS']
    ra = hdr['RA']
    dec=hdr['DEC']
    print 'correcting rv',date, time, ra, dec

    if os.path.exists('vh_'+target_to_rvcorrect):
        os.remove('vh_'+target_to_rvcorrect)
        print 'removing', 'vh_'+target_to_rvcorrect

    iraf.rv(_doprint=0)

    iraf.unlearn('rvcorrect')
    iraf.rvcorrect.setParam('files','')
    iraf.rvcorrect.setParam('images','')
    iraf.rvcorrect.setParam('header','yes')
    iraf.rvcorrect.setParam('input','no')
    iraf.rvcorrect.setParam('imupdate','yes')

    iraf.rvcorrect.setParam('observatory','mdm')
    iraf.rvcorrect.setParam('year',year)
    iraf.rvcorrect.setParam('month',month)
    iraf.rvcorrect.setParam('day',day)
    iraf.rvcorrect.setParam('ut',time)

    iraf.rvcorrect.setParam('ra',ra)
    iraf.rvcorrect.setParam('dec',dec)


    iraf.rvcorrect(Stdout=os.path.join(data_path_reduced+'/vhelio', target+"_rvcorrect.stdout"))
    rvfile = np.genfromtxt(os.path.join(data_path_reduced+'/vhelio', target+"_rvcorrect.stdout")) 
    vhelio = rvfile[2]*(-1.)

    #modify header
    iraf.unlearn('hedit')
    iraf.hedit.setParam('images',target_to_rvcorrect)
    iraf.hedit.setParam('fields','VHELIO')
    iraf.hedit.setParam('value',rvfile[2]) # should be wavelength calibrated?
    iraf.hedit.setParam('add','yes')
    iraf.hedit.setParam('verify','yes')
    iraf.hedit.setParam('show','yes')
    iraf.hedit()
  
    #shift image:
    iraf.unlearn('dopcor')
    iraf.dopcor.setParam('input',target_to_rvcorrect)
    iraf.dopcor.setParam('output','vh_'+target_to_rvcorrect)
    iraf.dopcor.setParam('redshift',vhelio) 
    iraf.dopcor.setParam('isvelocity','yes')
    iraf.dopcor.setParam('add','no')
    iraf.dopcor.setParam('dispersion','yes')
    iraf.dopcor.setParam('flux','no')
    iraf.dopcor.setParam('verbose','yes')
    iraf.dopcor()
    #dopcor Xari Xari.fin -27.42 isveloc+
    return None
#2456900.62700     0.00   -13.40    -1.93      0.247   -0.013  -13.637   11.473
#2456900.62700     0.00   -13.40    -1.93      0.247   -0.013  -13.637   11.473
#2456900.65405     0.00    -7.85     3.57      0.246   -0.012   -8.083   11.424
#2456900.65405     0.00    -7.85     3.57      0.246   -0.012   -8.083   11.424


parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/correct_rv.par') )
t = iraf.IrafTaskFactory(taskname="correct_rv", value=parfile, function=correct_rv)


def run_rvcorrect():
    flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()

    target_names = f7(targets)
 
    for i in range(0,len(target_names)):
        target_to_rvcorrect = target_names[i] +'_'+str(targets.count(target_names[i]))+'_med.log.cal.fits'
        if targets.count(target_names[i]) > 1:
            input_target = target_names[i] +'_'+str(targets.count(target_names[i]))+'_med.log.cal.fits'
            hdr = fits.getheader(input_target, 0) 
            correct_rv(hdr, target_to_rvcorrect, target_names[i])   
        if targets.count(target_names[i]) == 1:
            filename = object_files[i]
            input_target = os.path.join(data_path, filename+'s')
            hdr = fits.getheader(input_target, 0) 
            correct_rv(hdr, target_to_rvcorrect, target_names[i])
    return None


parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/run_rvcorrect.par') )
t = iraf.IrafTaskFactory(taskname="run_rvcorrect", value=parfile, function=run_rvcorrect)

def run_rvcorrect_1exp():
    flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()

    for filename in object_files[1:]:

        target_to_rvcorrect = 's_cr_f_'+filename.rstrip('.fit')+'.1d.cal.fits'        
        hdr = fits.getheader(target_to_rvcorrect, 0) 
        obje = hdr['OBJECT']
        correct_rv(hdr, target_to_rvcorrect, obje)
    #s_f_n3.0022.fits.1d.cal.fits - > vh_s_f_n3.0022.fits.1d.cal.fits
    return None

parfile = iraf.osfn(os.path.join(homedir,'/Users/iuliasimion/iraf/pars/run_rvcorrect_1exp.par') )
t = iraf.IrafTaskFactory(taskname="run_rvcorrect_1exp", value=parfile, function=run_rvcorrect_1exp)
