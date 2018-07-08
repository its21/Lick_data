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
import Fourier as F
from radec2deg import string2hours
from radec2deg import string2deg

'''
final version
- calculates the velocity and the errors
'''

folder_path = "/Users/iuliasimion/Observing2014/hac_spectra/allyson/data_reduction_example/"

#choose the path for your main folder:
data_path = os.path.expanduser(folder_path+"data/")
data_path_reduced = os.path.expanduser(folder_path+"data_reduced/")
#templates_path = os.path.expanduser(folder_path+"T3500-7250/") #for Allyson
templates_path = os.path.expanduser("/Users/iuliasimion/Observing2014/T3500-7250/")
plot_path = os.path.expanduser(folder_path+"plots/")
#results_path = os.path.expanduser(folder_path+"results/") #for Allyson
results_path = os.path.expanduser("/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/results/") 
#I take the results from here because it's where I save everything when I run other grids

c =  299792.458

'''
 remember to change the night in resolution_logw.py
'''

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

def chi_func(rv, w, f, dv, target, T, G, M,  N = 15):
	'''
	w, f, dv: wavelength and target spectrum, variance spectrum
	rv = the best fit velocity
	T,G,M are the temperature, logg, metallicity of the best fit template 
	N = degree of the polynomial that multiplies the template
	target = target name
	for this template you find the chi-squared distribution and from it you find the error 
	use this function to calculate the template fitting error on the velocity,
	from the chi squared distribution
	'''

	#feature scaling (-0.5<tw<0.5) : important because w ~ 1e3 and you have exponentials in the poly
	ww = w[((w<(np.log10(6300-15))) | (w>np.log10((6300+15)))) & ((w<np.log10((5577-9))) | (w>np.log10((5577+9)))) & ((w < np.log10((6868-9))) | (w > np.log10((6868+12)))) &  ((w < np.log10((6363-5))) | (w > np.log10((6363+5)))) ]
	ff = f[((w<(np.log10(6300-15))) | (w>np.log10((6300+15)))) & ((w<np.log10((5577-9))) | (w>np.log10((5577+9)))) & ((w < np.log10((6868-9))) | (w > np.log10((6868+12)))) &  ((w < np.log10((6363-5))) | (w > np.log10((6363+5)))) ]
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
	drv = 1
	drvs = np.arange(rvmin, rvmax, drv)

	#choose your templates grid
	#Read the best fit template
	if M < 10:
		try:
			filename = os.path.join(templates_path+'T_0'+str(T)+"/", "T0"+str(T)+'G'+str(G)+"M0"+str(M)+"V000K2ANWNVD01F.ASC")
			templ_file = np.loadtxt(filename)
		except:
			filename = os.path.join(templates_path+'T_0'+str(T)+"/", "T0"+str(T)+'G'+str(G)+"M0"+str(M)+"V000K2SNWNVD01F.ASC")
			templ_file = np.loadtxt(filename)
	else:
		try:
			filename = os.path.join(templates_path+'T_0'+str(T)+"/", "T0"+str(T)+'G'+str(G)+"M"+str(M)+"V000K2ANWNVD01F.ASC")
			templ_file = np.loadtxt(filename)
		except:
			filename = os.path.join(templates_path+'T_0'+str(T)+"/", "T0"+str(T)+'G'+str(G)+"M"+str(M)+"V000K2SNWNVD01F.ASC")
			templ_file = np.loadtxt(filename)

	fixpix = 0.00011
	logw_regular = np.arange(min(np.log10(tw)), max(np.log10(tw)), fixpix)
	#get your template on a regular grid
	templ = sc.interpolate.interp1d(np.log10(tw), templ_file)
	win = signal.gaussian(5, std = mean_sigma/fixpix)
	tf = signal.convolve(templ(logw_regular), win, mode='same') / sum(win)

	fi = sc.interpolate.interp1d(logw_regular+np.log10(1.0 + rv/c), tf)

	AT = np.transpose(A*fi(ww)/dvv)
	coeff_poly, residuals = scipy.linalg.lstsq(AT, ff/dvv)[0] , scipy.linalg.lstsq(AT, ff/dvv)[1]

	rank  = scipy.linalg.lstsq(AT, ff/dvv)[2]
	singular  = scipy.linalg.lstsq(AT, ff/dvv)[3]
	redchi = residuals/(len(ww)-N)
	chi = residuals

	return chi


def func(vel):
	wd, fd, dvd, target_name,T_best, logg_best, met_best = cPickle.load(open(os.path.join(data_path, "tmp.dat")))
	return chi_func(vel, wd, fd, dvd, target_name,int(T_best), int(logg_best), int(met_best), N = 15)


def best_templ(w, f, dv, target, T, G, M,  N = 10):
	#feature scaling (-0.5<tw<0.5) : important because w ~ 1e3 and you have exponentials in the poly

	ww = w[((w<(np.log10(6300-15))) | (w>np.log10((6300+15)))) & ((w<np.log10((5577-9))) | (w>np.log10((5577+9)))) & ((w < np.log10((6868-9))) | (w > np.log10((6868+12)))) &  ((w < np.log10((6363-5))) | (w > np.log10((6363+5))))  ]
	ff = f[((w<(np.log10(6300-15))) | (w>np.log10((6300+15)))) & ((w<np.log10((5577-9))) | (w>np.log10((5577+9)))) & ((w < np.log10((6868-9))) | (w > np.log10((6868+12)))) &  ((w < np.log10((6363-5))) | (w > np.log10((6363+5)))) ]
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
	drv =1.
	drvs = np.arange(rvmin, rvmax, drv)

	#choose your templates grid
	best_vel = []

	#Read the template
	if M < 10:
		filename = os.path.join(templates_path+'T_0'+str(T)+"/", "T0"+str(T)+'G'+str(G)+"M0"+str(M)+"V000K2ANWNVD01F.ASC")
	else:
		filename = os.path.join(templates_path+'T_0'+str(T)+"/", "T0"+str(T)+'G'+str(G)+"M"+str(M)+"V000K2ANWNVD01F.ASC")
	#print filename
	templ_file = np.loadtxt(filename)

	#NEED TO CHANGE THIS TO LOG  SCALE
	fixpix = 0.00011
	logw_regular = np.arange(min(np.log10(tw)), max(np.log10(tw)), fixpix)
	#get your template on a regular grid
	templ = sc.interpolate.interp1d(np.log10(tw), templ_file)
	win = signal.gaussian(5, std = mean_sigma/fixpix)
	tf = signal.convolve(templ(logw_regular), win, mode='same') / sum(win)

	redchi = np.zeros(len(drvs))
	chi = np.zeros(len(drvs))
	restrict = ((ww>(np.log10(6562-50))) & (ww<np.log10((6562+50)))) | ((ww>np.log10((4861-50))) & (ww<np.log10((4861+50)))) 

	for i, rv in enumerate(drvs):
		fi = sc.interpolate.interp1d(logw_regular+np.log10(1.0 + rv/c), tf)

		AT = np.transpose(A*fi(ww)/dvv)
		coeff_poly, residuals = scipy.linalg.lstsq(AT, ff/dvv)[0] , scipy.linalg.lstsq(AT, ff/dvv)[1]

		rank  = scipy.linalg.lstsq(AT, ff/dvv)[2]
		singular  = scipy.linalg.lstsq(AT, ff/dvv)[3]
		redchi[i] = residuals/(len(ww[restrict])-N)
		chi[i] = residuals


	velindx = np.argmin(chi)
	vel = drvs[velindx]
	templ = sc.interpolate.interp1d(logw_regular+np.log10(1.0 + vel/c), tf)

	hess = numdifftools.Hessian(func) # Hessian function
	hessmat = 0.5*hess(vel) 
	covariance_matrix = np.matrix(hessmat).I
	error = np.diagonal(abs(covariance_matrix))**.5

	print hess, hessmat, covariance_matrix, error

	#re-calculate poly for the best fit
	AT = np.transpose(A*templ(ww)/dvv)
	coeff_poly, residuals = scipy.linalg.lstsq(AT, ff/dvv)[0] , scipy.linalg.lstsq(AT, ff/dvv)[1]
	coeff_rev =  coeff_poly[::-1]
	p = np.poly1d(coeff_rev)
	#print np.poly1d(p)

	#plot chi distribution for the best fit template
	fig = plt.figure(figsize=(11.5, 7))
	fig.clf()
	plt.subplots_adjust(left=0.11, bottom=0.17, right=0.95, top=0.90, wspace=0.2)
	ax = fig.add_subplot(111)
	plt.title(target)
	plt.plot(drvs, redchi, '.', lw = 4)

	for item in ([ax.title,ax.yaxis.label,ax.xaxis.label]):
	    item.set_fontsize(40)

	for item in (ax.get_xticklabels() + ax.get_yticklabels()):
	    item.set_fontsize(30)
	ax.set_xlabel(r'$v$'+ ' ' +'$(km/s)$')
	ax.set_ylabel(r'$\chi_{red}^{2}$')

	fig.savefig(plot_path+target+'_chi.eps')

	#plot target and best fit template
	fig = plt.figure(figsize=(12.5, 7))
	fig.clf()
	plt.subplots_adjust(left=0.11, bottom=0.13, right=0.95, top=0.94, wspace=0.2)
	ax = fig.add_subplot(111)	
	plt.plot(w, f, '.',label = 'data excluded', lw = 2, color = 'green')
	plt.plot(ww[restrict], ff[restrict], '.',label = 'data fitted', lw = 2, color = 'blue')
	plt.plot(ww, p(dw_s)*templ(ww), '-',label = r'template: $v$ = '+str(int(vel)) +' '+' km/s ', lw = 2, color = 'red')
	plt.plot(ww, ff - p(dw_s)*templ(ww), '-',label = '(data-model)/model', lw = 0.5, color = 'magenta')
	plt.plot(ww, np.zeros(len(ww)), '--',lw = 0.5, color = 'black')

	if target == 'HAC11':
		ax.text(np.log10(6563-100), 120, r'$H_{\alpha}$',ha='left', va='top', fontsize = 30, color = 'blue')
		ax.text(np.log10(4861-80), 120, r'$H_{\beta}$',ha='left', va='top', fontsize = 30, color = 'blue')
		ax.text(np.log10(4340-100), 290, r'$H_{\gamma}$',ha='left', va='top', fontsize = 30, color = 'black')
		ax.text(np.log10(4192-200), 170, r'$H_{\delta}$',ha='left', va='top', fontsize = 30, color = 'black')
		#ax.text(np.log10(6300), 400, '[OI]',ha='left', va='top', fontsize = 30, color = 'black')
		#ax.text(np.log10(5577), 500, r'$[OI]$',ha='left', va='top', fontsize = 30, color = 'black')

	ax.set_title(target)
	ax.set_ylim(-70,1000)
	ax.set_xlabel('log$_{10}(\lambda)$')
	ax.set_ylabel('Intensity (counts)')
	#set size ticks for the main plot
	for item in ([ax.title,ax.yaxis.label,ax.xaxis.label]):
	    item.set_fontsize(22)

	for item in (ax.get_xticklabels() + ax.get_yticklabels()):
	    item.set_fontsize(18)

	plt.legend(loc = 2, fontsize = 20, frameon=False)
	#_template.eps for full spectrum
	
	fig.savefig(plot_path+target+'_template.eps')

	print 'min chi', min(redchi)

	return vel, error, min(redchi)


#read wavelengths template (features)
tw = np.loadtxt(os.path.join(templates_path, "LAMBDA_D1A.dat"))

#templates phase
#read norm_rv templates from Sesar et al. 2012
phase_Talpha, nrvalpha, phase_Tbeta, nrvbeta = [], [], [], []
norm_rv = np.genfromtxt('/Users/iuliasimion/Observing2014/RV_templates/RV_template.txt',dtype=[('mystring', '|S7'), ('myint', '<f8'), ('myfloat', '<f8')] )
for line in norm_rv:
	if line[0] == 'Halpha':
		phase_Talpha.append(line[1])
		nrvalpha.append(line[2])
	if line[0] == 'Hbeta':
		phase_Tbeta.append(line[1])
		nrvbeta.append(line[2])
import scipy as sc
templ_alpha = sc.interpolate.interp1d(phase_Talpha, nrvalpha)
templ_beta = sc.interpolate.interp1d(phase_Tbeta, nrvbeta)

def vGSR(vr,l,b):
    import scipy as sc
    vr=np.array(vr)
    dtor = sc.pi/180.
    vLSR = vr + 10.0*sc.cos(l*dtor)*sc.cos(b*dtor)+7.2*sc.sin(b*dtor)
    #+12.*sc.sin(l*dtor)*sc.cos(b*dtor)
    return vLSR + 225.2*sc.sin(l*dtor)*sc.cos(b*dtor)

'''
calculate the final error on the best fit template for all targets
'''

def plot_target():
	#choose your target
	flat_frames, bias_frames, arc_frames, arc_targets, object_files, targets = files()
	target_names = f7(targets)
	#for night 2 start j at 1 and object_files[1:]
	for j in range(1,len(target_names)):
		fd_new = []
		dvd_new = []
		jd = []
		for filename in object_files[1:-7]:
			filenn = 'vh_s_cr_f_'+filename.rstrip('.fit')+'.1d.cal.fits'  
			target_to_cc =os.path.join(data_path_reduced,filenn)     
			hdr = fits.getheader(target_to_cc, 0) 
			obje = hdr['OBJECT']
			if target_names[j] == obje:
				print target_names[j], filenn
				wd, fd = read_spec(target_to_cc)
				dw, dvd = read_spec(target_to_cc, var = 'yes') #variance
				h = pyfits.open(target_to_cc)[0].header
				jd.append(h['JD'])
				fd_new.append(fd)
				dvd_new.append(dvd)
				wd = np.log10(wd)
				dw = np.log10(dw)
			else:
				continue #print target_names[j],'not found'

		RA_tel = h['RA']
		DEC_tel = h['DEC']
		RA_tel = string2hours(RA_tel)*15.
		DEC_tel = string2deg(DEC_tel)

		dvd_new = np.array(dvd_new)
		fd_new = np.array(fd_new)
		fd_f, dvd_f = 0, 0
		for ss in range(len(fd_new)):
			fd_f += fd_new[ss]#+fd_new[1] +fd_new[2] 
			dvd_f += dvd_new[ss]**2.#+dvd_new[1]**2.+dvd_new[2]**2.
		restrict = ((wd<(np.log10(6562-50))) | (wd>np.log10((6562+50)))) & ((wd<np.log10((4861-50))) | (wd>np.log10((4861+50)))) 
		
		dvd_f[restrict] = 100*dvd_f[restrict]
		dvd_f[dvd_f==0]=1


		savefile = night+'_'+target_names[j]+'_templatefit.txt'
		values = np.loadtxt(os.path.join(results_path, savefile))

		dvd[dvd==0]=1
		chi = values[:,4]
		indx = np.argmin(chi)
		best_template = values[indx,:]

		T_best = best_template[0]
		logg_best = best_template[1]
		met_best = best_template[2]
		vel_best = best_template[3]
		chi_best = best_template[4]
		redchi_best = best_template[5]


		veloc, err, minchi = best_templ(wd, fd_f, np.sqrt(dvd_f), target_names[j],int(T_best), int(logg_best), int(met_best), N = 15)
		err_cal =np.sqrt(err**2.+sigma_vel_cal**2.)	
		name=target_names[j]

		# the light curves are from the CSS website and I saved them for each HAC candidate
		filen = os.path.expanduser(folder_path+"lightcurves_pickle/"+name+"_phot_ephem.dat")
		phases_ordered,Vmagphot,Vmagerrphot, ID_st, RA_st, DEC_st,l_st, b_st, Vmag_st, D_st, P_st, Amp_st,Red_st, ephemeris, MJD_found2, phase_test, phase_drake = cPickle.load(open(filen))

		#Vmagphot+Red_st if you want to fit for mags not extinction corrected

		
		f,p1,chi2=F.cutsFourier(phases_ordered,Vmagphot, LCcut=True,Gfill=True, ret = False) #you can give sigma = error on the magnitude (see Fourier.py) 
		
		vpuls, obs_mag_all, obs_phase_all = [], [], []

		plot_light  = True

		for ss in range(len(fd_new)):
			obs_MJD = jd[ss] - 2400000.5  
			obs_phase = F.phase(obs_MJD-ephemeris,P_st)
			obs_mag = f(obs_phase,p1) #magnitude at phase obs_phase
			obs_phase_all.append(obs_phase)
			obs_mag_all.append(obs_mag)

			if obs_phase > 0.947:
				obs_phase = 0.947
			T_alpha = templ_alpha(obs_phase)
			T_beta = templ_beta(obs_phase)
			Arv_alpha = 35.6* Amp_st +78.2
			Arv_beta = 42.1* Amp_st +51.1
			vpuls.append((Arv_alpha*T_alpha+Arv_beta*T_beta)/2.)
		
		if plot_light:		
			fig = plt.figure(figsize=(11.5, 7))
			plt.subplots_adjust(left=0.16, bottom=0.13, right=0.95, top=0.92, wspace=0.2)
			ax = fig.add_subplot(111)
			plt.errorbar(phases_ordered, Vmagphot,Vmagerrphot,capsize=0,fmt='.')
			plt.scatter(phases_ordered,Vmagphot, label = 'CSS data', color = 'blue') # Vmagphot is the magnitude , extinction corrected
			xx2 =np.linspace(-2,2.,1000)

			plt.plot(np.sort(phases_ordered),f(np.sort(phases_ordered),p1),  label = 'light curve fit',color = 'green', lw = 4)
			plt.plot(obs_phase_all,obs_mag_all, 'o',label = r'$\phi_{obs}$', color = 'red', ms = 20) # Vmagphot is the magnitude , extinction corrected

			textstr = 'P =' + str(round(P_st,3)) + ' days; $\eta$ = ' + str(ephemeris)
			ax.text(0.05, 0.90, textstr, transform=ax.transAxes, fontsize=22,
			                  verticalalignment='top')
		
			plt.title(target_names[j])
			ax.set_xlabel(r'$\phi$')
			ax.set_ylabel('$\mathrm{V (mag)}$')
			plt.xlim((0,1))

			#set size ticks for the main plot
			for item in ([ax.title,ax.yaxis.label,ax.xaxis.label]):
			    item.set_fontsize(32)

			for item in (ax.get_xticklabels() + ax.get_yticklabels()):
			    item.set_fontsize(26)

			plt.legend(loc = 4, fontsize = 22, frameon = False, numpoints= 1)
			
			fig.savefig(plot_path+target_names[j]+'.eps')


		#print 'vpuls',vpuls, np.mean(vpuls)
		vsys = vel_best - np.mean(vpuls)
		vgsr = vGSR(vsys, l_st, b_st)

		if obs_phase < 0.4:
			err_sys = 2
		if (obs_phase < 0.5) & (obs_phase > 0.4):
			err_sys = 0
		if (obs_phase > 0.5) & (obs_phase < 0.82):
			err_sys = -14+(14/0.5)*obs_phase
		if (obs_phase > 0.82):
			err_sys = -14+(14/0.5)*obs_phase
		err_tot = np.sqrt(err**2+sigma_vel_cal**2+err_sys**2)
		
		#save results
		myfile = open(folder_path+"results/observed_targets_"+night+'.txt','a') 
		myfile.write(night + '\t'+ name + '\t' + str(l_st) + '\t' + str(b_st) + '\t' +  str(RA_st) + '\t'  +str(DEC_st) + '\t' +\
		str(obs_phase) + '\t' +str(obs_mag) + '\t' + str(Amp_st) + '\t' + str(vel_best) + '\t' + str(minchi) + '\t' + str(err[0]) + '\t' + \
		str(Arv_alpha) + '\t' + str(T_alpha) + '\t' + str(Arv_beta) + '\t' + str(T_beta) + '\t' + str(vsys) + '\t' + \
		str(vgsr) + '\t' + str(err_tot[0]) + '\t' + str(Vmag_st) + '\t' + str(D_st) + '\t' + str(P_st) + '\t' + \
		str(Red_st) + '\t' + str(ephemeris) + '\t' +  str(RA_tel) + '\t' +  str(DEC_tel)+ '\t'  \
		   +  str(np.mean(vpuls)) + '\n') 
		myfile.close()
	
	return None

c =  299792.458

#these are the results for each night of #reload(resol)
#sigma_vel_cal, mean_sigma = resol.calibration_error()

sigmas = [12.723465477714441, 8.3742481968970246 , 9.7718080194882226, 14.172317339249838, 9.0600211435116105, 3.8823774213063063]
means = [-9.4653852149442552e-05,  -9.2634790369026325e-05, -8.8542384544618228e-05, -3.6228700237005092e-06, -3.8056416407812281e-06, -9.914964525885501e-05]

n = 2
#choose the night you want to reduce 'n1', 'n2', 'n3' , 'n4' , 'n5', 'n6'
sigma_vel_cal, mean_sigma = sigmas[n-1], means[n-1]
night = 'n'+str(n)

plot_target()


plt.clf()