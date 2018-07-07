import numpy as np
import matplotlib.pyplot as plt
import Fourier as F
import cPickle
import numpy as np
import glob
import os
import scipy.signal as signal
from astropy.io import fits
import pyfits, pylab
import asciitable
#      CSS ID            RA (deg)  Dec (deg)  V      Period   Amp   Npts  Dist   Red    MJD peak     CSS ID num  Prior ID
#======================================================================================================================
lmin, lmax, bmin, bmax = 0, 360, -90, 90
extents = (lmin, lmax, bmin, bmax)

#READ CSS data
data0 = asciitable.read('/Users/iuliasimion/IDLWorkspace81/DATA/RR/apj452476t1_mrt.txt', quotechar="'", names = ['ID', 'RAdeg', 'DEdeg', '<Vmag>', 'Period', 'A', 'Num', 'dh', 'AV', 'Eta','AID'])

#READ CSS data (second paper)
data1 = asciitable.read('/Users/iuliasimion/IDLWorkspace81/DATA/RR/apj461344t2_mrt.txt', quotechar="'", names = ['ID', 'RAdeg', 'DEdeg', '<V>', 'P', 'A', 'Num', 'dh', 'AV','IDalt'])

#READ MSL data
data2 = asciitable.read('/Users/iuliasimion/IDLWorkspace81/DATA/RR/apj461344t1_mrt.txt', quotechar="'", names = ['ID', 'RAdeg', 'DEdeg', '<Vmag>', 'P', 'A', 'Num', 'dh', 'AV', 'Eta','IDalt'])

#READ SDSS RR Lyrae
data3 = asciitable.read('/Users/iuliasimion/IDLWorkspace81/DATA/RR/apj452476t2_mrt.txt', quotechar="'", names = ['ID', 'u0mag', 'g0mag', 'r0mag', 'i0mag', 'z0mag', 'VGSR', '[Fe/H]'])


#read RRL_params table1 in paper 1 Drake et al. 2013
temp=open('/Users/iuliasimion/IDLWorkspace81/DATA/RR/apj452476t2_mrt.txt')
lines=temp.readlines()
temp.close()

ID_sdss, ID_sdss_ok, vgsr, feh , vgsr_ok= [], [],[],[], []
for i in lines:
    temp=i.split() 
    print temp
    ID_sdss.append(temp[0])
    if len(temp) == 6:
        vgsr.append(-999)
        feh.append(0)
    if len(temp) == 7:
        vgsr.append(-999)
        feh.append(temp[6])
    if len(temp) == 8:
        vgsr_ok.append(temp[6])
        ID_sdss_ok.append(temp[0])
        vgsr.append(temp[6])
        feh.append(temp[7])

vgsr_ok = np.array(vgsr_ok[2:])
ID_sdss_ok = np.array(ID_sdss_ok[2:])

#read RRL_params table1 in paper 1 Drake et al. 2013
temp=open('/Users/iuliasimion/Observing2014/RRL/data/RRL_params')
lines=temp.readlines()
temp.close()
CSS_ID, RA, DEC, Vmag, Period, Amp, Npts, Dist, Red, MJD_peak, ID =[], [], [], [], [], [], [], [], [], [], []

for i in lines:
    temp=i.split() 
    RA.append(float(temp[1]))
    DEC.append(float(temp[2])) ; Vmag.append(float(temp[3]))
    Period.append(float(temp[4])); Amp.append(float(temp[5]))
    Npts.append(float(temp[6])) ; Dist.append(float(temp[7]))
    Red.append(float(temp[8])); MJD_peak.append(float(temp[9]))
    ID.append(float(temp[10]))
    CSS_ID.append(temp[0])

#CSS_J000343.1-134943   0.92972 -13.82868  17.66  0.586691  0.77  166   24.98  0.038  1012001005120  
#read paper 2 -- MISSING MJD PEAK
#    CSS ID           RA (deg)  Dec (deg)  V      Period   Amp  Npts  Dist   Red     CS ID num  Prior ID
#======================================================================================================================
temp2 = open('RRL_params_new')
lines = temp2.readlines()
temp2.close()
#CSS_ID, RA, DEC, Vmag, Period, Amp, Npts, Dist, Red, MJD_peak, ID = [], [], [], [], [], [], [], [], [], [], []

for i in lines:
    temp2=i.split() 
    RA.append(float(temp2[1]))
    DEC.append(float(temp2[2])) ; Vmag.append(float(temp2[3]))
    Period.append(float(temp2[4])); Amp.append(float(temp2[5]))
    Npts.append(float(temp2[6])) ; Dist.append(float(temp2[7]))
    Red.append(float(temp2[8])); MJD_peak.append(float(0))
    ID.append(float(temp2[9]))
    CSS_ID.append(temp[0])

#ra dec to l b
from bovy_coords import radec_to_lb
l_b = np.empty([len(RA),2])
l_b = radec_to_lb(RA,DEC,degree=True)
lon = l_b[:,0]
lat = l_b[:,1]
#lon_r=lon_r-360*(lon_r>180)
lon = np.array(lon,dtype=float)
lat = np.array(lat,dtype=float)

Dist = np.array(Dist) ; RA = np.array(RA) ; DEC = np.array(DEC); Vmag = np.array(Vmag)
Red = np.array(Red) ; MJD_peak = np.array(MJD_peak) ; ID = np.array(ID) ; Amp = np.array(Amp)
Period = np.array(Period); CSS_ID = np.array(CSS_ID)
Vmag_ext = Vmag+Red

#select targets
mask = np.where( (lon > 28) & (lon < 55) & (lat > -45) & (lat < -20) & (Vmag_ext < 17.5) & (Dist < 20) & (Dist > 15))
D_targets = Dist[mask]
RA_targets = RA[mask]
DEC_targets = DEC[mask]
Vmag_targets = Vmag[mask]
ID_targets =  ID[mask]
lon_targets = lon[mask]
lat_targets = lat[mask]
P_targets = Period[mask]
ext_targets = Red[mask]
Vmag_ext_targets = Vmag_ext[mask]
MJD_targets = MJD_peak[mask]  #MJD is zero for targets in the second file
phases_targets = F.phase(MJD_targets,P_targets)
Amp_targets = Amp[mask]
Red_targets = Red[mask]
CSS_ID_targets = CSS_ID[mask]


b = -22.046
a = -30.876
c = -2.627

periodISH = [] 
aamplitudeSH = []
pperiodSH = []
dSH = []
for i in range(len(Amp_targets)):
    det = b**2 - 4*a*(c-np.array(Amp_targets[i]))
    if (det >= 0):
        aamplitudeSH.append(Amp_targets[i])
        pperiodSH.append(P_targets[i])
        dSH.append(D_targets[i])
        logperiodB = (-b+math.sqrt(b**2-4*a*(c-Amp_targets[i])))/(2*a) 
        logperiodC = (-b-math.sqrt(b**2-4*a*(c-Amp_targets[i])))/(2*a)
        periodB = 10**(logperiodB)
        periodC = 10**(logperiodC)
        if (periodB > periodC):
            periodISH.append(periodB)
        else:
            periodISH.append(periodC)


deltaPSH = np.array(pperiodSH)-(np.array(periodISH))

aamplitudeSH = np.array(aamplitudeSH)
pperiodSH = np.array(pperiodSH)
dSH =np.array(dSH)

#choose the OoI population
x_bound = 0.05481
wh = np.where(deltaPSH < x_bound)

#OoI population
aamplitudeSH_bound = aamplitudeSH[wh]
pperiodSH_bound = pperiodSH[wh]
dSH_bound = dSH[wh]
dSH_O1 = dSH_bound
dSH_O2 = dSH[deltaPSH > x_bound]
'''
if os.path.exists('/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/name_Otype.txt'):
    os.remove('/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/name_Otype.txt')
for star in range(len(RA_targets)):
    myfile = open(os.path.join('/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced/', 'name_Otype.txt'),'a') 
    myfile.write('HAC'+str(star) +'\t'+str(deltaPSH[star]) +'\t' \
    +str(RA_targets[star]) +'\t' +str(DEC_targets[star]) + '\n') #doesnt shift the first exposure
    myfile.close()

#save params in a file
for i in range(len(D_targets)):
    myfile = open('RRL_targets_phases.txt', 'a') 
    myfile.write(str(i)+'\t'+str("{0:.5f}".format(RA_targets[i]))+'\t'+ str("{0:.5f}".format(DEC_targets[i]))+'\t' + \
            str(Vmag_ext_targets[i])+'\t' + str(ext_targets[i]) + '\t' + str(D_targets[i])+'\t'+ \
            str(MJD_targets[i])+'\t' + str(P_targets[i]) + '\t' + str(phases_targets[i])+ '\t' + \
            str(lon_targets[i])+'\t' + str(lat_targets[i])+ '\t' + `ID_targets[i]` +'\n')
    myfile.close() 


#save params in a file
for i in range(len(D_targets)):
    myfile = open('RRL_radec.txt', 'a') 
    myfile.write('['+str(i)+']'+ ' '+str("{0:.5f}".format(RA_targets[i]))+' '+ str("{0:.5f}".format(DEC_targets[i])) +'\n')
    myfile.close() 
'''
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx


# 0.568743859261

#read the IDs from all the photometry files
for star in range(200):
    print star
    RA_st = RA_targets[star]
    DEC_st = DEC_targets[star]
    l_st = lon_targets[star]
    b_st = lat_targets[star]
    Vmag_st = Vmag_targets[star]
    ID_st = ID_targets[star]
    D_st = D_targets[star]
    P_st = P_targets[star]
    phases_st = phases_targets[star]
    Amp_st = Amp_targets[star]
    Red_st = Red_targets[star]
    ephemeris = MJD_targets[star]
    print 'ephemeris=', ephemeris

    #read the light curve
    for i in range(1, 42):
        ph = np.loadtxt('RRinput'+ str(i)+'.phot', delimiter = ',')
        IDphot = ph[:,0]
        #mask = np.nonzero(IDphot == ID_st)
        #print mask
        mask = IDphot == ID_st
        if (np.sum(mask) > 0):
            mask = IDphot == ID_st
            MJDphot = (ph[:,1])[mask]
            Vmagphot = np.array((ph[:,2])[mask], dtype = float)
            Vmagerrphot = (ph[:,3])[mask]
            RAphot = (ph[:,4])[mask]
            DECphot = (ph[:,5])[mask]

            phases = F.phase(MJDphot,P_st) #P i the period you get from table1 #MJDphot is the ephemeris #phases is the phase #gives u the fractional part of the division
            phases_ok = F.phase((MJDphot-ephemeris),P_st) #P i the period you get from table1 #MJDphot is the ephemeris #phases is the phase #gives u the fractional part of the division
            phase_drake = F.phase(ephemeris,P_st)

            f,p1,chi2=F.cutsFourier(phases,Vmagphot, LCcut=True,Gfill=True, ret = False) #you can give sigma = error on the magnitude (see Fourier.py)
            f_eval = f(phases,p1)
            xx = np.linspace(0,1,100000,endpoint=True)
            idxmin = np.argmin(f(xx,p1))
            phase_ephem = xx[idxmin]
            eph_index = find_nearest(phases, xx[idxmin])
            MJD_nearest = MJDphot[eph_index]
            phase_nearest = phases[eph_index]

            #phases_new = phases-xx[idxmin]
            #phases_new[phases_new<0]+=1

            #phmax = phases[f_eval == min(f_eval)]
            #ephemeris_calc = MJDphot[f_eval == min(f_eval)]
            #V_ephem = Vmagphot[f_eval == min(f_eval)]

            n = (MJD_nearest - phase_nearest*P_st)/P_st
            nn = np.floor(n)
            MJD_found2 = xx[idxmin]*P_st + nn*P_st
            phase_test = F.phase(MJD_found2,P_st)
            
            print 'HERE:::::', phase_nearest, phase_drake, phase_test, xx[idxmin], MJD_nearest, MJD_found2, ephemeris, ephemeris_calc


            #Vmag_st are the magnitudes of the star NOT extinction corrected.
            filen = '/Users/iuliasimion/Observing2014/OBSERVING/light_curves/lightcurves_pickle/'+'HAC'+str(star)+'_phot_ephem.dat'
            if os.path.exists(filen):
                os.remove(filen)
        
            print 'saving HAC'+str(star)

            fil=open('/Users/iuliasimion/Observing2014/OBSERVING/light_curves/lightcurves_pickle/'+'HAC'+str(star)+'_phot_ephem.dat','a')
            cPickle.dump((phases_ok,Vmagphot,Vmagerrphot, ID_st, RA_st, DEC_st,l_st, b_st, Vmag_st, D_st, P_st, Amp_st,Red_st, ephemeris, MJD_found2, phase_test, phase_drake), fil)
            fil.close()
         
            fig = plt.figure(figsize=(13.5, 7))
            plt.subplots_adjust(left=0.11, bottom=0.13, right=0.95, top=0.94, wspace=0.2)
            ax = fig.add_subplot(111)
            plt.errorbar(phases_ok, Vmagphot,Vmagerrphot,capsize=0,fmt='.')
            plt.scatter(phases_ok,Vmagphot, label = 'ephemeris Drake') # Vmagphot is the magnitude , extinction corrected
            xx2 =np.linspace(-2,2.,1000)
            plt.plot(xx2-phase_test,f(xx2,p1), color = 'green', lw = 4, label = 'ephemeris fitted')
            plt.plot(xx2-phase_drake,f(xx2,p1), color = 'blue', lw = 4, label = 'ephemeris Drake')

            #textstr = 'ID =' + `ID_st`
            #ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=15,
                              #verticalalignment='top')
            textstr = 'P =' + str(round(P_st,3)) + ' days; $\eta$ = ' + str(ephemeris)
            ax.text(0.05, 0.90, textstr, transform=ax.transAxes, fontsize=20,
                              verticalalignment='top')
            textstr = '$\phi_{Drake}$ =' + str(round(phase_drake,4)) + '; $\phi_{fitted}$ = ' + str(round(phase_test,4))

            ax.text(0.05, 0.83, textstr, transform=ax.transAxes, fontsize=20, verticalalignment='top')
            #textstr = 'Vmag =' + str(Vmag_st) + '; D = ' + str(D_st) + ' kpc'
            #ax.text(0.65, 0.8, textstr, transform=ax.transAxes, fontsize=15,
            #                  verticalalignment='top')

            plt.title('HAC'+str(star))
            ax.set_xlabel('$\mathrm{phase}$')
            ax.set_ylabel('$\mathrm{V (mag)}$')
            plt.xlim((0,1))

            #set size ticks for the main plot
            for item in ([ax.title,ax.yaxis.label,ax.xaxis.label]):
                item.set_fontsize(28)

            for item in (ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(20)

            plt.legend(loc = 1)
            #plt.show()
            #plt.ylim((Vmag_st+0.4,min(f_eval)-0.2))
            fig.savefig('/Users/iuliasimion/Observing2014/hac_spectra/allyson/reduced_plots/HAC'+str(star)+'.eps')



