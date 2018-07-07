import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import scipy.signal as signal
from astropy.io import fits
import asciitable
from bovy_coords import radec_to_lb

def phase(time,period):
    Phased = np.modf(time/period)[0]
    return Phased

temp=open('RRL_params_2018.txt')
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

ra = np.array(RA)
dec = np.array(DEC)
l_b = radec_to_lb(ra,dec,degree=True)
lon = l_b[:,0]
lat = l_b[:,1]
lon=lon-360*(lon>180)
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
Vmag_targets = Vmag_ext[mask]
Vmagpaper_targets = Vmag[mask]
ephemeris = MJD_peak[mask] 
#phases_targets = phase(MJD_targets,P_targets)
Amp_targets = Amp[mask]
Red_targets = Red[mask]
CSS_ID_targets = CSS_ID[mask]


#Vmag_st are the magnitudes of the star NOT extinction corrected.
filen1 = '/Users/iuliasimion/work/2018/observing/Lick/Lick_data/targets/targets_15_20_info'
if os.path.exists(filen1+'npy'):
    os.remove(filen1+'npy')
np.save(filen1, (RA_targets, DEC_targets, Vmag_targets, D_targets, P_targets, Amp_targets,Red_targets, ephemeris))

#save params in a file
for i in range(len(RA_targets)):
    myfile = open('obs_targets_info.txt', 'a') 
    myfile.write('HAC'+str(i)+'\t'+str("{0:.5f}".format(RA_targets[i]))+'\t'+ str("{0:.5f}".format(DEC_targets[i]))+'\t' + \
            str(Vmagpaper_targets[i])+'\t' + str(ext_targets[i]) + '\t' + str(D_targets[i])+'\t'+ \
            str(ephemeris[i])+'\t' + str(P_targets[i]) + '\t' + \
            str(lon_targets[i])+'\t' + str(lat_targets[i]) +'\n')
    myfile.close() 


#save params in a file
for i in range(len(D_targets)):
    myfile = open('obs_radec.txt', 'a') 
    myfile.write('HAC'+str(i)+ ' '+str("{0:.5f}".format(RA_targets[i]))+' '+ str("{0:.5f}".format(DEC_targets[i])) +'\n')
    myfile.close() 
