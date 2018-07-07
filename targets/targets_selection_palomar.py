import numpy as np
import matplotlib.pyplot as plt
import Fourier as F
import cPickle
#      CSS ID            RA (deg)  Dec (deg)  V      Period   Amp   Npts  Dist   Red    MJD peak     CSS ID num  Prior ID
#======================================================================================================================
#read RRL_params table1 in paper 1 Drake et al. 2013
temp=open('RRL_params')
lines=temp.readlines()
temp.close()
CSS_ID, RA, DEC, Vmag, Period, Amp, Npts, Dist, Red, MJD_peak, ID = [], [], [], [], [], [], [], [], [], [], []

for i in lines:
	temp=i.split() 
	RA.append(float(temp[1]))
	DEC.append(float(temp[2])) ; Vmag.append(float(temp[3]))
	Period.append(float(temp[4])); Amp.append(float(temp[5]))
	Npts.append(float(temp[6])) ; Dist.append(float(temp[7]))
	Red.append(float(temp[8])); MJD_peak.append(float(temp[9]))
	ID.append(float(temp[10]))

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
Period = np.array(Period)
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

import finder_chart_ps1
reload(finder_chart_ps1)
from finder_chart_ps1 import deg2hour
ra_h, dec_h = deg2hour(np.array(RA_targets), np.array(DEC_targets), sep=":")

vgsr = [] ; observedby = []
for i in range(225):
  vgsr.append('   ')
  observedby.append('   ')
vgsr[211] = -195.4; vgsr[197] = -25.0
vgsr[37] = -258.5; vgsr[152] = -59.4
vgsr[22] = 195.2; vgsr[27] =84.0
vgsr[186] = 60.5; vgsr[91] = -162.4
vgsr[26] = -255.8; vgsr[131] = -28.9
vgsr[11] = -174.7; vgsr[177] = 205.7
vgsr[9] = 130.2; vgsr[195] = -244.3
vgsr[96] = 20; vgsr[105] = -118.5
vgsr[102] = 204.7; vgsr[145] = 108.3
vgsr[51] = 202.8; vgsr[45] = -64.3
vgsr[52] = -252.4; vgsr[42] = 170.8
vgsr[156] = 115.4; vgsr[94] = -166
vgsr[191] = -21.9; vgsr[194] = -270.9
vgsr[21] = 20.7; vgsr[208] = 133.3
vgsr[166] = 60.2; vgsr[154] = -304.8
vgsr[120] = 101.5; vgsr[121] = -53.0
vgsr[187] = 159.3; vgsr[169] = -1.8
vgsr[31] = -261.0; vgsr[48] = 200.6
vgsr[115] = -42.1; vgsr[60] = -96.4
vgsr[158] = -200.6; vgsr[113] = -60.8
vgsr[71] = -281.0; vgsr[47] = -25.0
vgsr[126] = -104.3; vgsr[169] = 30.2
vgsr[82] = 71.9; vgsr[143] = 55.2
vgsr[184] = 5.3; 

osbserved = np.array([211, 197, 37, 152, 22,27,186,91,26,131,11,\
    177,9,195,96,105,102,145,\
    51,45,52,42,156,94,191,194,21,208,166,154,\
    120,121,187,169,31,48,115,60,158,113,71,47,\
    126,169,82,143]) 

for i in range(len(osbserved)):
    observedby[osbserved[i]] = 'MDM'

observedby[184] = 'SDSS'
savefile = True
if savefile:
    #save params in a file
    for i in range(len(D_targets)):
        myfile = open('RRL_targets_palomar.txt', 'a') 
        myfile.write('HAC'+str(i)+'\t'+ str(ra_h[i])+'\t'+ str(dec_h[i])) #+'\t' + \
                #'V = '+ str(Vmag_ext_targets[i])+'\t' + 'D = '+ str((D_targets[i])) +'\t' + \
                #str(vgsr[i])+'\t' +observedby[i] +'\n')
        myfile.close() 

'''
targetsample = False
if targetsample:
    target_names, ra, dec, comments = [], [], [], []
    for i in range(len(D_targets)):
        target_names.append('HAC'+str(i))
        ra.append(str(ra_h[i]))
        dec.append(str(dec_h[i]))
        comments.append('V = '+ str(Vmag_ext_targets[i])+ ' ' + 'D = '+ str((D_targets[i])) +' ' + \
                ' '+ str(vgsr[i]) + ' ' +observedby[i])
'''

