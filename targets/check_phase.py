import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import scipy.signal as signal
from astropy.io import fits
from bovy_coords import radec_to_lb
import ephem
from astropy.time import Time
import sys
import astropy.units as u

def phase(time,period):
    Phased = np.modf(time/period)[0]
    return Phased

def getalt(ra,dec, yr, mon, day, hr, minu, lon='-121:38:14', lat='37:20:35',  elev=1283):
    obs = ephem.Observer();
    obs.lon = lon # longitude 
    obs.lat = lat #latitude
    obs.elevation=elev;
    fb = ephem.FixedBody();
    fb._ra=np.deg2rad(ra);
    fb._dec=np.deg2rad(dec)
    obs.date=ephem.Date('%d/%02d/%d %d:%d:0'%(yr,mon,day,hr,minu))
    print(obs.date)
    fb.compute(obs);
    alt=np.rad2deg(1*(fb.alt))
    print(90-abs(alt))
    return alt

from astropy.time import Time
t = Time.now() 
cali_time = t - 7*u.hour #ok
RA, DEC, Vmag, D, P, Amp,Red, ephemeris = np.load('targets_15_20_info.npy')

ra = np.array(RA); dec = np.array(DEC)
alt = []
for i in range(len(ra)):
    al = getalt(ra[i],  dec[i], t.datetime.year, t.datetime.month, t.datetime.day, t.datetime.hour, t.datetime.minute)
    alt.append(al)
alt = np.array(alt)
airmass = 1./np.cos(np.deg2rad(90.-alt))
current_phase = np.round((phase(t.mjd-ephemeris, P)),2)
HAC_names = np.arange(0, len(ra),1)
wh = (airmass < 1.5) & (current_phase>0.2) & (current_phase<0.6)
print(HAC_names[wh], np.round(airmass[wh],2), current_phase[wh])
#print('HAC'+str(i), ': airmas =' , np.round(airmass,2), 'phase = ', current_phase)

