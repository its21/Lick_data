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
    #computes the altitude in degrees of a given object at the given utc time
    #ra dec in degress  
    #yr mon day hr minu in utc at the time of the observation
    # lon lat are in degrees and lon is positive to the East    
    obs = ephem.Observer();
    obs.lon = lon # longi tude 
    obs.lat = lat #latitude
    obs.elevation=elev;
    fb = ephem.FixedBody();
    fb._ra=np.deg2rad(ra);
    fb._dec=np.deg2rad(dec)
    obs.date=ephem.Date('%d/%02d/%d %d:%d:0'%(yr,mon,day,hr,minu))
    fb.compute(obs);
    alt=np.rad2deg(1*(fb.alt))
    return alt

from astropy.time import Time
t = Time.now() #+ 10*u.hour #ok
cali_time = t - 7*u.hour #ok
RA, DEC, Vmag, D, P, Amp,Red, ephemeris = np.load('targets_15_20_info.npy')

ra = np.array(RA); dec = np.array(DEC)
alt = []
for i in range(len(ra)):
    al = getalt(ra[i],  dec[i], t.datetime.year, t.datetime.month, t.datetime.day, t.datetime.hour, t.datetime.minute)
    i, print(al)
    alt.append(al)
alt = np.array(alt)
airmass = 1./np.cos(np.deg2rad(90.-alt))
current_phase = np.round((phase(t.mjd-ephemeris, P)),2)
HAC_names = np.arange(0, len(ra),1)
wh = (current_phase>0.2) & (current_phase<0.6)
#print(HAC_names[wh], airmass[wh], current_phase[wh])
for i in range(len(ra[wh])):
    print(HAC_names[wh][i], round(ra[wh][i],4),  round(dec[wh][i],4),'alt=',alt[wh][i] ,'airmas =' , np.round(airmass[wh],2)[i], 'phase = ', current_phase[wh][i])

