import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import scipy.signal as signal
from astropy.io import fits
import asciitable
from bovy_coords import radec_to_lb
import ephem
from astropy.time import Time
import sys
import astropy.units as u

def phase(time,period):
    Phased = np.modf(time/period)[0]
    return Phased

def getalt(ra,dec, yr, mon, day, hr, minu, lon='-121:38:14', lat='37:20:35',  elev=1283):
    #kitt peak: 31.9599° N, 111.5997° W, lon='-111:35:59', lat='31:57:12'

    #computes the altitude in degrees of a given object at the given utc time
    #ra dec in degress  
    #yr mon day hr minu in utc at the time of the observation
    # lon lat are in degrees and lon is positive to the East
    
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
t = Time.now() +10*u.hour#utc
cali_time = t - 7*u.hour #ok
RA, DEC, Vmag, D, P, Amp,Red, ephemeris = np.load('targets_15_20_info.npy')
ra = np.array(RA); dec = np.array(DEC)
alt = getalt(ra[0],  dec[0], t.datetime.year, t.datetime.month, t.datetime.day, t.datetime.hour, t.datetime.minute)
airmass = 1./np.cos(np.deg2rad(90.-alt))
print(np.round(airmass,2))

current_phase = np.round((phase(t.mjd-ephemeris, P)),2)
