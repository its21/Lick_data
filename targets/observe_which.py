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
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)


##############################################################################
# Import the packages necessary for finding coordinates and making
# coordinate transformations

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

def phase(time,period):
    Phased = np.modf(time/period)[0]
    return Phased

RA, DEC, Vmag, D, P, Amp,Red, ephemeris = np.load('targets_15_20_info.npy')
ra_hac = np.array(RA); dec_hac = np.array(DEC)
from astropy import units as u
c = SkyCoord(ra=ra_hac*u.degree, dec=dec_hac*u.degree, frame='icrs')
#can get other coordinates e.g. c.ra.hour  or c.to_string('dms')
#'10d41m04.488s 41d16m09.012s'
#c.to_string('hmsdms')
#or to galactic: c_icrs.galactic  

##############################################################################
# Use `astropy.coordinates.EarthLocation` to provide the location of telescope
# and set the time to UTC
#lon='-121:38:14', lat='37:20:35', elev=1283
#Shane_3m = EarthLocation(lat=37.2*u.deg, lon=-121*u.deg, height=1283*u.m)
lick = EarthLocation.of_site('lick')
#time = Time('2012-7-12 23:00:00') #in utc
utcoffset = +3*u.hour
time = Time.now() +utcoffset #utc

##############################################################################
# `astropy.coordinates.EarthLocation.get_site_names` and
# `~astropy.coordinates.EarthLocation.get_site_names` can be used to get
# locations of major observatories.
#
# Use `astropy.coordinates` to find the Alt, Az coordinates of M33 at as
# observed from Bear Mountain at 11pm on 2012 July 12.
hac_altaz = c.transform_to(AltAz(obstime=time,location=lick))
#hac_altaz = m33.transform_to(AltAz(obstime=time,location=lick))
#print("HAC Altitude = {0.alt:.2}".format(hac_altaz))
#print('HAC airmass = ',1./np.cos(np.deg2rad(90.-hac_altaz.alt.value)))
airmass = hac_altaz.secz

current_phase = np.round((phase(time.mjd-ephemeris, P)),2)
HAC_names = np.arange(0, len(ra_hac),1)
wh = (current_phase>0.2) & (current_phase<0.6)  & (airmass<2.5) & (airmass>0)
status = np.zeros(len(ra_hac))

observed_mdm = np.array([197,  37, 152, 22, 27, 186, 91, 26, 131, 11,\
	177, 9, 195, 96, 105, 102, 145, 51, 45, 52, 42, 156, 94, 191,\
	194, 21, 166, 154, 120, 121, 187, 169, 31, 48, 115, 60, 158, 113, 71, 47, 126, 169, 82, 143, 184])

observed_kast = np.array([41, 16, 72, 88, 98, 158, 189])
status[observed_mdm] = 1
status[observed_kast] = 2
#print(HAC_names[wh], airmass[wh], current_phase[wh])
for i in range(len(ra_hac[wh])):
    print('HAC '+str(HAC_names[wh][i]), status[wh][i],round(ra_hac[wh][i],4),  round(dec_hac[wh][i],4), 'Vmag=',Vmag[wh][i], 'dist=',D[wh][i],'airmass =' , np.round(airmass[wh],2)[i], 'phase = ', current_phase[wh][i])
