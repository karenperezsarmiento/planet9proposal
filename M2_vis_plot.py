import make_el_chart as mec
import astropy.coordinates as apc
import numpy as np
from astropy import units as u
from datetime import datetime
import os
from astropy.time import Time

### These lines don't get used later - they're just to play around.
times = ['1999-01-01T00:00:00.123456789', '2010-01-01T00:00:00']
t = Time(times, format='isot', scale='utc')
t.sidereal_time('apparent', 'greenwich')

#########################################################################
### Define the coordinates and name of your object:
#obj_ra = apc.Angle('2h00m00s')
#obj_dec= apc.Angle('-3d00m00s')
### GEMS
#obj_ra = apc.Angle('3h30m00s')
#obj_dec= apc.Angle('-27d49m00s')
#skyobj = apc.SkyCoord(obj_ra, obj_dec, equinox = 'J2000')
#target='GEMS'
#obj_ra = apc.Angle('2h42m46s')
#obj_dec= apc.Angle('-21d32m00s')
#skyobj = apc.SkyCoord(obj_ra, obj_dec, equinox = 'J2000')
#target='PKS0240-217'
obj_ra = apc.Angle('9h35m13.3s')
obj_dec= apc.Angle('0d47m46s')
skyobj = apc.SkyCoord(obj_ra, obj_dec, equinox = 'J2000')
target='eFEDS'
#date_obs  = datetime.strptime('01-12-2020 23:00:00', '%d-%m-%Y %H:%M:%S')
date_obs  = datetime.strptime('15-4-2022 23:00:00', '%d-%m-%Y %H:%M:%S')

elMin=25.0
elStr = str(int(elMin))
mydir = os.getcwd()    # If you want to specify a directory, do so here.

mec.plot_visibility(date_obs,skyobj,elMin=elMin,
                    mylabel=target,filename = target+date_obs.strftime('%d%b%Y')+"_Visibility_above"+elStr,
                    mydir=mydir)

#mec.plot_contours(elMin=elMin)
