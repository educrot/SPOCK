import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.io import ascii
from astropy.table import Table,Column, MaskedColumn
from scipy.interpolate import interp1d


import ETC

##########  Observation
mag_val     = 11.92 #star magnitude
mag_band    = 'J'   #Band (J or V)
spt         = 'M6'  #spectral type
filt        = 'r' #Filter in use
exp_t       = 100    #Exposure time
ADUpeak     = 33000 #counts in peak desired [ADU]
seeing      = 0.95  #effective seeing [arcsec]
airmass     = 1.2   #Airmass
moonphase   = 0.6   #Moon phase (0.=new moon, 1.=full moon)
irtf        = 0.8   #IRTF/UCS flux correction


a = ETC.etc(mag_val = 11.9,mag_band = 'J',spt = 'M6',filt = 'r',airmass = 1.2,moonphase = 0.6,irtf = 0.8,num_tel = 1,seeing = 0.95)

texp = a.exp_time_calculator(ADUpeak = 33000)


print()