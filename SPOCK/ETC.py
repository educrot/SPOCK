import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.io import ascii
from astropy.table import Table,Column, MaskedColumn
from scipy.interpolate import interp1d


class etc:
    def __init__(self,mag_val = None,mag_band = None,spt = None,filt = None,airmass = None,moonphase = None,irtf = None,\
                 num_tel = None,seeing = None,gain=None):

        self.mag_val = mag_val
        self.mag_band = mag_band
        self.spt = spt
        self.filt = filt
        self.airmass = airmass
        self.moonphase = moonphase
        self.irtf = irtf
        self.num_tel =num_tel
        self.seeing = seeing
        self.gain = gain
        self.airmass = airmass

        """
        Class to Make schedules for the target list, observatory, dtae_range and startegy indicated

        Parameters
        ----------
        target_list: list of all target with their basic info
        observatory: list of the observatory , and observatory/telescope if several telescope per observatory
        startegy: scheduling strategy, either 'continuous' or 'segmented'
        duration_segement: if strategy segmented chosen need to specify the duration of each segment
        nb_segments: if strategy segmented chosen, number of target per night

        Returns
        -------
        night blocks for the night and telescope according to the strategy selected

        """

        ########## CCD
        self.npix        = 2048 #Number of pixels x
        self.npiy        = 2048 #Number of pixels y
        self.mipi        = 13.5 #picelscale
        #self.gain        = 1.1  #Gain [el/ADU] # SSO
        self.tqe         = 243
        self.qefac       = 1.0  #coefficient
        self.d0          = 2e4  #dark_293K
        self.ron         = 10
        self.tccd        = -60  #CCD temperature [degree C]
        self.binning     = 1    #CCD binning
        self.tlost       = 10   #read-out & overhead time [s]

        ########## Optics
        self.m1dia       = 1030 #M1 free aperture [mm]
        self.m2dia       = 300  #M2 mechanical aperture [mm]
        self.focrat      = 8.0  #Focal ratio
        self.mloss       = 10   #system loss [%]

        ########## Observatory
        self.alti        = 2500 #Altitude [m]
        self.num_tel     = 1    #Number of telescopes

        ########## Light curve
        self.bin_lc      = 7.2  #Time bin for SNR [min]
        self.rednoise    = 0.   #Red noise[ppm]
        self.nsigma      = 5.   #

        ########## constants

        self.rsun        = 695800.
        self.rearth      = 6371.
        self.c           = 299792458.
        self.scinfac     = 0.09
        self.h           = 6.62607E-34
        self.e           = 2.71828

        self.c1_file="/Users/elsaducrot/spock/SPOCK/files_ETC/coating_1.dat"
        self.c1=ascii.read(self.c1_file, data_start=0)
        #plt.grid(True)
        #plt.xlabel("Wavelength [nm]")
        #plt.ylabel("Efficiency")
        #plt.plot(c1['col1'],c1['col2'])
        #plt.show()

        self.c2_file = "/Users/elsaducrot/spock/SPOCK/files_ETC/coating_2.dat"
        self.c2 = ascii.read(self.c2_file, data_start=0)
        # plt.grid(True)
        # plt.xlabel("Wavelength [nm]")
        # plt.ylabel("Efficiency")
        # plt.plot(c2['col1'],c2['col2'])
        # plt.show()

        self.ccd_file="/Users/elsaducrot/spock/SPOCK/files_ETC/ccd.dat"
        self.ccd=ascii.read(self.ccd_file, data_start=5)
        #plt.grid(True)
        #plt.xlabel("Wavelength [nm]")
        #plt.ylabel("Efficiency")
        #plt.plot(ccd['col1'],ccd['col2'])
        #plt.show()

        self.qet_file="/Users/elsaducrot/spock/SPOCK/files_ETC/qet.dat"
        self.qet=ascii.read(self.qet_file, data_start=0)
        #plt.grid(True)
        #plt.xlabel("Wavelength [nm]")
        #plt.ylabel("Efficiency")
        #plt.plot(qet['col1'],qet['col2'])
        #plt.show()

        self.window_file="/Users/elsaducrot/spock/SPOCK/files_ETC/window.dat"
        self.window=ascii.read(self.window_file, data_start=0)
        #plt.grid(True)
        #plt.xlabel("Wavelength [nm]")
        #plt.ylabel("Efficiency")
        #plt.plot(window['col1'],window['col2'])
        #plt.show()

        #System
        self.mloss=1-0.01*self.mloss
        self.surf=(np.pi*(self.m1dia*0.5*1.0e-3)**2)-(np.pi*(self.m2dia*0.5*1.0e-3)**2)
        self.eape=2*(self.surf/np.pi)**.5    #effective aperture

        self.c1['col2']=self.c1['col2']*self.mloss
        self.c2['col2'] = self.c2['col2']*self.mloss

        #CCD corrections
        if self.binning > 1:
            self.mipi=self.mipi*self.binning
            self.npix=self.npix/self.binning
            self.npiy=self.npiy/self.binning
            self.dark=self.dark*self.binning**2
        self.pixelscale = self.mipi * 206.265/(self.focrat*self.m1dia)
        self.fov = self.npix*self.pixelscale/60.

        #convert to percent
        self.rednoise=self.rednoise/1.0e6
        #convert from celsius to K
        self.tccd=self.tccd+273
        dt=(self.tccd-self.tqe)/125.

        #Dark current
        self.dark=(self.d0*122*(self.tccd-10)**3)*self.e**(-6400/(self.tccd-10))

        #Correction of ccd response for low temperatures
        self.qet['col2']=(1.-self.qet['col2'])*self.qefac*np.abs(dt)
        self.qet['col2']=(1.-self.qet['col2'])
        if dt < 0.:
            self.ccd['col2']=self.ccd['col2']*self.qet['col2']

        bg_file="/Users/elsaducrot/spock/SPOCK/files_ETC/background.dat"
        bg=ascii.read(bg_file, data_start=0)

        #plt.xlim(500,510)

        #plt.plot(bg['col1'],bg['col6'],bg['col1'],bg['col5'],bg['col1'],bg['col4'],bg['col1'],bg['col3'],bg['col1'],bg['col2'])
        moonph_sin=np.arcsin(self.moonphase)*360/(np.pi)
        #Moonphase dependence of background
        if moonph_sin>=0. and moonph_sin <45.:
            m1=0
            m2=1
            moonzero=0.
        elif  moonph_sin>=45. and moonph_sin <90.:
            m1=1
            m2=2
            moonzero=45.
        elif  moonph_sin>=90. and moonph_sin <135.:
            m1=2
            m2=3
            moonzero=90.
        elif  moonph_sin>=135. and moonph_sin <=180.:
            m1=3
            m2=4
            moonzero=135.
        else:
            print("moonphase unrealistic")
            print(5/0.)

        #Airmass dependence of background
        if self.airmass>=1. and self.airmass <1.5:
            l1=2
            l2=3
            airzero=1.
        elif self.airmass>=1.5 and self.airmass <2.:
            l1=3
            l2=4
            airzero=1.5
        elif self.airmass>=2. and self.airmass <2.5:
            l1=4
            l2=5
            airzero=2.
        elif self.airmass>=2.5 and self.airmass <3.:
            l1=5
            l2=6
            airzero=2.5
        else:
            print("airmass not supported (1-3)")
            print(5/0.)

        #only use first 5 entries of background.dat!!!
        exdif=np.array(bg['col'+str(int(l1))][:5]-bg['col'+str(int(l2))][:5])
        backe=bg['col'+str(int(l1))][:5]-exdif*((self.airmass-airzero)/0.5)
        exdif=backe[m1]-backe[m2]

        #
        self.back = Table(self.ccd)
        self.back['col2']=np.zeros(len(self.back['col2']))+backe[m1]-exdif*((moonph_sin-moonzero)/45.)
        #plt.grid(True)
        #plt.xlabel("Wavelength [nm]")
        #plt.ylabel("Background")
        #plt.plot(back['col1'],back['col2'])
        #plt.show()

        extind_file="/Users/elsaducrot/spock/SPOCK/files_ETC/extin.dat"
        self.extind=ascii.read(extind_file, data_start=0)

        #plt.plot(extind['col1'],extind['col6'],extind['col1'],extind['col5'],extind['col1'],extind['col4'], \
        #         extind['col1'],extind['col3'],extind['col1'],extind['col2'])#,extind['col1'],extin,'o')

        #select right extiction
        self.exdif=self.extind['col'+str(int(l1))]-self.extind['col'+str(int(l2))]
        self.extin=self.extind['col'+str(int(l1))]-self.exdif*((airmass-airzero)/0.5)

        #plt.grid(True)
        #plt.xlabel("Wavelength [nm]")
        #plt.ylabel("Extiction")
        #plt.plot(extind['col1'],extin)
        #plt.show()

        #available spectral types
        self.spectra=Table(names=['type','vmj','vref','rs','file'],
                      data=[['B0','B1','B3','B6','B8', \
                        'A0','A2','A3','A5', \
                        'F0','F2','F5','F8', \
                        'G0','G1','G2','G5','G8', \
                        'K0','K2','K5','K7', \
                        'M0','M2','M4','M5','M6','M7','M8','M9', \
                        'L2','L5','L8'], \
                        [-0.8,-0.73,-0.60,-0.46,-0.36, \
                        -0.16,-0.07,-0.02,0.09, \
                        0.37,0.48,0.67,0.79, \
                        1.03,1.055,1.08,1.25,1.32, \
                        1.46,1.81,2.22,2.71, \
                        2.96,3.58,4.42,5.39, 0.00, 0.00,0.00,0.00, \
                        0.00,0.00,0.00],\
                        [0.,0.,0.,0.,0., \
                        0.,0.,0.,0., \
                        0.,0.,0.,0., \
                        0.,0.,0.,0.,0., \
                        0.,0.,0.,0., \
                        0.,0.,0.,0.,7.09,9.78,9.91,9.54, \
                        13.41,12.83,13.25],
                        [7.4,6.5,4.8,3.7,3.0, \
                         2.4,2.15,2.0,1.7, \
                         1.5,1.4,1.3,1.2, \
                         1.1,1.05,1.,0.92,0.88, \
                         0.85,0.80,0.72,0.67, \
                         0.60,0.44,0.26,0.18,0.135,0.12,0.105,0.09, \
                         0.105,0.105,0.105],
                         ['b0_pickles.dat','b1_pickles.dat','b3_pickles.dat','b6_pickles.dat','b8_pickles.dat', \
                          'a0_pickles.dat','a2_pickles.dat','a3_pickles.dat','a5_pickles.dat', \
                          'f0_pickles.dat','f2_pickles.dat','f5_pickles.dat','f8_pickles.dat', \
                          'g0_ltt7379.dat','g1_pickles.dat','g2_pickles.dat','g5_pickles.dat','g8_pickles.dat', \
                          'k0_pickles.dat','k2_pickles.dat','k5_pickles.dat','k7_pickles.dat', \
                          'm0_pickles.dat','m2_pickles.dat','m4_pickles.dat','m5_pickles.dat', \
                          'm6_gl406.dat','m7_gj644c.dat','m8_vb10.dat','m9_den1048.dat', \
                          'l2_kelu1.dat','l5_2mass1507.dat','l8_den0255.dat']])
        #Changed vref from 10.23 0. for G0 standard star
        spt_sel=np.array(self.spectra['type'].data)

        #get spectral type information
        try:
            self.i = np.where(self.spt==spt_sel)[0][0]
        except:
            print("spectral type not in list")
            #print(5/0)

        #available spectra are in folder Spectra
        path='/Users/elsaducrot/spock/SPOCK/files_ETC/Spectra/'
        spec_file = os.path.join(path,self.spectra['file'][self.i])
        self.spec=ascii.read(spec_file, data_start=0)
        #plt.grid(True)
        #plt.xlabel("Wavelength [nm]")
        #plt.ylabel("Intensity")
        #plt.title(spectra['file'][i][:-4])
        #plt.plot(spec['col1'],spec['col2'])
        #plt.show()

        # available filters are in folder Filters, check available files
        path = '/Users/elsaducrot/spock/SPOCK/files_ETC/Filters/'
        files = []
        # r=root, d=directories, f = files
        for r, d, f in os.walk(path):
            for file in f:
                files.append(file)
        files = np.sort(files)
        # convert to lower case
        files_low = np.array([x.lower() if isinstance(x, str) else x for x in files])
        # check, if selceted filter is contained
        # in_filter_list=np.where(np.char.find(files_low, filt.lower()+'.dat')!=-1)
        try:
            in_filter_list = list(files_low).index(filt.lower() + '.dat')
            filter_file = os.path.join(path, files[in_filter_list])
            self.filter = ascii.read(filter_file, data_start=0)
            # plt.grid(True)
            # plt.xlabel("Wavelength [nm]")
            # plt.ylabel("Troughput")
            # plt.title(filt)
            # plt.plot(filter['col1'], filter['col2'])
            # plt.show()
        except ValueError:
            sys.exit("Filter not available")

        #get spectral type information
        try:
            i=np.where(self.spt==spt_sel)[0][0]
        except:
            print("spectral type not in list")
            print(5/0)

        #frequency instead of wavelength
        self.ener=self.h*self.c/(self.spec['col1']*1e-9)

        ####colour correction
        if self.mag_band == "V":
            #do nothing
            self.mag_val=self.mag_val+0.

        elif self.mag_band == "J":
            #apply V-J correction for spectral type
            self.mag_val=self.mag_val+self.spectra['vmj'][self.i]
        else:
            #ToDo add K-magnitude
            print("Band not implemented")
            print(5/0.)

        #Apply IRTF/UCS flux correction to photomeric standards of M-dwarfs
        #Standards were obtained with
        if len(np.where(self.spt==spt_sel[np.where(spt_sel=='M6')[0][0]:])[0]) > 0.0:
            corcal=self.irtf
        else:
            corcal=1.
        self.spec['col2']=self.spec['col2']*corcal

        #Correction factor, accounting the apparent magnitude of the target.
        self.corflux = 10**((self.spectra['vref'][i]-self.mag_val)/2.5)


    def peak_calculation(self,exp_t = 10):

        #Object independend effective troughput of the system
        self.effi=self.window['col2']*self.ccd['col2']*self.c1['col2']*self.c2['col2']*self.filter['col2']

        #Background collected by system
        back2=Table(self.back)
        back2['col2'] = np.array(self.back['col2'])*self.effi*exp_t*(5/1000.)*self.surf*self.pixelscale**2
        tback=sum(back2['col2'])    #Background [el/pixel]

        #print(tback)
        #Starlight collected by system
        spec2=Table(self.spec)
        spec2['col2']=(self.spec['col2']/self.ener)*self.corflux*self.extin*self.effi*exp_t*(5/1000.)*self.surf

        signal=sum(spec2['col2'])
        #print(signal)

        #light curve calculations
        bin_lc_c=self.bin_lc*60
        npbin=bin_lc_c/(exp_t+self.tlost)  #Nexp per bin
        earthtra=(self.rearth/self.spectra['rs'][self.i])**2    #earth transit depth for spectral type

        #sum (I do not really understand this)
        #lambsum=sum(spec2['col2'])
        #lambeff=sum(spec2['col1']*spec2['col2'])/lambsum

        seeing_p=self.seeing/self.pixelscale
        nape=np.pi*(2*seeing_p)**2    #Npixels in aperture

        tdark=self.dark*exp_t      #Dark [el\pixel]
        tdarkape=tdark*nape
        tbackape=tback*nape
        tronape=nape*self.ron**2
        peak=signal*0.66/(seeing_p)**2

        #scintilations (from model)
        scinti=self.scinfac*((self.eape*100.)**(-0.6666))*self.airmass**(1.75)
        scinti_alt=scinti*self.e**(-self.alti/8000.)
        scinti_exp=scinti_alt/np.sqrt(2*exp_t)
        scinti2=(scinti_exp*signal)**2
        snr = signal/np.sqrt(signal+tbackape+tronape+tdarkape+scinti2)

        #more (similar) telescopes used?
        if self.num_tel > 1:
            snr=snr*np.sqrt(self.num_tel)

        snrbin = snr*np.sqrt(npbin)       #snr for each bin in the light-curve
        errorbin = 1/snrbin
        errorbin_rn = np.sqrt(errorbin**2 + self.rednoise**2)
        #planet sensitivity
        sensi=np.sqrt(self.nsigma*errorbin_rn/earthtra)

        print("Peak [ADU]:\t", peak/self.gain)
        print("Sky [ADU]:\t", tbackape/self.gain)
        peak_ADU = peak/self.gain
        sky_gain = tbackape/self.gain
        return

    def exp_time_calculator(self,ADUpeak = None):
        ##########  Observation
        # ADUpeak = 33000  # counts in peak desired [ADU]
        # seeing = 0.95  # effective seeing [arcsec]
        # airmass = 1.2  # Airmass
        # num_tel = 1

        peak=ADUpeak*self.gain

        #Object independend effective troughput of the system
        self.effi = self.window['col2']*self.ccd['col2']*self.c1['col2']*self.c2['col2']*self.filter['col2']

        #Starlight collected by system
        spec2=Table(self.spec)
        spec2['col2']=(self.spec['col2']/self.ener)*self.corflux*self.extin*self.effi*(5/1000.)*self.surf
        signal=sum(spec2['col2'])

        #calculate exposure time from desired peak value
        seeing_p=self.seeing/self.pixelscale
        #peak=self.exp_t*signal*0.66/(seeing_p)**2

        exp_t=peak*(seeing_p)**2/(signal*0.66)
        # print("Exposure [s]:\t",exp_t)

        #Background collected by system
        back2=Table(self.back)
        back2['col2'] = np.array(self.back['col2'])*self.effi*(5/1000.)*self.surf*self.pixelscale**2
        tback=sum(back2['col2'])*exp_t    #Background [el/pixel]

        #light curve calculations
        bin_lc_c=self.bin_lc*60
        npbin=bin_lc_c/(exp_t+self.tlost)  #Nexp per bin

        nape=np.pi*(2*seeing_p)**2    #Npixels in aperture
        earthtra=(self.rearth/self.spectra['rs'][self.i])**2    #earth transit depth for spectral type

        #sum (I do not really understand this)

        lambsum=sum(spec2['col2'])
        lambeff=sum(spec2['col1']*spec2['col2'])/lambsum

        tdark=self.dark*exp_t      #Dark [el\pixel]
        tdarkape=tdark*nape
        tbackape=tback*nape
        tronape=nape*self.ron**2

        #scintilations
        scinti=self.scinfac*((self.eape*100.)**(-0.6666))*self.airmass**(1.75)
        scinti_alt=scinti*self.e**(-self.alti/8000.)
        scinti_exp=scinti_alt/np.sqrt(2*exp_t)
        scinti2=(scinti_exp*signal)**2
        snr = signal/np.sqrt(signal+tbackape+tronape+tdarkape+scinti2)

        #more (similar) telescopes used?
        if self.num_tel > 1:
            snr=snr*np.sqrt(self.num_tel)

        snrbin = snr*np.sqrt(npbin)       #snr for each bin in the light-curve
        errorbin = 1/snrbin
        errorbin_rn = np.sqrt(errorbin**2 + self.rednoise**2)
        #planet sensitivity
        sensi=np.sqrt(self.nsigma*errorbin_rn/earthtra)

        # print("Peak [ADU]:\t", peak/self.gain)
        # print("Sky [ADU]:\t", tbackape/self.gain)
        peak_ADU = peak/self.gain
        sky_ADU =  tbackape/self.gain
        return exp_t, peak_ADU,sky_ADU
