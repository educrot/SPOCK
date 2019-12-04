#!/usr/bin/python
import os
import shutil
from astropy.table import Table
from astroplan import Observer
from astropy.time import Time
from SPOCK.test_txtfiles import startup, startup_no_flats, Path_txt_files, flatexo_gany, flatexo_io, flatexo_euro, first_target_offset, flatexo_artemis_morning, flatexo_artemis_evening, startup_artemis,flatexo_saintex
from SPOCK.test_txtfiles import first_target,target, flatdawn, biasdark, shutdown, flatexo_calli, flatdawn_no_flats, target_no_DONUTS, target_offset, biasdark_comete, flatdawn_artemis
from astropy.coordinates import SkyCoord, get_sun, AltAz, EarthLocation
from astropy import units as u

#initialisation
filt={}
texp={}
index={}
ra1={}
dec1={}
ra2={}
dec2={}
ra3={}
dec3={}


def make_np(t_now,nb_jours,tel):

	telescope=tel
	print(t_now,nb_jours,telescope)
	t0=Time(t_now)
	dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

	for nb_day in range(0,nb_jours):
		print(telescope)
		t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso
		Path='./DATABASE'
		p=os.path.join(Path,str(telescope),'Plans_by_date',str(t_now))
		print(p)
		if not os.path.exists(p):
			os.makedirs(p)
		# if not os.path.exists(p2):
		#     os.makedirs(p2)

		#Read night_blocks file
		scheduler_table=Table.read('./DATABASE/' + str(telescope) +'/night_blocks_'+ str(telescope) +'_' + str(t_now)+'.txt', format='ascii')
		name=scheduler_table['target']
		date_start=scheduler_table['start time (UTC)']
		date_end=scheduler_table['end time (UTC)']
		ra1=scheduler_table['ra (h)']
		ra2=scheduler_table['ra (m)']
		ra3=scheduler_table['ra (s)']
		dec1=scheduler_table['dec (d)']
		dec2=scheduler_table['dec (m)']
		dec3=scheduler_table['dec (s)']
		config_mask=scheduler_table['configuration']

		scheduler_table.add_index('target')
		try:
			index_to_delete=scheduler_table.loc['TransitionBlock'].index
			scheduler_table.remove_row(index_to_delete)
		except KeyError:
			print(scheduler_table)

		config_filled=config_mask#.filled()


		for i in range(0,len(scheduler_table)):
			if name[i]!='TransitionBlock':
				config=config_filled[i].split(',')
				config_tup=tuple(config)
				for item in config_tup:

					# if item.find('index=',1)!=-1:
					#     item=item.replace('{','')
					#     item=item.replace('}','')
					#     item=item.replace('\'','')
					#     index[i]=item.replace('index=','')
					if item.find('filt=',1)!=-1:
						item=item.replace('{','')
						item=item.replace('}','')
						item=item.replace('\'','')
						item=item.replace('\"','')
						filt[i]=item.replace('filt=','')
						filt[i]=filt[i].replace(' ','')

					if item.find('texp=',1)!=-1:
						item=item.replace('{','')
						item=item.replace('}','')
						item=item.replace('\'','')
						texp[i]=item.replace('texp=','')
						texp[i]=texp[i].replace(' ','')
				print('avant :',filt[i])
				if filt[i]=='z' or filt[i]=='g' or filt[i]=='g' or filt[i]=='i' or filt[i]=='r':
					print('ici')
					a=filt[i]
					filt[i]=a+'\''
				print('apres :',filt[i],texp[i])

		location = EarthLocation.from_geodetic(-70.40300000000002*u.deg, -24.625199999999996*u.deg,2635.0000000009704*u.m)
		paranal = Observer(location=location, name="paranal", timezone="UTC")
		t=Time(t_now)
		sun_set_jd=paranal.sun_set_time(t,which='next')
		sun_rise_jd=paranal.sun_rise_time(t+1,which='next')
		sun_set=Time(sun_set_jd,format='jd',scale='utc')
		sun_rise=Time(sun_rise_jd+1,format='jd',scale='utc')

		#Defaults parameters
		autofocus=None
		waitlimit=600
		afinterval=60
		count='5000'

		#Set & rise time between
		sun_set=paranal.sun_set_time(t,which='next')
		twilight_evening_nautic_jd=paranal.twilight_evening_nautical(t,which='nearest')#,horizon=-6*u.deg)
		twilight_evening_civil_jd=paranal.twilight_evening_civil(t,which='nearest')#,horizon=0*u.deg)
		twilight_evening_astro_jd=paranal.twilight_morning_astronomical(t,which='nearest')#,horizon=0*u.deg)

		sun_rise=paranal.sun_rise_time(t,which='next')
		twilight_morning_nautic_jd=paranal.twilight_morning_nautical(t,which='nearest')#,horizon=-6*u.deg)
		twilight_morning_civil_jd=paranal.twilight_morning_civil(t,which='nearest')#,horizon=0*u.deg)
		twilight_morning_astro_jd=paranal.twilight_morning_astronomical(t,which='nearest')#,horizon=0*u.deg)

		twilight_evening_between_civil_nautic=Time((twilight_evening_nautic_jd.jd+twilight_evening_civil_jd.jd)/2,format='jd',scale='utc')
		twilight_morning_between_civil_nautic=Time((twilight_morning_nautic_jd.jd+twilight_morning_civil_jd.jd)/2,format='jd',scale='utc')
		# if telescope.find('Artemis'):
		#location_SNO = EarthLocation.from_geodetic(28.2999988*u.deg, -16.50583131*u.deg, 2390*u.m)
		#teide = Observer(location=location_SNO, name="SNO", timezone="UTC")
		location_SNO = EarthLocation.from_geodetic(-16.50583131*u.deg, 28.2999988*u.deg, 2390*u.m)
		teide = Observer(location=location_SNO, name="SNO", timezone="UTC")

	#	teide = Observer.at_site('Roque de los Muchachos')
		sun_set_teide=teide.sun_set_time(t,which='next')
		sun_rise_teide=teide.sun_rise_time(t+1,which='next')
		location_saintex = EarthLocation.from_geodetic(-115.48694444444445*u.deg, 31.029166666666665*u.deg, 2829.9999999997976*u.m)
		san_pedro = Observer(location=location_saintex, name="saintex", timezone="UTC")
		#san_pedro = Observer.at_site("Observatorio Astronomico Nacional, San Pedro Martir")
		sun_set_san_pedro=san_pedro.sun_set_time(t+1,which='next')
		sun_rise_san_pedro=san_pedro.sun_rise_time(t+1,which='next')

		#Making plans
		Path=Path_txt_files(telescope)
		print(Path)
		if telescope.find('Europa') is not -1:
			startup(t_now,name[0],sun_set.iso,date_start[0],Path)
		if telescope.find('Ganymede') is not -1:
			startup(t_now,name[0],sun_set.iso,date_start[0],Path)
		if telescope.find('Io') is not -1:
			startup(t_now,name[0],sun_set.iso,date_start[0],Path)
		if telescope.find('Callisto') is not -1:
			startup(t_now,name[0],sun_set.iso,date_start[0],Path)
		if telescope.find('Artemis') is not -1:
			startup_artemis(t_now,name[0],sun_set_teide.iso,date_start[0],Path)
		if telescope.find('Saint-Ex') is not -1:
			startup(t_now,name[0],sun_set_san_pedro.iso,date_start[0],Path)
		for i,nam in enumerate(name):
			print(i,len(name))
			if nam!='TransitionBlock':
				if (len(name)==2):
					if i==0:
						first_target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],name[i+1],Path)
					if i==0 and telescope.find('Ganymede') is not -1:
						first_target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],name[i+1],Path)
					# if i==1:
					# 	target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
					# 	flatdawn(t_now,date_end[i],sun_rise.iso,Path)
					if i==1 and telescope.find('Europa') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn(t_now,date_end[i],sun_rise.iso,Path)
					if i==1 and telescope.find('Callisto') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn(t_now,date_end[i],sun_rise.iso,Path)
					if i==1 and telescope.find('Io') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn(t_now,date_end[i],sun_rise.iso,Path)
					if i==1 and telescope.find('Ganymede') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn(t_now,date_end[i],sun_rise.iso,Path)
					if i==1 and telescope.find('Artemis') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i].replace('\'',''),texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn_artemis(t_now,date_end[i],sun_rise_teide.iso,Path)
					if i==1 and telescope.find('Saint-Ex') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn(t_now,date_end[i],sun_rise_san_pedro.iso,Path)
				else:
					if i==(len(name)-1) and telescope.find('Europa') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn(t_now,date_end[i],sun_rise.iso,Path)
					if i==(len(name)-1) and telescope.find('Callisto') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn(t_now,date_end[i],sun_rise.iso,Path)
					if i==(len(name)-1) and telescope.find('Io') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn(t_now,date_end[i],sun_rise.iso,Path)
					if i==(len(name)-1) and telescope.find('Ganymede') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn(t_now,date_end[i],sun_rise.iso,Path)
					if i==(len(name)-1) and telescope.find('Artemis') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i].replace('\'',''),texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn_artemis(t_now,date_end[i],sun_rise_teide.iso,Path)
					if i==(len(name)-1) and telescope.find('Saint-Ex') is not -1:
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path)
						flatdawn(t_now,date_end[i],sun_rise_san_pedro.iso,Path)
					if i<(len(name)-1):
						target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],name[i+1],Path)

		print('str(filt)',str(filt))
		if telescope.find('Callisto') is not -1:
			print('CALLISTO')
			flatexo_calli(Path,t_now,str(filt),nbu=3,nbB=3,nbz=3,nbV=3,nbr=3,nbi=3,nbg=3,nbIz=7,nbExo=3,nbClear=3)
		if telescope.find('Ganymede') is not -1:
			flatexo_gany(Path,t_now,str(filt),nbOIII=3,nbHa=3,nbSII=3,nbz=3,nbr=3,nbi=3,nbg=3,nbIz=7,nbExo=3,nbClear=3)
		if telescope.find('Io') is not -1:
			flatexo_io(Path,t_now,str(filt),nbu=3,nbHa=3,nbRc=3,nbz=3,nbr=3,nbi=3,nbg=3,nbIz=7,nbExo=3,nbClear=3)
		if telescope.find('Europa') is not -1:
			flatexo_euro(Path,t_now,str(filt),nbRc=3,nbB=3,nbz=3,nbV=3,nbr=3,nbi=3,nbg=3,nbIz=7,nbExo=3,nbClear=3)
		if telescope.find('Artemis') is not -1:
			flatexo_artemis_evening(Path,t_now,str(filt),nbu=3,nbz=3,nbr=3,nbi=3,nbg=3,nbIz=7,nbExo=3,nbClear=3)
			flatexo_artemis_morning(Path,t_now,str(filt),nbu=3,nbz=3,nbr=3,nbi=3,nbg=3,nbIz=7,nbExo=3,nbClear=3)

		if telescope.find('Saint-Ex') is not -1:
			flatexo_saintex(Path,t_now,str(filt),nbu=3,nbz=3,nbr=3,nbi=3,nbg=3,nbIz=7,nbExo=3,nbClear=3)

		if telescope.find('Artemis') is not -1:
			biasdark_comete(t_now,Path)
		else:
			biasdark(t_now,Path)

		p2=os.path.join('./DATABASE',str(telescope),str(t_now))
		shutil.make_archive(p2, 'zip', p)

