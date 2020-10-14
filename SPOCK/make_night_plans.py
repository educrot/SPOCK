#!/usr/bin/python
import os
import shutil
from astropy.table import Table
from astroplan import Observer,FixedTarget
from astropy.time import Time
from SPOCK.test_txtfiles import startup, startup_no_flats, Path_txt_files, flatexo_gany, flatexo_io, flatexo_euro, first_target_offset, flatexo_artemis_morning, flatexo_artemis_evening, startup_artemis,flatexo_saintex
from SPOCK.test_txtfiles import first_target,target, flatdawn, biasdark, shutdown, flatexo_calli, flatdawn_no_flats, target_no_DONUTS, target_offset, biasdark_comete, flatdawn_artemis
from astropy.coordinates import SkyCoord, get_sun, AltAz, EarthLocation
from astropy import units as u
import pandas as pd
pd.set_option('display.max_columns', 50)
import numpy as np
from astroplan.utils import time_grid_from_range

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


def make_scheduled_table(telescope,day_of_night):
    Path = './DATABASE'
    scheduled_table = None
    day_of_night = Time(day_of_night)
    try:
        os.path.exists(os.path.join(Path, telescope,
                                    'night_blocks_' + telescope + '_' + day_of_night.tt.datetime[0].strftime(
                                        "%Y-%m-%d") + '.txt'))
        print('INFO: Path exists and is: ', os.path.join(Path, telescope, 'night_blocks_' + telescope + '_' +
                                                         day_of_night.tt.datetime[0].strftime(
                                                             "%Y-%m-%d") + '.txt'))
    except TypeError:
        os.path.exists(os.path.join(Path, telescope,
                                    'night_blocks_' + telescope + '_' + day_of_night.tt.datetime.strftime(
                                        "%Y-%m-%d") + '.txt'))
        print('INFO: Path exists and is: ', os.path.join(Path, telescope,
                                                         'night_blocks_' + telescope + '_' + day_of_night.tt.datetime.strftime(
                                                             "%Y-%m-%d") + '.txt'))
    except NameError:
        print('INFO: no input night_block for this day')
    except FileNotFoundError:
        print('INFO: no input night_block for this day')

    if not (scheduled_table is None):
        return scheduled_table
    else:
        try:
            scheduled_table = Table.read(os.path.join(Path, telescope,
                                                           'night_blocks_' + telescope + '_' +
                                                           day_of_night.tt.datetime[0].strftime(
                                                               "%Y-%m-%d") + '.txt'), format='ascii')
            return scheduled_table
        except TypeError:
            scheduled_table = Table.read(os.path.join(Path, telescope,
                                                           'night_blocks_' + telescope + '_' + day_of_night.tt.datetime.strftime(
                                                               "%Y-%m-%d") + '.txt'), format='ascii')
            return scheduled_table


def dome_rotation(day_of_night,telescope):
    scheduled_table = make_scheduled_table(telescope,day_of_night)
    location = EarthLocation.from_geodetic(-70.40300000000002 * u.deg, -24.625199999999996 * u.deg,
                                           2635.0000000009704 * u.m)
    paranal = Observer(location=location, name="paranal", timezone="UTC")
    dur_dome_rotation = 2 / 60 / 24   # 5min
    number_of_targets = len(scheduled_table['target'])

    if number_of_targets == 1:
        old_end_time = scheduled_table['end time (UTC)'][0]

        start_dome_rot = Time((Time(scheduled_table['start time (UTC)'][0],format='iso').jd +
                              Time(scheduled_table['end time (UTC)'][0],format='iso').jd)/2,format='jd')

        end_dome_rot = Time(start_dome_rot.jd + dur_dome_rotation, format='jd')

        dur_first_block = Time(start_dome_rot.jd - Time(scheduled_table['start time (UTC)'][0],format='iso').jd,format='jd').jd
        dur_second_block = Time(Time(old_end_time).jd - end_dome_rot.jd,format='jd').jd

        coords = SkyCoord(str(int(scheduled_table['ra (h)'][0])) + 'h' + str(
            int(scheduled_table['ra (m)'][0])) + 'm' + str(
            round(scheduled_table['ra (s)'][0], 3)) + 's' + \
                          ' ' + str(int(scheduled_table['dec (d)'][0])) + 'd' + str(
            abs(int(scheduled_table['dec (m)'][0]))) + \
                          'm' + str(abs(round(scheduled_table['dec (s)'][0], 3))) + 's').transform_to(
            AltAz(obstime=start_dome_rot, location=paranal.location))
        coords_dome_rotation = SkyCoord(alt=coords.alt, az=(coords.az.value - 180) * u.deg, obstime=start_dome_rot,
                                        frame='altaz', location=paranal.location)
        if (coords.alt.value < 50):
            print('WARNING: not possible at that time because of altitude constraint, adding 20 degrees altitude')
            coords_dome_rotation = SkyCoord(alt=coords.alt + 20 * u.deg, az=(coords.az.value - 180) * u.deg,
                                            obstime=start_dome_rot, frame='altaz', location=paranal.location)

        target = FixedTarget(coord=SkyCoord(ra=coords_dome_rotation.icrs.ra.value * u.degree,
                                            dec=coords_dome_rotation.icrs.dec.value * u.degree),
                             name='dome_rot')

        scheduled_table.add_row([target.name, start_dome_rot.iso,  end_dome_rot.iso,
                           dur_dome_rotation * 24 * 60,
                            target.coord.ra.hms[0],
                            target.coord.ra.hms[1],  target.coord.ra.hms[2],
                            target.coord.dec.dms[0],  target.coord.dec.dms[1],
                            target.coord.dec.dms[2],  "{\'filt=I+z\', \'texp=10\'}"])

        scheduled_table['end time (UTC)'][0] = start_dome_rot.iso
        scheduled_table['duration (minutes)'][0] = dur_first_block* 24 * 60

        scheduled_table.add_row([scheduled_table['target'][0], end_dome_rot.iso,  old_end_time,
                           dur_second_block* 24 * 60,
                            scheduled_table['ra (h)'][0],
                            scheduled_table['ra (m)'][0],  scheduled_table['ra (s)'][0],
                            scheduled_table['dec (d)'][0],  scheduled_table['dec (m)'][0],
                            scheduled_table['dec (s)'][0],  scheduled_table['configuration'][0]])

        df = scheduled_table.to_pandas()
        df['target'][2] = df['target'][2] + '_2'

        scheduled_table = Table.from_pandas(df)
        scheduled_table.sort('start time (UTC)')


    if number_of_targets > 1:
        start_dome_rot = Time(scheduled_table['end time (UTC)'][0],format='iso')

        end_dome_rot = Time(start_dome_rot.jd + dur_dome_rotation, format='jd')

        coords = SkyCoord(str(int(scheduled_table['ra (h)'][0])) + 'h' + str(
            int(scheduled_table['ra (m)'][0])) + 'm' + str(
            round(scheduled_table['ra (s)'][0], 3)) + 's' + \
                          ' ' + str(int(scheduled_table['dec (d)'][0])) + 'd' + str(
            abs(int(scheduled_table['dec (m)'][0]))) + \
                          'm' + str(abs(round(scheduled_table['dec (s)'][0], 3))) + 's').transform_to(
            AltAz(obstime=start_dome_rot, location=paranal.location))
        coords_dome_rotation = SkyCoord(alt=coords.alt, az=(coords.az.value - 180) * u.deg, obstime=start_dome_rot,
                                        frame='altaz', location=paranal.location)
        if (coords.alt.value < 50):
            print('WARNING: not possible at that time because of altitude constraint, adding 20 degrees altitude')
            coords_dome_rotation = SkyCoord(alt=coords.alt + 20 * u.deg, az=(coords.az.value - 180) * u.deg,
                                            obstime=start_dome_rot,frame='altaz', location=paranal.location)

        target = FixedTarget(coord=SkyCoord(ra=coords_dome_rotation.icrs.ra.value * u.degree,
                                            dec=coords_dome_rotation.icrs.dec.value * u.degree),
                             name='dome_rot')

        scheduled_table.add_row([target.name, start_dome_rot.iso,  end_dome_rot.iso,
                           dur_dome_rotation * 24 * 60,
                            target.coord.ra.hms[0],
                            target.coord.ra.hms[1],  target.coord.ra.hms[2],
                            target.coord.dec.dms[0],  target.coord.dec.dms[1],
                            target.coord.dec.dms[2],  "{\'filt=I+z\', \'texp=10\'}"])

        scheduled_table['start time (UTC)'][1] = end_dome_rot.iso

        scheduled_table.sort('start time (UTC)')

    return scheduled_table


def make_np(t_now,nb_jours,tel):
    telescope=tel
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').tt.datetime.strftime("%Y-%m-%d")
        Path='./DATABASE'
        p=os.path.join(Path,str(telescope),'Plans_by_date',str(t_now))
        if not os.path.exists(p):
            os.makedirs(p)

        scheduler_table = Table.read('./DATABASE/' + str(telescope) +'/night_blocks_'+ str(telescope) +'_' + str(t_now)+'.txt', format='ascii')

        if (tel == 'Io') or tel == ('Europa') or (tel == 'Ganymede') or (tel == 'Callisto'):
            scheduler_table = dome_rotation(telescope=tel, day_of_night=t_now) # Intentional dome rotation to
            # avoid technical pb on Callisto with dome

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
            print()

        config_filled=config_mask#.filled()

        for i in range(0,len(scheduler_table)):
            if name[i]!='TransitionBlock':
                config=config_filled[i].split(',')
                config_tup=tuple(config)
                for item in config_tup:
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
                if filt[i]=='z' or filt[i]=='g' or filt[i]=='g' or filt[i]=='i' or filt[i]=='r':
                    a=filt[i]
                    filt[i]=a+'\''


        autofocus=None
        waitlimit=600
        afinterval=60
        count='5000'

        location = EarthLocation.from_geodetic(-70.40300000000002*u.deg, -24.625199999999996*u.deg,2635.0000000009704*u.m)
        paranal = Observer(location=location, name="paranal", timezone="UTC")
        t=Time(t_now)
        sun_set =paranal.sun_set_time(t,which='next')
        sun_rise =paranal.sun_rise_time(t,which='next')
        location_SNO = EarthLocation.from_geodetic(-16.50583131*u.deg, 28.2999988*u.deg, 2390*u.m)
        teide = Observer(location=location_SNO, name="SNO", timezone="UTC")
        sun_set_teide=teide.sun_set_time(t,which='next')
        sun_rise_teide=teide.sun_rise_time(t+1,which='next')
        location_saintex = EarthLocation.from_geodetic(-115.48694444444445*u.deg, 31.029166666666665*u.deg, 2829.9999999997976*u.m)
        san_pedro = Observer(location=location_saintex, name="saintex", timezone="UTC")
        sun_set_san_pedro=san_pedro.sun_set_time(t+1,which='next')
        sun_rise_san_pedro=san_pedro.sun_rise_time(t+1,which='next')

        Path=Path_txt_files(telescope)
        if telescope.find('Europa') is not -1:
            startup(t_now,name[0],sun_set.iso,date_start[0],Path,telescope)
        if telescope.find('Ganymede') is not -1:
            startup(t_now,name[0],sun_set.iso,date_start[0],Path,telescope)
        if telescope.find('Io') is not -1:
            startup(t_now,name[0],sun_set.iso,date_start[0],Path,telescope)
        if telescope.find('Callisto') is not -1:
            startup(t_now,name[0],sun_set.iso,date_start[0],Path,telescope)
        if telescope.find('Artemis') is not -1:
            startup_artemis(t_now,name[0],sun_set_teide.iso,date_start[0],Path)
        if telescope.find('Saint-Ex') is not -1:
            startup(t_now,name[0],sun_set_san_pedro.iso,date_start[0],Path,telescope)
        for i,nam in enumerate(name):
            if nam!='TransitionBlock':
                if (len(name)>=2):
                    if i==0:
                        first_target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],name[i+1],Path,telescope)
                    if i==0 and telescope.find('Ganymede') is not -1:
                        first_target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],name[i+1],Path,telescope)
                    if i==0 and telescope.find('Artemis') is not -1:
                        filt[i] = filt[i].replace('\'', '')
                        first_target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],name[i+1],Path,telescope='Artemis')
                    if i==0 and telescope.find('Saint-Ex') is not -1:
                        filt[i] = filt[i].replace('\'', '')
                        first_target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],name[i+1],Path,telescope='Saint-Ex')


                    if i==(len(name)-1) and telescope.find('Europa') is not -1:
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope)
                        flatdawn(t_now,date_end[i],sun_rise.iso,Path,telescope)
                    if i==(len(name)-1) and telescope.find('Callisto') is not -1:
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope)
                        flatdawn(t_now,date_end[i],sun_rise.iso,Path,telescope)
                    if i==(len(name)-1) and telescope.find('Io') is not -1:
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope)
                        flatdawn(t_now,date_end[i],sun_rise.iso,Path,telescope)
                    if i==(len(name)-1) and telescope.find('Ganymede') is not -1:
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope)
                        flatdawn(t_now,date_end[i],sun_rise.iso,Path,telescope)
                    if i==(len(name)-1) and telescope.find('Artemis') is not -1:
                        filt[i] = filt[i].replace('\'', '')
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,
                               filt[i].replace('\'',''),texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope='Artemis')
                        flatdawn_artemis(t_now,date_end[i],sun_rise_teide.iso,Path)
                    if i==(len(name)-1) and telescope.find('Saint-Ex') is not -1:
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope='Saint-Ex')
                        flatdawn(t_now,date_end[i],sun_rise_san_pedro.iso,Path,telescope)

                    if i<(len(name)-1):
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],name[i+1],Path,telescope=telescope)

                else:
                    if i==(len(name)-1) and telescope.find('Europa') is not -1:
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope)
                        flatdawn(t_now,date_end[i],sun_rise.iso,Path,telescope)
                    if i==(len(name)-1) and telescope.find('Callisto') is not -1:
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope)
                        flatdawn(t_now,date_end[i],sun_rise.iso,Path,telescope)
                    if i==(len(name)-1) and telescope.find('Io') is not -1:
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope)
                        flatdawn(t_now,date_end[i],sun_rise.iso,Path,telescope)
                    if i==(len(name)-1) and telescope.find('Ganymede') is not -1:
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope)
                        flatdawn(t_now,date_end[i],sun_rise.iso,Path,telescope)
                    if i==(len(name)-1) and telescope.find('Artemis') is not -1:
                        filt[i] = filt[i].replace('\'', '')
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i].replace('\'',''),texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope='Artemis')
                        flatdawn_artemis(t_now,date_end[i],sun_rise_teide.iso,Path)
                    if i==(len(name)-1) and telescope.find('Saint-Ex') is not -1:
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],None,Path,telescope='Saint-Ex')
                        flatdawn(t_now,date_end[i],sun_rise_san_pedro.iso,Path,telescope)
                    if i<(len(name)-1):
                        target(t_now,nam,date_start[i],date_end[i],waitlimit,afinterval, autofocus,count,filt[i],texp[i],ra1[i],ra2[i],ra3[i],dec1[i],dec2[i],dec3[i],name[i+1],Path,telescope=telescope)

        if telescope.find('Callisto') is not -1:
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
            flatexo_saintex(Path,t_now,str(filt),nbu=3,nbz=3,nbr=3,nbi=3,nbg=3,nbIz=9,nbExo=3,nbClear=3)

        biasdark(t_now,Path,telescope)
        p2=os.path.join('./DATABASE',str(telescope),'Zip_files',str(t_now))
        shutil.make_archive(p2, 'zip', p)

