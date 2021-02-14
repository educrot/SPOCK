#!/anaconda3/bin/python3.6
from astropy.time import Time
from astropy.table import unique,Table
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astroplan.utils import time_grid_from_range
from astroplan import FixedTarget, AltitudeConstraint, MoonSeparationConstraint,AtNightConstraint,AirmassConstraint,TimeConstraint
from astroplan.plots import dark_style_sheet,plot_airmass
from astroplan import Observer,moon_illumination
from astroplan.periodic import EclipsingSystem
from astroplan.constraints import is_event_observable
from colorama import Fore, Style
from docx import Document
from docx.shared import *
from datetime import datetime
from eScheduler.spe_schedule import SPECULOOSScheduler, Schedule, ObservingBlock,Transitioner
import matplotlib.pyplot as plt
import numpy as np
import os
import pkg_resources
import pandas as pd
import requests
import SPOCK.upload_night_plans as SPOCKunp #import upload_np_calli, upload_np_gany, upload_np_io, upload_np_euro,upload_np_artemis,upload_np_ts,upload_np_tn
from SPOCK.make_night_plans import make_np
import SPOCK.long_term_scheduler as SPOCKLT
import subprocess
import sys
import shutil
import SPOCK.ETC as ETC
from SPOCK import pwd_appcs,pwd_HUB,user_portal,pwd_portal,pwd_appcs,pwd_SNO_Reduc1,user_chart_studio,pwd_chart_studio,path_spock
import yaml

# pwd_appcs,pwd_HUB,user_portal,pwd_portal,pwd_appcs,pwd_SNO_Reduc1,user_chart_studio,pwd_chart_studio,path_spock = SPOCKLT._get_files()
dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day
constraints = [AltitudeConstraint(min=24*u.deg), AtNightConstraint()] #MoonSeparationConstraint(min=30*u.deg)

def charge_observatories(Name):
    observatories = []
    #Oservatories
    if 'SSO' in str(Name):
        location = EarthLocation.from_geodetic(-70.40300000000002*u.deg, -24.625199999999996*u.deg,2635.0000000009704*u.m)
        observatories.append(Observer(location=location, name="SSO", timezone="UTC"))

    if 'SNO' in str(Name):
        location_SNO = EarthLocation.from_geodetic(-16.50583131*u.deg, 28.2999988*u.deg, 2390*u.m)
        observatories.append(Observer(location=location_SNO, name="SNO", timezone="UTC"))

    if 'Saint-Ex' in str(Name):
        location_saintex = EarthLocation.from_geodetic(-115.48694444444445*u.deg, 31.029166666666665*u.deg, 2829.9999999997976*u.m)
        observatories.append(Observer(location=location_saintex, name="Saint-Ex", timezone="UTC"))

    if 'TS_La_Silla' in str(Name):
        location_TSlasilla = EarthLocation.from_geodetic(-70.73000000000002*u.deg, -29.25666666666666*u.deg, 2346.9999999988418*u.m)
        observatories.append(Observer(location=location_TSlasilla, name="TS_La_Silla", timezone="UTC"))

    if 'TN_Oukaimeden' in str(Name):
        location_TNOuka = EarthLocation.from_geodetic( -7.862263*u.deg, 31.20516*u.deg,2751*u.m)
        observatories.append(Observer(location=location_TNOuka, name="TN_Oukaimeden", timezone="UTC"))

    if 'Munich' in str(Name):
        location_munich= EarthLocation.from_geodetic(48.2*u.deg, -11.6*u.deg, 600*u.m)
        observatories.append(Observer(location=location_munich, name="Munich", timezone="UTC"))

    return observatories

def target_list_good_coord_format(path_target_list):

    """
    Give target corrdinates in ICRS format (used for astropy.coordinates SkyCoord function)

    Parameters
    ----------
    path_target_list: path on your computer toward the target list, by default take the one on the Cambridge server

    Returns
    -------
    targets: targets list with the following format : [<FixedTarget "Sp0002+0115" at SkyCoord (ICRS): (ra, dec) in deg (0.52591667, 1.26003889)>,


    """
    df = pd.read_csv(path_target_list, delimiter=' ')
    target_table_spc = Table.from_pandas(df)
    targets = [FixedTarget(coord=SkyCoord(ra=target_table_spc['RA'][i]* u.degree,dec=target_table_spc['DEC'][i] * u.degree),name=target_table_spc['Sp_ID'][i]) for i in range(len(target_table_spc['RA']))]
    return targets


class Schedules:

    def __init__(self):
        self.target_list = None
        self.telescopes = []
        self.telescope =  []
        self.start_end_range = None  # date_range
        self.day_of_night = None
        self.target_table_spc = []
        self.constraints = None
        self.SS1_night_blocks = None
        self.scheduled_table = None
        self.SS1_night_blocks_old  = None
        self.scheduled_table_sorted = None

    def load_parameters(self, input_file_short_term=None, nb_observatory=None,obs_name = None, start_end=None):
        if input_file_short_term == None:
            self.observatory = charge_observatories(obs_name)[0]
            self.target_list = path_spock + '/target_lists/target_list_special.txt'
            self.target_list_follow_up = path_spock + '/target_lists/target_transit_follow_up.txt'
            self.constraints = [AtNightConstraint.twilight_nautical()]
            df = pd.read_csv(self.target_list, delimiter=' ')
            df_followup = pd.read_csv(self.target_list_follow_up, delimiter=' ')
            self.target_table_spc = Table.from_pandas(df)
            self.target_table_spc_follow_up = Table.from_pandas(df_followup)
            self.targets = target_list_good_coord_format(self.target_list)
            self.targets_follow_up = target_list_good_coord_format(self.target_list_follow_up)
            if start_end is not  None:
                self.start_end_range = Time(start_end)
        else:
            with open(input_file_short_term, "r") as f:
                Inputs = yaml.load(f, Loader=yaml.FullLoader)
                self.target_list = path_spock + '/target_lists/target_list_special.txt'
                self.target_list_follow_up = path_spock + '/target_lists/target_transit_follow_up.txt'
                self.use = Inputs['use']
                df = pd.DataFrame.from_dict(Inputs['observatories'])
                self.observatory = charge_observatories(df[nb_observatory]['name'])[0]
                self.telescopes = df[nb_observatory]['telescopes']
                self.telescope = self.telescopes[0]
                self.day_of_night = Time(Inputs['day_of_night'])  # ,Time(Inputs['date_range'][1])]
                try:
                    self.day_of_night = self.day_of_night[0]
                except TypeError:
                    self.day_of_night = self.day_of_night
                self.start_end_range = Time(Inputs['start_end_range'])
                if self.start_end_range[1] <= self.start_end_range[0]:
                    sys.exit(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' end date inferior to start date')
                self.constraints = [AtNightConstraint.twilight_nautical()]
                df = pd.read_csv(self.target_list, delimiter=' ')
                df_followup = pd.read_csv(self.target_list_follow_up, delimiter=' ')
                self.target_table_spc = Table.from_pandas(df)
                self.target_table_spc_follow_up = Table.from_pandas(df_followup)
                self.targets = target_list_good_coord_format(self.target_list)
                self.targets_follow_up = target_list_good_coord_format(self.target_list_follow_up)

    def night_duration(self, day):
        '''

        :param day: day str format '%y%m%d HH:MM:SS.sss'
        :return:
        '''
        dt_1day = Time('2018-01-02 00:00:00', scale='tcg') - Time('2018-01-01 00:00:00', scale='tcg')
        return ((Time(self.observatory.twilight_morning_nautical(day + dt_1day, which='nearest'))) \
                - (Time(self.observatory.twilight_evening_nautical(day, which='next'))))

    def monitoring(self,input_name,airmass_max,time_monitoring):
        idx_first_target = int(np.where((self.target_table_spc['Sp_ID'] == input_name))[0])
        dur_mon_target = time_monitoring * u.minute
        constraints_monitoring_target = [AltitudeConstraint(min=25 * u.deg), MoonSeparationConstraint(min=25 * u.deg),\
                                         AirmassConstraint(max = airmass_max,boolean_constraint=True), \
                                         TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.day_of_night, which='next'))), \
                                                        (Time(self.observatory.twilight_morning_nautical(self.day_of_night+1, which='nearest'))))]
        if self.target_table_spc['texp_spc'][idx_first_target] == 0:
            self.target_table_spc['texp_spc'][idx_first_target]= self.exposure_time(input_name=self.target_table_spc['Sp_ID'][idx_first_target])
        blocks = []
        a = ObservingBlock(self.targets[idx_first_target], dur_mon_target, -1, constraints=constraints_monitoring_target,\
                           configuration={'filt=' + str(self.target_table_spc['Filter_spc'][idx_first_target]),'texp=' + str(self.target_table_spc['texp_spc'][idx_first_target])})
        blocks.append(a)
        transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
        seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night+1)
        sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_monitoring_target, observer=self.observatory,transitioner=transitioner)
        sequen_scheduler_SS1(blocks, seq_schedule_SS1)
        self.SS1_night_blocks = seq_schedule_SS1.to_table()  # Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
        return self.SS1_night_blocks

    def dome_rotation(self):
         # not necessary
        dur_dome_rotation = 5 / 60 / 24 * u.day  # 5min

        sun_set = Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next'))
        sun_rise = Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest'))

        time_resolution = 1*u.hour
        time_grid = time_grid_from_range([sun_set, sun_rise], time_resolution=time_resolution)
        dom_rot_possible = False

        self.make_scheduled_table()
        end_times = Time(self.scheduled_table['end time (UTC)'][:]).iso

        for i in range(1,len(time_grid)-1):

            start = time_grid[i]
            print(start.iso)
            end = start + dur_dome_rotation
            idx = np.where((start < end_times[:]))[0]

            coords = SkyCoord(str(int(self.scheduled_table['ra (h)'][idx][0])) +  'h' +   str(int(self.scheduled_table['ra (m)'][idx][0])) + 'm' + str(round(self.scheduled_table['ra (s)'][idx][0],3)) + 's' +\
                                  ' ' + str(int(self.scheduled_table['dec (d)'][idx][0])) +  'd' +   str(abs(int(self.scheduled_table['dec (m)'][idx][0]))) +\
                              'm' + str(abs(round(self.scheduled_table['dec (s)'][idx][0],3))) + 's').transform_to(AltAz(obstime=start,location=self.observatory.location))
            coords_dome_rotation = SkyCoord(alt=coords.alt, az=(coords.az.value - 180) * u.deg, obstime=start,
                                            frame='altaz', location=self.observatory.location)
            if (coords.alt.value > 60) or (coords.alt.value < 30):
                dom_rot_possible = False
                print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' not possible at that time because of altitude constraint')# = ' + str(round(coords_dome_rotation.alt.value,3)))
                #sys.exit(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' No dome rotation possible at that time, altitude = ' + str(coords.alt.value) + '°  (altitude <30° or altitude >60°)')

            else:
                target = FixedTarget(coord=SkyCoord(ra= coords_dome_rotation.icrs.ra.value * u.degree, dec= coords_dome_rotation.icrs.dec.value * u.degree),
                        name='dome_rot')
                dom_rot_possible = True

                df = pd.DataFrame({'target': target.name, 'start time (UTC)': start.iso, 'end time (UTC)': end.iso,
                              'duration (minutes)': dur_dome_rotation.value * 24 * 60, 'ra (h)': target.coord.ra.hms[0],
                              'ra (m)': target.coord.ra.hms[1], 'ra (s)': target.coord.ra.hms[2],
                              'dec (d)': target.coord.dec.dms[0], 'dec (m)': target.coord.dec.dms[1],
                              'dec (s)': target.coord.dec.dms[2], 'configuration': "{\'filt=I+z\', \'texp=10\'}"}, index=[0])
                self.SS1_night_blocks = Table.from_pandas(df)
                return self.SS1_night_blocks

                break

        if dom_rot_possible == False:
            sys.exit(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' No Dom rotation possible that night')

    def special_target_with_start_end(self, input_name):
        start = self.start_end_range[0]
        end = self.start_end_range[1]

        if (start >= Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso) or (end <= Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next')).iso):
            sys.exit('WARNING : Start time (or End time) is not on the same day.')

        dur_obs_both_target = (self.night_duration(self.day_of_night) / (2 * u.day)) * 2 * u.day
        constraints_special_target = [AltitudeConstraint(min=25 * u.deg), MoonSeparationConstraint(min=25 * u.deg),
                                      TimeConstraint(start, end)]
        idx_to_insert_target = int(np.where((self.target_table_spc['Sp_ID'] == input_name))[0])
        if self.target_table_spc['texp_spc'][idx_to_insert_target] == 0 \
                or self.target_table_spc['texp_spc'][idx_to_insert_target] == "00":
            self.target_table_spc['texp_spc'][idx_to_insert_target] = \
                int(self.exposure_time(input_name=self.target_table_spc['Sp_ID'][idx_to_insert_target],
                                       target_list=self.target_table_spc))
        blocks = []
        a = ObservingBlock(self.targets[idx_to_insert_target], dur_obs_both_target, -1,
                           constraints=constraints_special_target,
                           configuration={'filt=' + str(self.target_table_spc['Filter_spc'][idx_to_insert_target]),
                                          'texp=' + str(self.target_table_spc['texp_spc'][idx_to_insert_target])})
        blocks.append(a)
        transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
        seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night + 1)
        sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_special_target, observer=self.observatory,
                                                  transitioner=transitioner)
        sequen_scheduler_SS1(blocks, seq_schedule_SS1)
        self.SS1_night_blocks = seq_schedule_SS1.to_table()  # Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
        return self.SS1_night_blocks

    def special_target(self,input_name):
        dur_obs_both_target=(self.night_duration(self.day_of_night)/(2*u.day))*2*u.day
        constraints_special_target=[AltitudeConstraint(min=24*u.deg),MoonSeparationConstraint(min=30*u.deg), \
        TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next'))), \
                       (Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest'))))]
        idx_first_target = int(np.where((self.target_table_spc['Sp_ID']==input_name))[0])
        if int(self.target_table_spc['texp_spc'][idx_first_target]) == 0:
            self.target_table_spc['texp_spc'][idx_first_target] = \
                self.exposure_time(input_name=self.target_table_spc['Sp_ID'][idx_first_target],
                                   target_list=self.target_table_spc)
        blocks=[]
        a = ObservingBlock(self.targets[idx_first_target],dur_obs_both_target,-1,constraints= constraints_special_target,\
                         configuration={'filt=' + str(self.target_table_spc['Filter_spc'][idx_first_target]),
                                        'texp=' + str(self.target_table_spc['texp_spc'][idx_first_target])})
        blocks.append(a)
        transitioner = Transitioner(slew_rate= 11*u.deg/u.second)
        seq_schedule_SS1 = Schedule(self.day_of_night,self.day_of_night+1)
        sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_special_target, observer = self.observatory,transitioner=transitioner)
        sequen_scheduler_SS1(blocks,seq_schedule_SS1)
        self.SS1_night_blocks=seq_schedule_SS1.to_table()
        if len(self.SS1_night_blocks) == 0:
            print('WARNING : Impossible to schedule target ' + input_name + ' at this time range and/or with those constraints')
        return self.SS1_night_blocks

    def transit_follow_up(self,follow_up_list=None,name=None):
        self.targets = self.targets_follow_up
        if follow_up_list is None:
            follow_up_list = self.target_list_follow_up

        df = pd.read_csv(follow_up_list,delimiter= ' ')
        constraints = [AltitudeConstraint(min=25 * u.deg), AtNightConstraint(),
                       TimeConstraint(Time(self.day_of_night),Time(self.day_of_night+1))]

        if name is not None:
            i = int(np.where((df['Sp_ID']==name))[0])
            blocks = []
            epoch = Time(df['T0'][i], format='jd')
            period = df['P'][i] * u.day
            duration = df['W'][i] * u.day
            oot_time = duration.value * 1.5 * u.day
            T0_err_transit = df['T0_err'][i]
            P_err_transit = df['P_err'][i]
            W_err_transit = df['W_err'][i]
            target_transit = EclipsingSystem(primary_eclipse_time=epoch, orbital_period=period, duration=duration,
                                             name=df['Sp_ID'][i])
            print(Fore.GREEN + Fore.GREEN + Fore.GREEN + 'INFO: ' + Fore.BLACK + ' ' + Fore.BLACK + Fore.BLACK + str(df['Sp_ID'][i]) + ' next transit: ')
            print(Time(target_transit.next_primary_eclipse_time(self.day_of_night, n_eclipses=1)).iso)
            timing_to_obs_jd = Time(target_transit.next_primary_eclipse_time(self.day_of_night, n_eclipses=1)).jd
            target = target_list_good_coord_format(follow_up_list)
            n_transits = 1
            try:
                ing_egr = target_transit.next_primary_ingress_egress_time(self.day_of_night, n_eclipses=n_transits)
            except ValueError:
                print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' No transit of ', df['Sp_ID'][i], ' on the period chosen')

            observable = is_event_observable(constraints, self.observatory, target, times_ingress_egress=ing_egr)
            if np.any(observable):
                err_T0_neg = timing_to_obs_jd[0] - (np.round((timing_to_obs_jd[0] - epoch.jd) / period.value, 1) *
                                                    (period.value - P_err_transit) + (epoch.jd - T0_err_transit))
                err_T0_pos = (np.round((timing_to_obs_jd[0] - epoch.jd) / period.value, 1) *
                              (period.value + P_err_transit) + (epoch.jd + T0_err_transit)) - timing_to_obs_jd[0]
                start_transit = Time(ing_egr[0][0].value - err_T0_neg - oot_time.value - W_err_transit,
                                     format='jd')  # - T0_err_transit - W_err_transit  - oot_time.value/24 - (timing_to_obs_jd[0] - epoch.jd) / period.value * P_err_transit,format='jd')
                # print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit', start_transit.iso)
                end_transit = Time(ing_egr[0][1].value + err_T0_pos + oot_time.value + W_err_transit,
                                   format='jd')  # + T0_err_transit + W_err_transit + oot_time.value/24 + (timing_to_obs_jd[0] - epoch.jd) / period.value * P_err_transit ,format='jd')
                # print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit', end_transit.iso)
                dur_obs_transit_target = (end_transit - start_transit).value * 1. * u.day
                if (end_transit > Time(
                        self.observatory.twilight_morning_nautical(self.day_of_night + 1, which='nearest'))) \
                        or (start_transit < Time(
                    self.observatory.twilight_evening_nautical(self.day_of_night, which='next'))):
                    if (Time(ing_egr[0][1]) < Time(
                            self.observatory.twilight_morning_nautical(self.day_of_night + 1, which='nearest'))) \
                            and (Time(ing_egr[0][0]) > Time(
                        self.observatory.twilight_evening_nautical(self.day_of_night, which='next'))):
                        constraints_transit_target = [AltitudeConstraint(min=25 * u.deg),
                                                      TimeConstraint(start_transit, end_transit),
                                                      MoonSeparationConstraint(min=25 * u.deg),
                                                      AtNightConstraint(max_solar_altitude=-12 * u.deg)]
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ', str(df['Sp_ID'][i]), ' : ', start_transit.iso)
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit of ', str(df['Sp_ID'][i]), ' : ', end_transit.iso)
                        print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Time out of transit not optimal.')
                        idx_first_target = int(np.where((df['Sp_ID'] == df['Sp_ID'][i]))[0])
                        if df['texp_spc'][idx_first_target] == 0:
                            df['texp_spc'][idx_first_target] = self.exposure_time(
                                input_name=df['Sp_ID'][idx_first_target], target_list=self.target_table_spc_follow_up)
                        a = ObservingBlock(target[idx_first_target], dur_obs_transit_target, -1,
                                           constraints=constraints_transit_target,
                                           configuration={'filt=' + str(df['Filter_spc'][idx_first_target]),
                                                          'texp=' + str(df['texp_spc'][idx_first_target])})
                        blocks.append(a)
                        transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
                        seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night + 1)
                        sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_transit_target,
                                                                  observer=self.observatory,
                                                                  transitioner=transitioner)
                        sequen_scheduler_SS1(blocks, seq_schedule_SS1)

                    else:
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ', str(df['Sp_ID'][i]), ' : ', Time(ing_egr[0][0]).iso)
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit of ', str(df['Sp_ID'][i]), ' : ', Time(ing_egr[0][1]).iso)
                        print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Transit not full.')
                        constraints_transit_target = [AltitudeConstraint(min=25 * u.deg),
                                                      MoonSeparationConstraint(min=25 * u.deg),
                                                      TimeConstraint(Time(ing_egr[0][0]), Time(ing_egr[0][1]))]
                        # continue
                        idx_first_target = int(np.where((df['Sp_ID'] == df['Sp_ID'][i]))[0])
                        if df['texp_spc'][idx_first_target] == 0:
                            df['texp_spc'][idx_first_target] = self.exposure_time(
                                input_name=df['Sp_ID'][idx_first_target], target_list=self.target_table_spc_follow_up)
                        a = ObservingBlock(target[idx_first_target], dur_obs_transit_target, -1,
                                           constraints=constraints_transit_target,
                                           configuration={'filt=' + str(df['Filter_spc'][idx_first_target]),
                                                          'texp=' + str(df['texp_spc'][idx_first_target])})
                        blocks.append(a)
                        transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
                        seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night + 1)
                        sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_transit_target,
                                                                  observer=self.observatory,
                                                                  transitioner=transitioner)
                        sequen_scheduler_SS1(blocks, seq_schedule_SS1)

                else:
                    constraints_transit_target = [AltitudeConstraint(min=25 * u.deg),
                                                  MoonSeparationConstraint(min=25 * u.deg),
                                                  TimeConstraint(start_transit, end_transit)]
                    idx_first_target = int(np.where((df['Sp_ID'] == df['Sp_ID'][i]))[0])
                    if df['texp_spc'][idx_first_target] == 0:
                        df['texp_spc'][idx_first_target] = self.exposure_time(input_name=df['Sp_ID'][idx_first_target],
                                                                              target_list=self.target_table_spc_follow_up)
                    a = ObservingBlock(target[idx_first_target], dur_obs_transit_target, -1,
                                       constraints=constraints_transit_target,
                                       configuration={'filt=' + str(df['Filter_spc'][idx_first_target]),
                                                      'texp=' + str(df['texp_spc'][idx_first_target])})
                    blocks.append(a)
                    transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
                    seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night + 1)
                    sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_transit_target,
                                                              observer=self.observatory,
                                                              transitioner=transitioner)
                    sequen_scheduler_SS1(blocks, seq_schedule_SS1)
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ', str(df['Sp_ID'][i]), ' : ', Time(ing_egr[0][0]).iso)
                    print(Fore.GREEN + Fore.GREEN + 'INFO: ' + Fore.BLACK + ' ' + Fore.BLACK + ' end_transit of ', str(df['Sp_ID'][i]), ' : ', Time(ing_egr[0][1]).iso)
                    print(Fore.GREEN + Fore.GREEN + 'INFO: ' + Fore.BLACK + ' ' + Fore.BLACK + ' Transit is expected to be full.')
                    if len(seq_schedule_SS1.to_table()['target']) == 0:
                        print(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' The moon is too closed for the transit to be observed')

                if len(seq_schedule_SS1.to_table()['target']) != 0:
                    self.SS1_night_blocks = seq_schedule_SS1.to_table()  # Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
                    return self.SS1_night_blocks
            else:
                print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' no transit of ', df['Sp_ID'][i], ' this day')

        else:
            for i in range(len(df['Sp_ID'])):
                blocks = []
                epoch = Time(df['T0'][i], format='jd')
                period = df['P'][i] * u.day
                duration = df['W'][i] * u.day
                oot_time = duration.value * 1.5 * u.day
                T0_err_transit = df['T0_err'][i]
                P_err_transit = df['P_err'][i]
                W_err_transit = df['W_err'][i]
                target_transit = EclipsingSystem(primary_eclipse_time = epoch, orbital_period=period, duration=duration, name=df['Sp_ID'][i])
                print(Fore.GREEN + Fore.GREEN + Fore.GREEN + 'INFO: ' + Fore.BLACK + ' ' + Fore.BLACK + Fore.BLACK + str(df['Sp_ID'][i]) + ' next transit: ', Time(target_transit.next_primary_eclipse_time(self.day_of_night, n_eclipses=1)).iso)
                timing_to_obs_jd = Time(target_transit.next_primary_eclipse_time(self.day_of_night, n_eclipses=1)).jd
                target = target_list_good_coord_format(follow_up_list)
                n_transits = 1
                try:
                    ing_egr = target_transit.next_primary_ingress_egress_time(self.day_of_night,n_eclipses=n_transits)
                except ValueError:
                    print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' No transit of ', df['Sp_ID'][i], ' on the period chosen')

                observable = is_event_observable(constraints,self.observatory, target, times_ingress_egress=ing_egr)
                if np.any(observable):
                    err_T0_neg = timing_to_obs_jd[0] - (np.round((timing_to_obs_jd[0] - epoch.jd) / period.value,1) *
                                                        (period.value - P_err_transit) + (epoch.jd - T0_err_transit ))
                    err_T0_pos = (np.round((timing_to_obs_jd[0] - epoch.jd) / period.value,1) *
                                  (period.value + P_err_transit) + (epoch.jd + T0_err_transit )) - timing_to_obs_jd[0]
                    start_transit = Time(ing_egr[0][0].value - err_T0_neg - oot_time.value - W_err_transit,format='jd') #- T0_err_transit - W_err_transit  - oot_time.value/24 - (timing_to_obs_jd[0] - epoch.jd) / period.value * P_err_transit,format='jd')
                    #print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit', start_transit.iso)
                    end_transit = Time(ing_egr[0][1].value + err_T0_pos + oot_time.value + W_err_transit,format='jd') #+ T0_err_transit + W_err_transit + oot_time.value/24 + (timing_to_obs_jd[0] - epoch.jd) / period.value * P_err_transit ,format='jd')
                    #print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit', end_transit.iso)
                    dur_obs_transit_target = (end_transit - start_transit).value * 1. * u.day
                    if (end_transit > Time(self.observatory.twilight_morning_nautical(self.day_of_night + 1, which='nearest'))) \
                            or (start_transit < Time(self.observatory.twilight_evening_nautical(self.day_of_night, which='next'))):
                        if (Time(ing_egr[0][1]) < Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')))\
                                and (Time(ing_egr[0][0]) > Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next'))):
                            constraints_transit_target = [AltitudeConstraint(min=25 * u.deg),TimeConstraint(start_transit, end_transit),
                                                          MoonSeparationConstraint(min=25 * u.deg),AtNightConstraint(max_solar_altitude=-12*u.deg)]
                            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ', str(df['Sp_ID'][i]), ' : ', start_transit.iso)
                            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit of ', str(df['Sp_ID'][i]), ' : ', end_transit.iso)
                            print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Time out of transit not optimal.')
                            idx_first_target = int(np.where((df['Sp_ID'] == df['Sp_ID'][i]))[0])
                            if df['texp_spc'][idx_first_target] == 0:
                                df['texp_spc'][idx_first_target] = self.exposure_time(
                                    input_name=df['Sp_ID'][idx_first_target],target_list=self.target_table_spc_follow_up)
                            a = ObservingBlock(target[idx_first_target], dur_obs_transit_target, -1,
                                               constraints=constraints_transit_target,
                                               configuration={'filt=' + str(df['Filter_spc'][idx_first_target]),
                                                              'texp=' + str(df['texp_spc'][idx_first_target])})
                            blocks.append(a)
                            transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
                            seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night + 1)
                            sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_transit_target,
                                                                      observer=self.observatory,
                                                                      transitioner=transitioner)
                            sequen_scheduler_SS1(blocks, seq_schedule_SS1)

                        else:
                            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ', str(df['Sp_ID'][i]), ' : ', Time(ing_egr[0][0]).iso)
                            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit of ', str(df['Sp_ID'][i]), ' : ', Time(ing_egr[0][1]).iso)
                            print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Transit not full.')
                            constraints_transit_target = [AltitudeConstraint(min=25 * u.deg),
                                                          MoonSeparationConstraint(min=25 * u.deg),
                                                          TimeConstraint(Time(ing_egr[0][0]),Time(ing_egr[0][1]))]
                            #continue
                            idx_first_target = int(np.where((df['Sp_ID'] == df['Sp_ID'][i]))[0])
                            if df['texp_spc'][idx_first_target] == 0:
                                df['texp_spc'][idx_first_target] = self.exposure_time(
                                    input_name=df['Sp_ID'][idx_first_target],target_list=self.target_table_spc_follow_up)
                            a = ObservingBlock(target[idx_first_target], dur_obs_transit_target, -1,
                                               constraints=constraints_transit_target,
                                               configuration={'filt=' + str(df['Filter_spc'][idx_first_target]),
                                                              'texp=' + str(df['texp_spc'][idx_first_target])})
                            blocks.append(a)
                            transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
                            seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night + 1)
                            sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_transit_target,
                                                                      observer=self.observatory,
                                                                      transitioner=transitioner)
                            sequen_scheduler_SS1(blocks, seq_schedule_SS1)

                    else:
                        constraints_transit_target = [AltitudeConstraint(min=25 * u.deg),
                                                      MoonSeparationConstraint(min=25 * u.deg),
                                                      TimeConstraint(start_transit, end_transit)]
                        idx_first_target = int(np.where((df['Sp_ID'] == df['Sp_ID'][i]))[0])
                        if df['texp_spc'][idx_first_target] == 0:
                            df['texp_spc'][idx_first_target] = self.exposure_time(input_name=df['Sp_ID'][idx_first_target],
                                                                                  target_list=self.target_table_spc_follow_up)
                        a = ObservingBlock(target[idx_first_target], dur_obs_transit_target, -1,
                                           constraints=constraints_transit_target,
                                           configuration={'filt=' + str(df['Filter_spc'][idx_first_target]),
                                                          'texp=' + str(df['texp_spc'][idx_first_target])})
                        blocks.append(a)
                        transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
                        seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night+1)
                        sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_transit_target, observer=self.observatory,
                                                                  transitioner=transitioner)
                        sequen_scheduler_SS1(blocks, seq_schedule_SS1)
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ',str(df['Sp_ID'][i]), ' : ',Time(ing_egr[0][0]).iso)
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit of ',str(df['Sp_ID'][i]), ' : ',Time(ing_egr[0][1]).iso)
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' Transit is expected to be full.')
                        if len(seq_schedule_SS1.to_table()['target']) == 0:
                            print(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' The moon is too closed for the transit to be observed')

                    if len(seq_schedule_SS1.to_table()['target']) != 0:
                        self.SS1_night_blocks = seq_schedule_SS1.to_table()  # Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
                        return self.SS1_night_blocks
                else:
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' no transit of ', df['Sp_ID'][i],' this day')

    def planification(self):
        end_scheduled_table = pd.DataFrame(columns=['target','start time (UTC)','end time (UTC)','duration (minutes)','ra (h)','ra (m)','ra (s)',\
                                                    'dec (d)','dec (m)','dec (s)','configuration'])
        try:
            self.SS1_night_blocks['target'][0]
        except TypeError:
            sys.exit(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' No block to insert ')
        if self.SS1_night_blocks['target'][0] in self.scheduled_table['target']:
            sys.exit(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' This target is already scheduled this day')
        for i in range(len(self.scheduled_table['target'])):
            print(i)
            try:
                self.SS1_night_blocks['target'][0]
            except IndexError:
                sys.exit(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' No block to insert ')
            except TypeError:
                sys.exit(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' No block to insert ')
            #if self.SS1_night_blocks['target'][0] == self.scheduled_table['target'][i]:
            #    sys.exit(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' The target ' + str(self.SS1_night_blocks['target'][0]) + 'that you wich to insert is already scheduled for this night !')
            end_beFore_cut = self.scheduled_table['end time (UTC)'][i]
            start_beFore_cut = self.scheduled_table['start time (UTC)'][i]

            if self.SS1_night_blocks is None :
                sys.exit('WARNING : No block to insert !')

            if not (self.scheduled_table_sorted is None):
                print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' Several transits this night')
                if (np.any(str(self.SS1_night_blocks['target'][0]).find('Trappist-1'))==0) and (np.any(str(self.scheduled_table_sorted['target'][0]).find('Trappist-1'))==0):
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' two TRAPPIST-1 planets this night!')
                    self.SS1_night_blocks['start time (UTC)'][0] = min(self.SS1_night_blocks['start time (UTC)'][0],self.scheduled_table_sorted['start time (UTC)'][0])
                    self.SS1_night_blocks['end time (UTC)'][0] = max(self.SS1_night_blocks['end time (UTC)'][0],self.scheduled_table_sorted['end time (UTC)'][0])
                    self.SS1_night_blocks['duration (minutes)'][0] = (Time(self.SS1_night_blocks['end time (UTC)'][0]) - Time(self.SS1_night_blocks['start time (UTC)'][0])).value*24*60
                    end_scheduled_table.append({'target':self.SS1_night_blocks['target'][0],\
                                                'start time (UTC)':self.SS1_night_blocks['start time (UTC)'][0],\
                                                'end time (UTC)':self.SS1_night_blocks['end time (UTC)'][0],\
                                                'duration (minutes)':self.SS1_night_blocks['duration (minutes)'][0],\
                                                'ra (h)':self.SS1_night_blocks['ra (h)'][0],\
                                                'ra (m)':self.SS1_night_blocks['ra (m)'][0],\
                                                'ra (s)':self.SS1_night_blocks['ra (s)'][0],\
                                                'dec (d)':self.SS1_night_blocks['dec (d)'][0],\
                                                'dec (m)':self.SS1_night_blocks['dec (m)'][0],\
                                                'dec (s)':self.SS1_night_blocks['dec (s)'][0],\
                                                'configuration':self.SS1_night_blocks['configuration'][0]},ignore_index=True)
            #situation 1
            if (self.SS1_night_blocks['start time (UTC)'][0] <= Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next')).iso):
                #case 1 # if new_block starts beFore twilight sun set
                print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' situation 1')
                self.SS1_night_blocks['start time (UTC)'][0] = Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next')).iso

            if (self.SS1_night_blocks['start time (UTC)'][0] < start_beFore_cut) and \
                    (self.SS1_night_blocks['start time (UTC)'][0] < end_beFore_cut): # if new_block starts beFore exiting_block starts AND beFore existing_block ends
                # situation 2
                if (self.SS1_night_blocks['end time (UTC)'][0] > start_beFore_cut) and (self.SS1_night_blocks['end time (UTC)'][0] < end_beFore_cut): #case 2 # if new_block ends beFore the end of existing block
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' situation 2')
                    self.scheduled_table['start time (UTC)'][i] = self.SS1_night_blocks['end time (UTC)'][0]
                    self.scheduled_table['duration (minutes)'][i] = (Time(self.scheduled_table['end time (UTC)'][i])-Time(self.scheduled_table['start time (UTC)'][i])).value*24*60

                # situation 3
                if (self.SS1_night_blocks['end time (UTC)'][0] <= start_beFore_cut):
                    self.scheduled_table['duration (minutes)'][i] = (Time(self.scheduled_table['end time (UTC)'][i])-Time(self.scheduled_table['start time (UTC)'][i])).value*24*60
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' situation 3, no change made to initial schedule')

                # situation 4
                elif (self.SS1_night_blocks['end time (UTC)'][0] >= end_beFore_cut) and \
                        (self.SS1_night_blocks['end time (UTC)'][0] <= Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso) :
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' situation 4')
                    self.scheduled_table[i] = self.SS1_night_blocks[0] #a way  to erase self.scheduled_table block

                # situation 5
                elif (self.SS1_night_blocks['end time (UTC)'][0] >= Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso) :
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' situation 5')
                    self.SS1_night_blocks['end time (UTC)'][0] = Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso
                    self.scheduled_table[i] = self.SS1_night_blocks[0] #a way  to erase self.scheduled_table block

            if (self.SS1_night_blocks['start time (UTC)'][0] >= start_beFore_cut):

                # situation 6
                if (self.SS1_night_blocks['start time (UTC)'][0] <= end_beFore_cut) and (self.SS1_night_blocks['end time (UTC)'][0] <= end_beFore_cut):
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' situation 6')
                    self.scheduled_table['end time (UTC)'][i] =  self.SS1_night_blocks['start time (UTC)'][0]
                    self.scheduled_table['duration (minutes)'] = (Time(self.scheduled_table['end time (UTC)'][i]) - Time(self.scheduled_table['start time (UTC)'][i])).value * 24*60

                    end_scheduled_table = end_scheduled_table.append({'target': self.scheduled_table['target'][i] + '_2', \
                                                                      'start time (UTC)':
                                                                          self.SS1_night_blocks['end time (UTC)'][0], \
                                                                      'end time (UTC)': end_beFore_cut, \
                                                                      'duration (minutes)': (Time(end_beFore_cut) - Time(
                                                                          self.SS1_night_blocks['end time (UTC)'][
                                                                              0])).value * 24 * 60, \
                                                                      'ra (h)': self.scheduled_table['ra (h)'][i], \
                                                                      'ra (m)': self.scheduled_table['ra (m)'][i], \
                                                                      'ra (s)': self.scheduled_table['ra (s)'][i], \
                                                                      'dec (d)': self.scheduled_table['dec (d)'][i], \
                                                                      'dec (m)': self.scheduled_table['dec (m)'][i], \
                                                                      'dec (s)': self.scheduled_table['dec (s)'][i], \
                                                                      'configuration':
                                                                          self.scheduled_table['configuration'][i]},
                                                                     ignore_index=True)

                if (self.SS1_night_blocks['end time (UTC)'][0] >= end_beFore_cut):

                    # situation 7
                    if (self.SS1_night_blocks['end time (UTC)'][0] >= Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso):
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' situation 7')
                        self.scheduled_table['end time (UTC)'][i] = self.SS1_night_blocks['start time (UTC)'][0]
                        self.scheduled_table['duration (minutes)'] = (Time(self.scheduled_table['end time (UTC)'][i]) - Time(self.scheduled_table['start time (UTC)'][i])).value * 24 * 60

                        self.SS1_night_blocks['end time (UTC)'][0] = Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso
                        self.SS1_night_blocks['duration (minutes)'][0] = (Time(self.SS1_night_blocks['end time (UTC)'][0])  - Time(self.SS1_night_blocks['start time (UTC)'][0])).value * 24*60

                    elif (self.SS1_night_blocks['end time (UTC)'][0] <= Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso):
                        # situation 8
                        if (self.SS1_night_blocks['start time (UTC)'][0] >= end_beFore_cut):
                            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' situation 8, no change made to initial schedule')
                        # situation 9
                        else:
                            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' situation 9')
                            self.scheduled_table['end time (UTC)'][i] = self.SS1_night_blocks['start time (UTC)'][0]
                            self.scheduled_table['duration (minutes)'] = (Time(self.scheduled_table['end time (UTC)'][i]) - Time(self.scheduled_table['start time (UTC)'][i])).value * 24 * 60

                # situation 10
                if (self.SS1_night_blocks['start time (UTC)'][0] >= end_beFore_cut):
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' situation 10, no change made to initial schedule')

            end_scheduled_table = end_scheduled_table.append({'target': self.SS1_night_blocks['target'][0], \
                                        'start time (UTC)': self.SS1_night_blocks['start time (UTC)'][0], \
                                        'end time (UTC)': self.SS1_night_blocks['end time (UTC)'][0], \
                                        'duration (minutes)': self.SS1_night_blocks['duration (minutes)'][0], \
                                        'ra (h)': self.SS1_night_blocks['ra (h)'][0], \
                                        'ra (m)': self.SS1_night_blocks['ra (m)'][0], \
                                        'ra (s)': self.SS1_night_blocks['ra (s)'][0], \
                                        'dec (d)': self.SS1_night_blocks['dec (d)'][0], \
                                        'dec (m)': self.SS1_night_blocks['dec (m)'][0], \
                                        'dec (s)': self.SS1_night_blocks['dec (s)'][0], \
                                        'configuration': self.SS1_night_blocks['configuration'][0]},
                                       ignore_index=True)
            end_scheduled_table = end_scheduled_table.append({'target': self.scheduled_table['target'][i], \
                                        'start time (UTC)': self.scheduled_table['start time (UTC)'][i], \
                                        'end time (UTC)': self.scheduled_table['end time (UTC)'][i], \
                                        'duration (minutes)': self.scheduled_table['duration (minutes)'][i], \
                                        'ra (h)': self.scheduled_table['ra (h)'][i], \
                                        'ra (m)': self.scheduled_table['ra (m)'][i], \
                                        'ra (s)': self.scheduled_table['ra (s)'][i], \
                                        'dec (d)': self.scheduled_table['dec (d)'][i], \
                                        'dec (m)': self.scheduled_table['dec (m)'][i], \
                                        'dec (s)': self.scheduled_table['dec (s)'][i], \
                                        'configuration': self.scheduled_table['configuration'][i]},
                                       ignore_index=True)

        end_scheduled_table = Table.from_pandas(end_scheduled_table)
        end_scheduled_table = unique(end_scheduled_table, keys='target')
        end_scheduled_table.sort(keys='start time (UTC)')
        idx_too_short_block = np.where((end_scheduled_table['duration (minutes)'] <= 3))

        if idx_too_short_block:
            for i in list(idx_too_short_block[0]):
                idx_for_target_2 = np.where((end_scheduled_table['target'] == end_scheduled_table['target'][i] + '_2'))[0]
                if idx_for_target_2.size != 0:
                    end_scheduled_table['target'][idx_for_target_2[0]] = end_scheduled_table['target'][i]

            list_too_short_reverse = list(idx_too_short_block[0])
            list_too_short_reverse.reverse() # important to reverse if you want to delete lines correctly
            for i in list_too_short_reverse:
                end_scheduled_table.remove_row(i)

        self.scheduled_table_sorted = end_scheduled_table
        # ***** Update target lists on server *****
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", path_spock + '/target_lists/target_list_special.txt',
                          'speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/spock_files/target_lists/'])
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", path_spock + '/target_lists/target_transit_follow_up.txt',
                          'speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/spock_files/target_lists/'])

        return self.scheduled_table_sorted

    def make_scheduled_table(self):
        Path= path_spock + '/DATABASE'
        try:
            os.path.exists(os.path.join(Path,self.telescope,'night_blocks_' + self.telescope + '_' +
                                        self.day_of_night.tt.datetime[0].strftime("%Y-%m-%d")+ '.txt'))
            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' Local path exists and is: ',os.path.join(Path,self.telescope,'night_blocks_' + self.telescope + '_'+
                                                    self.day_of_night.tt.datetime[0].strftime("%Y-%m-%d") + '.txt'))
        except TypeError:
            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' Local path does not exist yet ')
        except NameError:
            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' no input night_block for this day')
        except FileNotFoundError:
            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' no input night_block for this day')

        if not (self.scheduled_table is None):
            return self.scheduled_table
        else:
            try:
                self.scheduled_table = read_night_block(telescope=self.telescope,
                                                            day=self.day_of_night.tt.datetime.strftime("%Y-%m-%d"))
                return self.scheduled_table
            except TypeError:
                self.scheduled_table = read_night_block(telescope=self.telescope,
                                                            day=self.day_of_night.tt.datetime.strftime("%Y-%m-%d"))
                return self.scheduled_table

    def make_night_block(self):

        Path=  path_spock + '/night_blocks_propositions/'
        if not (self.scheduled_table_sorted is None):
            if isinstance(self.scheduled_table_sorted,pd.DataFrame):
                self.scheduled_table_sorted = self.scheduled_table_sorted.set_index('target')
                try:
                    self.scheduled_table_sorted.drop('TransitionBlock', inplace=True)
                except KeyError:
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' no transition block')
                panda_table = self.scheduled_table_sorted
                try:
                    panda_table.to_csv(os.path.join(Path,'night_blocks_' + self.telescope + '_' +  self.day_of_night.tt.datetime[0].strftime("%Y-%m-%d") + '.txt'),sep=' ')
                except TypeError:
                    panda_table.to_csv(os.path.join(Path, 'night_blocks_' + self.telescope + '_' +self.day_of_night.tt.datetime.strftime("%Y-%m-%d") + '.txt'),sep=' ')
            else:
                try:
                    self.scheduled_table_sorted.add_index('target')
                    index_to_delete = self.scheduled_table_sorted.loc['TransitionBlock'].index
                    self.scheduled_table_sorted.remove_row(index_to_delete)
                except KeyError:
                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' no transition block')
                panda_table = self.scheduled_table_sorted.to_pandas()
                try:
                    panda_table.to_csv(os.path.join(Path,'night_blocks_' + self.telescope + '_' +  self.day_of_night.tt.datetime[0].strftime("%Y-%m-%d") + '.txt'),sep=' ')
                except TypeError:
                    panda_table.to_csv(os.path.join(Path, 'night_blocks_' + self.telescope + '_' +self.day_of_night.tt.datetime.strftime("%Y-%m-%d") + '.txt'),sep=' ')

    def exposure_time(self, input_name,target_list,day=None):
        i = np.where((target_list['Sp_ID'] == input_name))[0]
        if day is None:
            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' Not using moon phase in ETC')
        # moon_phase = round(moon_illumination(Time(day.iso, out_subfmt='date')), 2)
        try:
            spectral_type = round(float(target_list['SpT'][i].data.data[0]))
        except NotImplementedError:
            try:
                spectral_type = round(float(target_list['SpT'][i]))
            except ValueError:
                spectral_type = target_list['SpT'][i].data.data[0]
        if not isinstance(spectral_type, str):
            if spectral_type <= 9:
                spt_type = 'M' + str(spectral_type)
            if spt_type == 'M3':
                    spt_type = 'M2'
            elif (spectral_type == 12) or (spectral_type == 15) or (
                    int(spectral_type) == 18):
                spt_type = 'M' + str(spectral_type - 10)
            elif spectral_type == 10:
                spt_type = 'M9'
            elif spectral_type == 11:
                spt_type = 'L2'
            elif spectral_type == 13:
                spt_type = 'L2'
            elif spectral_type == 14:
                spt_type = 'L5'
            elif spectral_type > 14:
                spt_type = 'L8'
        else:
            spt_type = spectral_type
        filt_ = str(target_list['Filter_spc'][i])
        if (filt_ == 'z\'') or (filt_ == 'r\'') or (filt_ == 'i\'') or (filt_ == 'g\''):
            filt_ = filt_.replace('\'', '')
        filters = ['I+z', 'z', 'i', 'r','g']
        filt_idx = 0
        filt_ = filters[filt_idx]
        if self.telescope == 'Saint-Ex':
            if target_list['J'][i] != 0.:
                a = (ETC.etc(mag_val=target_list['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.1,
                             moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=3.48, temp_ccd=-70,
                             observatory_altitude=2780))
            if (target_list['J'][i] == 0.) and (target_list['V'][i] != 0.) :
                a = (ETC.etc(mag_val=target_list['V'][i], mag_band='V', spt=spt_type, filt=filt_, airmass=1.1,
                             moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=3.48, temp_ccd=-70,
                             observatory_altitude=2780))
            else:
                sys.exit('ERROR: You must precise Vmag or Jmag for this target')
            texp = a.exp_time_calculator(ADUpeak=35000)[0]

        elif self.telescope == 'TS_La_Silla' or self.telescope == 'TN_Oukaimeden':
            if target_list['J'][i] != 0.:
                a = (ETC.etc(mag_val=target_list['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.1, \
                             moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
            if (target_list['J'][i] == 0.) and (target_list['V'][i] != 0.):
                a = (ETC.etc(mag_val=target_list['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.1, \
                             moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
            else:
                sys.exit('ERROR: You must precise Vmag or Jmag for this target')
            texp = a.exp_time_calculator(ADUpeak=50000)[0]
            print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Don\'t forget to  calculate exposure time for TRAPPIST observations!!')

        elif self.telescope == 'Artemis':
            if target_list['J'][i] != 0.:
                a = (ETC.etc(mag_val=target_list['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.1, \
                             moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
            if (target_list['J'][i] == 0.) and (target_list['V'][i] != 0.):
                a = (ETC.etc(mag_val=target_list['V'][i], mag_band='V', spt=spt_type, filt=filt_, airmass=1.1, \
                             moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
            else:
                sys.exit('ERROR: You must precise Vmag or Jmag for this target')
            texp = a.exp_time_calculator(ADUpeak=45000)[0]

        elif self.telescope == 'Io' or self.telescope == 'Europa' \
                or self.telescope == 'Ganymede' or self.telescope == 'Callisto':
            if target_list['J'][i] != 0.:
                a = (ETC.etc(mag_val=target_list['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.1, \
                             moonphase=0.5, irtf=0.8, num_tel=1, seeing=0.7, gain=1.1))
            if (target_list['J'][i] == 0.) and (target_list['V'][i] != 0.):
                a = (ETC.etc(mag_val=target_list['V'][i], mag_band='V', spt=spt_type, filt=filt_, airmass=1.1, \
                             moonphase=0.5, irtf=0.8, num_tel=1, seeing=0.7, gain=1.1))
            else:
                sys.exit('ERROR: You must precise Vmag or Jmag for this target')
            texp = a.exp_time_calculator(ADUpeak=45000)[0]

        while texp < 10:
            print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Change filter to avoid saturation!!')
            filt_idx += 1
            if (filt_idx >= 3):
                print(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' You have to defocus in we want to observe this target')
                texp = 10.0001
            filt_ = filters[filt_idx]
            if self.telescope == 'Saint-Ex':
                if target_list['J'][i] != 0.:
                    a = (ETC.etc(mag_val=target_list['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.1,
                                 moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=3.48, temp_ccd=-70,
                                 observatory_altitude=2780))
                if (target_list['J'][i] == 0.) and (target_list['V'][i] != 0.):
                    a = (ETC.etc(mag_val=target_list['V'][i], mag_band='V', spt=spt_type, filt=filt_, airmass=1.1,
                                 moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=3.48, temp_ccd=-70,
                                 observatory_altitude=2780))
                else:
                    sys.exit('ERROR: You must precise Vmag or Jmag for this target')
                texp = a.exp_time_calculator(ADUpeak=35000)[0]

            elif self.telescope == 'Artemis':
                target_list['Filter_spc'][i] = filt_
                if target_list['J'][i] != 0.:
                    a = (ETC.etc(mag_val=target_list['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.1,\
                                 moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
                if (target_list['J'][i] == 0.) and (target_list['V'][i] != 0.):
                    a = (ETC.etc(mag_val=target_list['V'][i], mag_band='V', spt=spt_type, filt=filt_, airmass=1.1,\
                                 moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
                else:
                    sys.exit('ERROR: You must precise Vmag or Jmag for this target')
                texp = a.exp_time_calculator(ADUpeak=45000)[0]

            elif self.telescope == 'TS_La_Silla' or self.telescope == 'TN_Oukaimeden':
                target_list['Filter_spc'][i] = filt_
                if target_list['J'][i] != 0.:
                    a = (ETC.etc(mag_val=target_list['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.1,\
                                 moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
                if (target_list['J'][i] == 0.) and (target_list['V'][i] != 0.):
                    a = (ETC.etc(mag_val=target_list['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.1,\
                                 moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
                else:
                    sys.exit('ERROR: You must precise Vmag or Jmag for this target')
                texp = a.exp_time_calculator(ADUpeak=50000)[0]
                print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Don\'t forget to  calculate exposure time for TRAPPIST observations!!')

            elif self.telescope == 'Io' or self.telescope == 'Europa' or self.telescope == 'Ganymede' or self.telescope == 'Callisto':
                target_list['Filter_spc'][i] = filt_
                if target_list['J'][i] != 0.:
                    a = (ETC.etc(mag_val=target_list['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.1,\
                                 moonphase=0.5, irtf=0.8, num_tel=1, seeing=0.7, gain=1.1))
                if (target_list['J'][i] == 0.) and (target_list['V'][i] != 0.):
                    a = (ETC.etc(mag_val=target_list['V'][i], mag_band='V', spt=spt_type, filt=filt_, airmass=1.1,\
                                 moonphase=0.5, irtf=0.8, num_tel=1, seeing=0.7, gain=1.1))
                else:
                    sys.exit('ERROR: You must precise Vmag or Jmag for this target')
                texp = a.exp_time_calculator(ADUpeak=45000)[0]

        return int(texp)

    def visibility_plot(self):
        day = self.day_of_night
        delta_midnight = Time(np.linspace(self.observatory.twilight_evening_nautical(day,which='next').jd, self.observatory.twilight_morning_nautical(day+1,which='nearest').jd, 100),format='jd')
        for i in range(len(self.scheduled_table_sorted)):
            idx = np.where((self.target_table_spc['Sp_ID'] ==  self.scheduled_table_sorted['target'][i]))[0]
            plot_airmass(self.targets[int(idx)], self.observatory, delta_midnight, style_sheet=dark_style_sheet)
            t = Time(self.scheduled_table_sorted['start time (UTC)'][i])
            plt.vlines(t.iso, 3, 1, color='r')
        plt.legend(shadow=True, loc=2)


def visibility_plot(day,observatory,night_block):
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day,which='next').jd, observatory.twilight_morning_nautical(day+1,which='nearest').jd, 100),format='jd')
    dec = str(int(float(night_block['dec (d)'][0]))) + 'd' + str(int(float(night_block['dec (m)'][0]))) + 'm' + str(int(float(night_block['dec (s)'][0]))) + 's'
    ra = str(int(float(night_block['ra (h)'][0]))) + 'h' + str(int(float(night_block['ra (m)'][0]))) + 'm' + str(int(float(night_block['ra (s)'][0]))) + 's'
    for i in range(len(night_block)):
        plot_airmass(SkyCoord(ra=ra,dec=dec,frame='icrs'), observatory, delta_midnight, style_sheet=dark_style_sheet)
        t = Time(night_block['start time (UTC)'][i])
        plt.vlines(t.iso, 3, 1,linestyle='-', color='r',alpha = 0.7)
        plt.legend(shadow=True, loc=2)

def save_schedule(input_file,nb_observatory,save,over_write,date_range,telescope):
    if input_file is None:
        telescope = telescope
        date_range = date_range
        try:
            day = date_range[0]
        except TypeError:
            day = date_range
    else:
        with open(input_file, "r") as f:
            Inputs = yaml.load(f, Loader=yaml.FullLoader)
            df = pd.DataFrame.from_dict(Inputs['observatories'])
            observatory = charge_observatories(df[nb_observatory]['name'])[0]
            date_range = Time(Inputs['day_of_night'])
            telescope = df[nb_observatory]['telescopes'][0]
            day = date_range[0]
    if save:
        source = path_spock + '/night_blocks_propositions/' + 'night_blocks_' + telescope + '_' +  day.tt.datetime.strftime("%Y-%m-%d") + '.txt'
        destination = path_spock + '/DATABASE/' + telescope + '/'
        destination_2 = path_spock + '/DATABASE/' + telescope + '/' + 'Archive_night_blocks/'
        if over_write:
            dest = shutil.copy(source, destination)
            dest2 = shutil.copy(source, destination_2)
            print(Fore.GREEN + 'INFO:  ' + Fore.BLACK + '\"' + source + '\"' + ' has been over-written to ' + '\"' +  destination + '\"' )
            print(Fore.GREEN + 'INFO:  ' + Fore.BLACK + '\"' + source + '\"' + ' has been over-written to ' + '\"' + destination_2 + '\"')
        if not over_write:
            try:
                dest = shutil.move(source, destination)
                dest2 = shutil.move(source, destination_2)
                print(Fore.GREEN + 'INFO:  ' + Fore.BLACK + '\"' +  source + '\"' +  ' has been copied to ' + '\"' + destination + '\"' )
                print(Fore.GREEN + 'INFO:  ' + Fore.BLACK + '\"' + source + '\"' + ' has been copied to ' + '\"' + destination_2 + '\"')
            except shutil.Error:
                print(Fore.GREEN + 'INFO:  ' + Fore.BLACK + '\"' + destination + 'night_blocks_' + telescope + '_' +  day.tt.datetime.strftime("%Y-%m-%d") + '.txt' + '\"' +  ' already exists')
    if not save:
        print(Fore.GREEN + 'INFO:  ' + Fore.BLACK + ' Those plans have not been saved')

def make_plans(day, nb_days, telescope):
    make_np(day, nb_days, telescope)

def upload_plans(day, nb_days, telescope):
    if telescope.find('Callisto') is not -1:
        SPOCKunp.upload_np_calli(day, nb_days)
    if telescope.find('Ganymede') is not -1:
        SPOCKunp.upload_np_gany(day, nb_days)
    if telescope.find('Io') is not -1:
        SPOCKunp.upload_np_io(day, nb_days)
    if telescope.find('Europa') is not -1:
        SPOCKunp.upload_np_euro(day, nb_days)
    if telescope.find('Artemis') is not -1:
        SPOCKunp.upload_np_artemis(day, nb_days)
    if telescope.find('TS_La_Silla') is not -1:
        SPOCKunp.upload_np_ts(day, nb_days)
    if telescope.find('TN_Oukaimeden') is not -1:
        SPOCKunp.upload_np_tn(day, nb_days)
    # ------------------- update archive date by date plans folder  ------------------

    path_gant_chart = os.path.join(path_spock + '/SPOCK_Figures/Preview_schedule.html')
    path_database_home = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/Preview_schedule.html')
    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' Path local \'Gant chart\' = ', path_gant_chart)
    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' Path database = \'Gant chart\' = ',  path_database_home)
    subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r", path_gant_chart, path_database_home ])
    path_gant_chart_masterfile = os.path.join('/Users/elsaducrot/spock_2/SPOCK_Figures/spock_stats_masterfile.csv')
    path_database_home_masterfile = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/spock_stats_masterfile.csv')
    subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r", path_gant_chart_masterfile, path_database_home_masterfile])


def read_night_block(telescope, day):
    day_fmt = Time(day, scale='utc', out_subfmt='date').tt.datetime.strftime("%Y-%m-%d")
    path_local = path_spock + '/DATABASE/' + telescope+'/Archive_night_blocks/night_blocks_'+\
                 telescope+'_'+day_fmt+'.txt'

    if os.path.exists(path_local):
        day_fmt = Time(day, scale='utc', out_subfmt='date').tt.datetime.strftime("%Y-%m-%d")
        scheduler_table = Table.read(
            path_spock + '/DATABASE/' + str(telescope) + '/Archive_night_blocks' + '/night_blocks_' + str(telescope) + '_' + str(day_fmt) + '.txt',
            format='ascii')
    else:
        nightb_url = "http://www.mrao.cam.ac.uk/SPECULOOS/"+telescope+'/schedule/Archive_night_blocks/night_blocks_'+\
                     telescope+'_'+day_fmt+'.txt'
        nightb = requests.get(nightb_url, auth=(user_portal, pwd_portal))

        if nightb.status_code == 404:
            sys.exit(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' No plans on the server for this date')
        else:
            open(path_local, 'wb').write(nightb.content)
            scheduler_table = pd.read_csv(path_local, delimiter=' ',
                                          skipinitialspace=True, error_bad_lines=False)


    return scheduler_table

def make_docx_schedule(observatory,telescope, date_range, name_operator):
    if not os.path.exists(path_spock + '/TRAPPIST_schedules_docx'):
        os.makedirs(path_spock + '/TRAPPIST_schedules_docx')
    df_speculoos = pd.read_csv(path_spock + '/target_lists/speculoos_target_list_v6.txt', delimiter=' ')
    df_follow_up = pd.read_csv(path_spock + '/target_lists/target_transit_follow_up.txt', delimiter=' ')
    df_special = pd.read_csv(path_spock + '/target_lists/target_list_special.txt', delimiter=' ')

    df_follow_up['nb_hours_surved'] = [0]*len(df_follow_up)
    df_follow_up['nb_hours_threshold'] = [0] * len(df_follow_up)
    df_special['nb_hours_surved'] = [0] * len(df_special)
    df_special['nb_hours_threshold'] = [0] * len(df_special)

    df_pandas = pd.DataFrame({'Sp_ID':pd.concat([df_speculoos['Sp_ID'],df_follow_up['Sp_ID'],df_special['Sp_ID']]),
        'RA': pd.concat([df_speculoos['RA'],df_follow_up['RA'],df_special['RA']]),
        'DEC': pd.concat([df_speculoos['DEC'],df_follow_up['DEC'],df_special['DEC']]),
        'J': pd.concat([df_speculoos['J'],df_follow_up['J'],df_special['J']]),
        'SpT': pd.concat([df_speculoos['SpT'], df_follow_up['SpT'], df_special['SpT']]),
        'nb_hours_surved': pd.concat([df_speculoos['nb_hours_surved'],df_follow_up['nb_hours_surved'],df_special['nb_hours_surved']]),
        'nb_hours_threshold': pd.concat([df_speculoos['nb_hours_threshold'],df_follow_up['nb_hours_threshold'],df_special['nb_hours_threshold']])})
    df_pandas = df_pandas.drop_duplicates()
    df = Table.from_pandas(df_pandas)
    nb_day_date_range = date_range_in_days(date_range)
    doc = Document()
    par = doc.add_paragraph()
    par_format = par.paragraph_format
    par_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    par_format.space_beFore = Pt(0)
    par_format.space_after = Pt(6)
    run = par.add_run(observatory.name)
    run.bold = True
    font = run.font
    font.size = Pt(16)
    font.color.rgb = RGBColor(0, 0, 0)
    par = doc.add_paragraph()
    par_format = par.paragraph_format
    par_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    par_format.space_beFore = Pt(0)
    par_format.space_after = Pt(12)
    run = par.add_run('Schedule from ' + Time(date_range[0].iso,out_subfmt='date').iso + ' to ' + Time(date_range[1].iso,out_subfmt='date').iso)
    run.bold = True
    font = run.font
    font.size = Pt(16)
    font.color.rgb = RGBColor(0, 0, 0)
    par = doc.add_paragraph()
    par_format = par.paragraph_format
    par_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    par_format.space_beFore = Pt(0)
    par_format.space_after = Pt(12)
    run = par.add_run(
        '(Total time = 0hr, technical loss = 0hr, weather loss = 0hr, Exotime = 0hr, cometime = 0hr,   chilean time = 0hr)')
    run.bold = True
    font = run.font
    font.size = Pt(12)
    font.color.rgb = RGBColor(255, 0, 0)
    par = doc.add_paragraph()
    par_format = par.paragraph_format
    par_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    par_format.space_beFore = Pt(0)
    par_format.space_after = Pt(20)
    run = par.add_run(name_operator)
    run.italic = True
    font = run.font
    font.size = Pt(12)
    font.color.rgb = RGBColor(0, 0, 0)
    par = doc.add_paragraph()
    par_format = par.paragraph_format
    par_format.space_beFore = Pt(16)
    par_format.space_after = Pt(0)

    for i in range(nb_day_date_range):

        date = date_range[0] + i
        table_schedule = read_night_block(telescope, date)
        sun_set = observatory.sun_set_time(date, which='next').iso
        sun_rise = observatory.sun_rise_time(date, which='next').iso
        moon_illum = int(round(moon_illumination(date) * 100, 0)) * u.percent
        civil_twilights = [Time(observatory.twilight_evening_civil(date, which='next')).iso,
                           Time(observatory.twilight_morning_civil(date + 1, which='nearest')).iso]
        nautic_twilights = [Time(observatory.twilight_evening_nautical(date, which='next')).iso,
                            Time(observatory.twilight_morning_nautical(date + 1, which='nearest')).iso]
        astro_twilights = [Time(observatory.twilight_evening_astronomical(date, which='next')).iso,
                           Time(observatory.twilight_morning_astronomical(date + 1, which='nearest')).iso]
        start_night = table_schedule['start time (UTC)'][0]
        end_night = np.array(table_schedule['end time (UTC)'])[-1]
        night_duration = round((Time(nautic_twilights[1]) - Time(nautic_twilights[0])).jd * 24, 3) * u.hour

        run = par.add_run('Night starting on the ' + Time(date, out_subfmt='date').value)
        run.bold = True
        run.underline = True
        font = run.font
        font.size = Pt(12)
        font.color.rgb = RGBColor(0, 0, 0)
        par = doc.add_paragraph()
        par_format = par.paragraph_format
        par_format.space_beFore = Pt(0)
        par_format.space_after = Pt(0)
        run = par.add_run('Moon illumination: ' + str(moon_illum))
        run.italic = True
        font = run.font
        font.size = Pt(12)
        font.color.rgb = RGBColor(0, 0, 0)
        par = doc.add_paragraph()
        par_format = par.paragraph_format
        par_format.space_beFore = Pt(0)
        par_format.space_after = Pt(0)
        run = par.add_run(
            'Sunset - Sunrise: ' + '{:02d}'.format(Time(sun_set, out_subfmt='date_hm').datetime.hour) + 'h' + '{:02d}'.format(Time(sun_set, out_subfmt='date_hm').datetime.minute) +\
            '  / ' + '{:02d}'.format(Time(sun_rise, out_subfmt='date_hm').datetime.hour) + 'h' + '{:02d}'.format(Time(sun_rise, out_subfmt='date_hm').datetime.minute))
        run.italic = True
        font = run.font
        font.size = Pt(12)
        font.color.rgb = RGBColor(0, 0, 0)
        par = doc.add_paragraph()
        par_format = par.paragraph_format
        par_format.space_beFore = Pt(0)
        par_format.space_after = Pt(0)
        run = par.add_run(
            'Civil/Naut./Astro. twilights: ' + \
            '{:02d}'.format(Time(civil_twilights[0], out_subfmt='date_hm').datetime.hour) + 'h' + '{:02d}'.format(Time(civil_twilights[0], out_subfmt='date_hm').datetime.minute)  +\
            '-' + '{:02d}'.format(Time(civil_twilights[1], out_subfmt='date_hm').datetime.hour) + 'h' + '{:02d}'.format(Time(civil_twilights[1], out_subfmt='date_hm').datetime.minute)  +\
            ' / ' + '{:02d}'.format(Time(nautic_twilights[0], out_subfmt='date_hm').datetime.hour) + 'h' + '{:02d}'.format(Time(nautic_twilights[0], out_subfmt='date_hm').datetime.minute) +\
            '-' + '{:02d}'.format(Time(nautic_twilights[1], out_subfmt='date_hm').datetime.hour) + 'h' + '{:02d}'.format(Time(nautic_twilights[1], out_subfmt='date_hm').datetime.minute) +\
            '  / ' + '{:02d}'.format(Time(astro_twilights[0], out_subfmt='date_hm').datetime.hour) + 'h' + '{:02d}'.format(Time(astro_twilights[0], out_subfmt='date_hm').datetime.minute) +\
            '-' + '{:02d}'.format(Time(astro_twilights[1], out_subfmt='date_hm').datetime.hour) + 'h' + '{:02d}'.format(Time(astro_twilights[1], out_subfmt='date_hm').datetime.minute))
        run.italic = True
        font = run.font
        font.size = Pt(12)
        font.color.rgb = RGBColor(0, 0, 0)
        par = doc.add_paragraph()
        par_format = par.paragraph_format
        par_format.space_beFore = Pt(0)
        par_format.space_after = Pt(0)
        run = par.add_run('Start-end of night (Naut. twil.): ' + '{:02d}'.format(Time(start_night).datetime.hour) + 'h' + '{:02d}'.format(Time(start_night).datetime.minute) +\
                          ' to ' + '{:02d}'.format(Time(end_night).datetime.hour) + 'h' + '{:02d}'.format(Time(end_night).datetime.minute))
        run.italic = True
        font = run.font
        font.size = Pt(12)
        font.color.rgb = RGBColor(0, 0, 0)
        par = doc.add_paragraph()
        par_format = par.paragraph_format
        par_format.space_beFore = Pt(0)
        par_format.space_after = Pt(3)
        run = par.add_run('Night duration (Naut. twil.): ' + str(night_duration))
        run.italic = True
        font = run.font
        font.size = Pt(12)
        font.color.rgb = RGBColor(0, 0, 0)

        for i in range(len(table_schedule)):
            trappist_planets = ['Trappist-1b','Trappist-1c','Trappist-1d','Trappist-1e',
                                'Trappist-1f','Trappist-1g','Trappist-1h']

            if any(table_schedule['target'][i] == p for p in trappist_planets):
                idx_target = np.where((df['Sp_ID'] == 'Sp2306-0502'))[0]
            else:
                idx_target = np.where((df['Sp_ID'] == table_schedule['target'][i]))[0]

            start_time_target = table_schedule['start time (UTC)'][i]
            end_time_target = table_schedule['end time (UTC)'][i]
            config = table_schedule['configuration'][i]
            try:
                coords = SkyCoord(ra=df['RA'][idx_target].data.data[0] * u.deg, dec=df['DEC'][idx_target].data.data[0] * u.deg)
            except IndexError:
                break
            dist_moon = '34'

            par = doc.add_paragraph()
            par_format = par.paragraph_format
            par_format.space_beFore = Pt(0)
            par_format.space_after = Pt(0)
            run = par.add_run(
                'From ' + '{:02d}'.format(Time(start_time_target, out_subfmt='date_hm').datetime.hour) + 'h' + '{:02d}'.format(Time(start_time_target, out_subfmt='date_hm').datetime.minute) +\
                ' to ' + '{:02d}'.format(Time(end_time_target, out_subfmt='date_hm').datetime.hour) + 'h' + '{:02d}'.format(Time(end_time_target, out_subfmt='date_hm').datetime.minute) +\
                ' : ' + str(table_schedule['target'][i]))
            run.bold = True
            font = run.font
            font.size = Pt(12)
            font.color.rgb = RGBColor(0, 0, 0)
            par = doc.add_paragraph()
            par_format = par.paragraph_format
            par_format.space_beFore = Pt(0)
            par_format.space_after = Pt(0)
            run = par.add_run('  Note: Prio_target                                         ')
            font = run.font
            font.size = Pt(10)
            font.color.rgb = RGBColor(0, 0, 0)
            par = doc.add_paragraph()
            par_format = par.paragraph_format
            par_format.space_beFore = Pt(0)
            par_format.space_after = Pt(0)
            run = par.add_run(
                '  SPECULOOS : ' + str(df['nb_hours_surved'][idx_target].data.data[0]*u.hour) + ' of obs over ' + str(
                    df['nb_hours_threshold'][idx_target].data.data[0]*u.hour))
            font = run.font
            font.size = Pt(10)
            font.color.rgb = RGBColor(0, 0, 0)
            par = doc.add_paragraph()
            par_format = par.paragraph_format
            par_format.space_beFore = Pt(0)
            par_format.space_after = Pt(0)
            run = par.add_run('Jmag= ' + str(df['J'][idx_target].data.data[0]) + ',  SpT= ' + str(
                df['SpT'][idx_target].data[0]))  # + ', Moon at ' + str(dist_moon))
            font = run.font
            font.size = Pt(12)
            font.color.rgb = RGBColor(0, 0, 0)
            par = doc.add_paragraph()
            par_format = par.paragraph_format
            par_format.space_beFore = Pt(0)
            par_format.space_after = Pt(3)
            run = par.add_run(' RA = ' + str('{:02d}'.format(int(coords.ra.hms[0]))) + " " + str('{:02d}'.format(int(coords.ra.hms[1]))) + " " + str('{:05.3f}'.format(round(coords.ra.hms[2], 3))) + \
                ', DEC = ' + str('{:02d}'.format(int(coords.dec.dms[0]))) + " " + str('{:02d}'.format(int(abs(coords.dec.dms[1])))) + " " + str('{:05.3f}'.format(round(abs(coords.dec.dms[2]), 3))) + ', ' + str(config[2:-2]).replace('\'',' '))
            font = run.font
            font.size = Pt(12)
            font.color.rgb = RGBColor(0, 0, 0)
            par = doc.add_paragraph()
            par_format = par.paragraph_format
            par_format.space_beFore = Pt(16)
            par_format.space_after = Pt(0)

    font = run.font
    font.size = Pt(12)
    font.color.rgb = RGBColor(0, 0, 0)
    if telescope == 'TS_La_Silla':
        doc.save(path_spock + '/TRAPPIST_schedules_docx/TS_' + Time(date_range[0],out_subfmt='date').value.replace('-','') + '_to_' +\
                 Time(date_range[1],out_subfmt='date').value .replace('-','')  +'.docx')
    if telescope == 'TN_Oukaimeden':
        doc.save(path_spock + '/TRAPPIST_schedules_docx/TN_' + Time(date_range[0],out_subfmt='date').value.replace('-','') + '_to_' +\
                 Time(date_range[1],out_subfmt='date').value .replace('-','')  +'.docx')

def date_range_in_days(date_range):
    date_format = "%Y-%m-%d %H:%M:%S.%f"
    date_start = datetime.strptime(date_range[0].value, date_format)
    date_end = datetime.strptime(date_range[1].value, date_format)
    date_range_in_days = (date_end - date_start).days
    return date_range_in_days


def prediction(obs_name,name,ra,dec,timing,period,duration, start_date,ntr):
    start_date = Time(start_date)

    constraints = [AltitudeConstraint(min=25 * u.deg), AtNightConstraint()]

    target_transit = EclipsingSystem(primary_eclipse_time=Time(timing, format='jd'),
                                     orbital_period=period*u.day, duration=duration*u.day,
                                     name=name)

    mid_transit_timing = Time(target_transit.next_primary_eclipse_time(start_date, n_eclipses=ntr)).iso
    mid_transit_timing_jd = Time(target_transit.next_primary_eclipse_time(start_date, n_eclipses=ntr)).jd

    ing_egr = target_transit.next_primary_ingress_egress_time(start_date, n_eclipses=ntr)

    target = [FixedTarget(coord=SkyCoord(ra=ra * u.degree, dec=dec * u.degree),
                    name=name)]
    #observatory = charge_observatories(obs_name)[0]

    observable_SSO = is_event_observable(constraints,charge_observatories('SSO')[0], target,
                                         times_ingress_egress=ing_egr)
    observable_SNO = is_event_observable(constraints, charge_observatories('SNO')[0], target,
                                         times_ingress_egress=ing_egr)
    observable_SaintEx = is_event_observable(constraints, charge_observatories('Saint-Ex')[0], target,
                                         times_ingress_egress=ing_egr)
    observable_TN_Oukaimeden = is_event_observable(constraints, charge_observatories('TN_Oukaimeden')[0], target,
                                         times_ingress_egress=ing_egr)
    observable_TS_La_Silla = is_event_observable(constraints, charge_observatories('TS_La_Silla')[0], target,
                                         times_ingress_egress=ing_egr)


    df = pd.DataFrame(data=ing_egr.sort().iso, index=[name]*len(ing_egr),columns=["Ingress", "Egress"])

    df['mid transit'] = mid_transit_timing

    df['mid transit JD'] = np.round(mid_transit_timing_jd,3)

    df['Observable SSO'] = observable_SSO[0]

    df['Observable SNO'] = observable_SNO[0]

    df['Observable Saint-Ex'] = observable_SaintEx[0]

    df['Observable TN Oukaimeden'] = observable_TN_Oukaimeden[0]

    df['Observable TS La Silla'] = observable_TS_La_Silla[0]

    return df



    #
    # for i in range(len(df['Sp_ID'])):
    #     blocks = []
    #     epoch = Time(df['T0'][i], format='jd')
    #     period = df['P'][i] * u.day
    #     duration = df['W'][i] * u.day
    #     oot_time = duration.value * 1.5 * u.day
    #     T0_err_transit = df['T0_err'][i]
    #     P_err_transit = df['P_err'][i]
    #     W_err_transit = df['W_err'][i]
    #     target_transit = EclipsingSystem(primary_eclipse_time = epoch, orbital_period=period, duration=duration, name=df['Sp_ID'][i])
    #     print(Fore.GREEN + Fore.GREEN + 'INFO: ' + Fore.BLACK + ' ' + Fore.BLACK + str(df['Sp_ID'][i]) + ' next transit: ', Time(target_transit.next_primary_eclipse_time(self.day_of_night, n_eclipses=1)).iso)
    #     timing_to_obs_jd = Time(target_transit.next_primary_eclipse_time(self.day_of_night, n_eclipses=1)).jd
    #     target = target_list_good_coord_format(follow_up_list)
    #     n_transits = 1
    #     try:
    #         ing_egr = target_transit.next_primary_ingress_egress_time(self.day_of_night,n_eclipses=n_transits)
    #     except ValueError:
    #         print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' No transit of ', df['Sp_ID'][i], ' on the period chosen')
    #
    #     observable = is_event_observable(constraints,self.observatory, target, times_ingress_egress=ing_egr)
    #     if np.any(observable):
    #         err_T0_neg = timing_to_obs_jd[0] - (np.round((timing_to_obs_jd[0] - epoch.jd) / period.value,1) *
    #                                             (period.value - P_err_transit) + (epoch.jd - T0_err_transit ))
    #         err_T0_pos = (np.round((timing_to_obs_jd[0] - epoch.jd) / period.value,1) *
    #                       (period.value + P_err_transit) + (epoch.jd + T0_err_transit )) - timing_to_obs_jd[0]
    #         start_transit = Time(ing_egr[0][0].value - err_T0_neg - oot_time.value - W_err_transit,format='jd') #- T0_err_transit - W_err_transit  - oot_time.value/24 - (timing_to_obs_jd[0] - epoch.jd) / period.value * P_err_transit,format='jd')
    #         #print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit', start_transit.iso)
    #         end_transit = Time(ing_egr[0][1].value + err_T0_pos + oot_time.value + W_err_transit,format='jd') #+ T0_err_transit + W_err_transit + oot_time.value/24 + (timing_to_obs_jd[0] - epoch.jd) / period.value * P_err_transit ,format='jd')
    #         #print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit', end_transit.iso)
    #         dur_obs_transit_target = (end_transit - start_transit).value * 1. * u.day
    #         if (end_transit > Time(self.observatory.twilight_morning_nautical(self.day_of_night + 1, which='nearest'))) \
    #                 or (start_transit < Time(self.observatory.twilight_evening_nautical(self.day_of_night, which='next'))):
    #             if (Time(ing_egr[0][1]) < Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')))\
    #                     and (Time(ing_egr[0][0]) > Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next'))):
    #                 constraints_transit_target = [AltitudeConstraint(min=25 * u.deg),TimeConstraint(start_transit, end_transit),
    #                                               MoonSeparationConstraint(min=25 * u.deg),AtNightConstraint(max_solar_altitude=-12*u.deg)]
    #                 print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ', str(df['Sp_ID'][i]), ' : ', start_transit.iso)
    #                 print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit of ', str(df['Sp_ID'][i]), ' : ', end_transit.iso)
    #                 print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Time out of transit not optimal.')
    #                 idx_first_target = int(np.where((df['Sp_ID'] == df['Sp_ID'][i]))[0])
    #                 if df['texp_spc'][idx_first_target] == 0:
    #                     df['texp_spc'][idx_first_target] = self.exposure_time(
    #                         input_name=df['Sp_ID'][idx_first_target],target_list=self.target_table_spc_follow_up)
    #                 a = ObservingBlock(target[idx_first_target], dur_obs_transit_target, -1,
    #                                    constraints=constraints_transit_target,
    #                                    configuration={'filt=' + str(df['Filter_spc'][idx_first_target]),
    #                                                   'texp=' + str(df['texp_spc'][idx_first_target])})
    #                 blocks.append(a)
    #                 transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
    #                 seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night + 1)
    #                 sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_transit_target,
    #                                                           observer=self.observatory,
    #                                                           transitioner=transitioner)
    #                 sequen_scheduler_SS1(blocks, seq_schedule_SS1)
    #
    #             else:
    #                 print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ', str(df['Sp_ID'][i]), ' : ', Time(ing_egr[0][0]).iso)
    #                 print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit of ', str(df['Sp_ID'][i]), ' : ', Time(ing_egr[0][1]).iso)
    #                 print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Transit not full.')
    #                 constraints_transit_target = [AltitudeConstraint(min=25 * u.deg),
    #                                               MoonSeparationConstraint(min=25 * u.deg),
    #                                               TimeConstraint(Time(ing_egr[0][0]),Time(ing_egr[0][1]))]
    #                 #continue
    #                 idx_first_target = int(np.where((df['Sp_ID'] == df['Sp_ID'][i]))[0])
    #                 if df['texp_spc'][idx_first_target] == 0:
    #                     df['texp_spc'][idx_first_target] = self.exposure_time(
    #                         input_name=df['Sp_ID'][idx_first_target],target_list=self.target_table_spc_follow_up)
    #                 a = ObservingBlock(target[idx_first_target], dur_obs_transit_target, -1,
    #                                    constraints=constraints_transit_target,
    #                                    configuration={'filt=' + str(df['Filter_spc'][idx_first_target]),
    #                                                   'texp=' + str(df['texp_spc'][idx_first_target])})
    #                 blocks.append(a)
    #                 transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
    #                 seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night + 1)
    #                 sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_transit_target,
    #                                                           observer=self.observatory,
    #                                                           transitioner=transitioner)
    #                 sequen_scheduler_SS1(blocks, seq_schedule_SS1)
    #
    #         else:
    #             constraints_transit_target = [AltitudeConstraint(min=25 * u.deg),
    #                                           MoonSeparationConstraint(min=25 * u.deg),
    #                                           TimeConstraint(start_transit, end_transit)]
    #             idx_first_target = int(np.where((df['Sp_ID'] == df['Sp_ID'][i]))[0])
    #             if df['texp_spc'][idx_first_target] == 0:
    #                 df['texp_spc'][idx_first_target] = self.exposure_time(input_name=df['Sp_ID'][idx_first_target],
    #                                                                       target_list=self.target_table_spc_follow_up)
    #             a = ObservingBlock(target[idx_first_target], dur_obs_transit_target, -1,
    #                                constraints=constraints_transit_target,
    #                                configuration={'filt=' + str(df['Filter_spc'][idx_first_target]),
    #                                               'texp=' + str(df['texp_spc'][idx_first_target])})
    #             blocks.append(a)
    #             transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
    #             seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night+1)
    #             sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_transit_target, observer=self.observatory,
    #                                                       transitioner=transitioner)
    #             sequen_scheduler_SS1(blocks, seq_schedule_SS1)
    #             print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ',str(df['Sp_ID'][i]), ' : ',Time(ing_egr[0][0]).iso)
    #             print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit of ',str(df['Sp_ID'][i]), ' : ',Time(ing_egr[0][1]).iso)
    #             print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' Transit is expected to be full.')
    #             if len(seq_schedule_SS1.to_table()['target']) == 0:
    #                 print(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' The moon is too closed for the transit to be observed')
    #
    #         if len(seq_schedule_SS1.to_table()['target']) != 0:
    #             self.SS1_night_blocks = seq_schedule_SS1.to_table()  # Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
    #             return self.SS1_night_blocks
    #     else:
    #         print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' no transit of ', df['Sp_ID'][i],' this day')
