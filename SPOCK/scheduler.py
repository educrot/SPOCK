from astroplan import FixedTarget, Observer, AltitudeConstraint, MoonSeparationConstraint,\
    AtNightConstraint, is_observable, observability_table, moon_illumination, TimeConstraint
from astroplan.periodic import EclipsingSystem
from astroplan.constraints import is_event_observable
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.table import Table
from astropy.time import Time
from colorama import Fore
from datetime import date, datetime, timedelta
from eScheduler.spe_schedule import SPECULOOSScheduler, Schedule, ObservingBlock,Transitioner
import functools
import numpy as np
import os
import pandas as pd
from SPOCK import path_spock
import SPOCK.short_term_scheduler as SPOCKST
import SPOCK.ETC as ETC
import sys
from tqdm.auto import tqdm


def first_elem_list(x):
    """ first element of list

    Parameters
    ----------
    x : list

    Returns
    -------
    float
        first element of list

    """
    a = x[0]
    return a


def last_elem_list(x):
    """ last element of list

    Parameters
    ----------
    x : list

    Returns
    -------
    float
        last element

    """
    a = x[-1]
    return a


def target_list_good_coord_format(path_target_list):

    """ Give target corrdinates in ICRS format (used for astropy.coordinates SkyCoord function)

    Parameters
    ----------
    path_target_list: str
        path on your computer toward the target list, by default take the one on the Cambridge server

    Returns
    -------
    targets: astropy.FixedTarget
        targets list with the following format : [<FixedTarget "Sp0002+0115" at SkyCoord (ICRS):
        (ra, dec) in deg (0.52591667, 1.26003889)>,


    """
    df = pd.read_csv(path_target_list, delimiter=' ')
    target_table_spc = Table.from_pandas(df)
    targets = [FixedTarget(coord=SkyCoord(ra=target_table_spc['RA'][i]* u.degree,
                                          dec=target_table_spc['DEC'][i] * u.degree),
                           name=target_table_spc['Sp_ID'][i]) for i in range(len(target_table_spc['RA']))]
    return targets


def charge_observatories(name):
    """

    Parameters
    ----------
    name : str
        name of the observatory (ex: 'SSO')

    Returns
    -------
    astroplan.observer
        all info on observatory loaded

    """
    observatories = []
    # Observatories
    if 'SSO' in str(name):
        location = EarthLocation.from_geodetic(-70.40300000000002*u.deg, -24.625199999999996*u.deg,
                                               2635.0000000009704*u.m)
        observatories.append(Observer(location=location, name="SSO", timezone="UTC"))

    if 'SNO' in str(name):
        location_sno = EarthLocation.from_geodetic(-16.50583131*u.deg, 28.2999988*u.deg, 2390*u.m)
        observatories.append(Observer(location=location_sno, name="SNO", timezone="UTC"))

    if 'Saint-Ex' in str(name):
        location_saintex = EarthLocation.from_geodetic(-115.48694444444445*u.deg, 31.029166666666665*u.deg,
                                                       2829.9999999997976*u.m)
        observatories.append(Observer(location=location_saintex, name="Saint-Ex", timezone="UTC"))

    if 'TS_La_Silla' in str(name):
        location_tslasilla = EarthLocation.from_geodetic(-70.73000000000002*u.deg, -29.25666666666666*u.deg,
                                                         2346.9999999988418*u.m)
        observatories.append(Observer(location=location_tslasilla, name="TS_La_Silla", timezone="UTC"))

    if 'TN_Oukaimeden' in str(name):
        location_tnouka = EarthLocation.from_geodetic(-7.862263*u.deg, 31.20516*u.deg, 2751*u.m)
        observatories.append(Observer(location=location_tnouka, name="TN_Oukaimeden", timezone="UTC"))

    if 'Munich' in str(name):
        location_munich = EarthLocation.from_geodetic(48.2*u.deg, -11.6*u.deg, 600*u.m)
        observatories.append(Observer(location=location_munich, name="Munich", timezone="UTC"))

    return observatories


def reverse_observability(observatory, targets, constraints, time_ranges):
    """
    Reverse observability table, rows become columns
    Parameters
    ----------
    observatory : str
        observatory chosen
    targets : list of astropy.FixedTarget
        target list on the FixedTarget() format from astroplan
    constraints : astroplan.constraints
        general constraints for a target to be scheduled
    time_ranges : list
        List of astropy.Time range with start and end times

    Returns
    -------
    reverse_df1 : astropy.table
        observability table with no NaN value (0 instead) inverse with targets as columns
        and elements of time_ranges (months) as rows
    """
    start_fmt = Time(time_ranges[0][0].iso, out_subfmt ='date').iso
    end_fmt = Time(time_ranges[len(time_ranges)-1][1].iso, out_subfmt='date').iso

    if os.path.exists(path_spock + '/SPOCK_files/reverse_Obs_' + str(observatory.name) + '_' + start_fmt + '_' +
                      end_fmt + '_' + str(len(targets)) + '.csv'):
        name_file = path_spock + '/SPOCK_files/reverse_Obs_' + str(observatory.name) + '_' + start_fmt + '_' + \
                    end_fmt + '_' + str(len(targets)) + '.csv'
        reverse_df1 = pd.read_csv(name_file, delimiter=',')
        return reverse_df1
    else:
        tables_observability = list(map((lambda x: observability(x, time_ranges[x], observatory, targets, constraints)),
                                        range(0, len(time_ranges))))
        df = list(map((lambda x: pd.DataFrame(tables_observability[x].to_pandas())), range(0, len(time_ranges))))
        a = functools.reduce(functools.partial(pd.merge, how='outer', on='target name'), df)
        df = a.replace(to_replace=float('NaN'), value=0.0)
        df1 = df.set_index('target name')
        reverse_df1 = df1.T
        reverse_df1.to_csv(path_spock + '/SPOCK_files/reverse_Obs_' + str(observatory.name) + '_' + start_fmt + '_' +
                           end_fmt + '_' + str(len(targets)) + '.csv', sep=',')
        return reverse_df1


def observability(j, time_range, observatory, targets, constraints):
    """ Give a table with the observability score for each target of targets
    regarding the constraints given and for all ranges of time_range

    Parameters
    ----------
        j : list
            [start end], element of time range (that is to say, month of the year)
        time_range : list
            of astropy.Time range with start and end times
        observatory : astroplan.observer
            observatory chosen
        targets : astropy
            target list on the FixedTarget() format from astroplan
        constraints : astroplan.constraint
            general constraints for a target to be shceduled

    Returns
    -------
        Observability table: astropy.table.Table
             12 columns (target name and each element of time_range, e.g months),
        rows= nb of targets
    """
    targets_observable = []
    observable = np.argwhere(is_observable(constraints, observatory, targets, time_range=time_range))

    # WARNING: Need to be replace by a np.where, but I can't find how
    [targets_observable.append(targets[int(obs)]) for obs in observable]

    table_observability = observability_table(constraints, observatory, targets_observable, time_range=time_range,
                                              time_grid_resolution=0.5*u.hour)  # compute a second time is observable
    table_observability.remove_column('always observable')
    table_observability.remove_column('ever observable')
    table_observability.rename_column('fraction of time observable', 'Month' + str(j))

    return table_observability


def month_option(target_name, reverse_df1):
    """ create a list of the best month for oservation for each target

    Parameters
    ----------
    target_name : str
        name of the target
    reverse_df1 : astropy.table
        observability table with no NaN value (0 instead) inversed with targets as columns
        and elements of time_ranges (months) as rows

    Returns
    -------
    month : list
        a list with the best month to observe the target
    month_2nd_option : list
        same but for the second best month
    months_3rd_option : list
        same but for the third best month
    months_4th_option : list
        same but for the fourth best month
    months_5th_option : list
        same but for the fiveth best month


    Remarks
    -------
        the 2nd, 3rd etc choices are here in case the target list is not long enough to give
        solutions for all the telescopes each months, allows to avoid blancks in observations
    """

    try:
        months = reverse_df1[str(target_name)].idxmax()
        months_2nd_option = reverse_df1[str(target_name)].nlargest(2, keep='first').index[1]
        months_3rd_option = reverse_df1[str(target_name)].nlargest(3, keep='first').index[2]
        months_4th_option = reverse_df1[str(target_name)].nlargest(4, keep='first').index[3]
        months_5th_option= reverse_df1[str(target_name)].nlargest(5, keep='first').index[4]

    except KeyError:
        months = 0
        months_2nd_option = 0
        months_3rd_option = 0
        months_4th_option = 0
        months_5th_option = 0
    return [months, months_2nd_option, months_3rd_option, months_4th_option, months_5th_option]


class schedules:
    """
    Class to Make schedules for the target list, observatory, date_range and startegy indicated

    """

    def __init__(self):
        self.Altitude_constraint = 25
        self.Moon_constraint = 30
        self.observatory_name = 'TS_La_Silla'
        self.telescope = 'TS_La_Silla'
        self.observatory = None
        self.date_range = None
        self.target_list = path_spock + '/target_lists/target_transit_follow_up_mathilde.txt'
        self.constraints = [AtNightConstraint()]
        self.reverse_df1 = None
        self.priority_table_ranked = None
        self.priority_table = None
        self.index_prio = None
        self.priority_table_ranked_by_day = []
        self.idx_first_target = None
        self.first_target = None
        self.moon_and_visibility_constraint_table = None
        self.observability_table_day = None
        self.ntr = 10
        self.day_of_night = None
        self.night_block = None
        self.next_transit_window = None
        self.list_next_transit_window = []

    @property
    def target_table_spc(self):
        df = pd.read_csv(self.target_list, delimiter=' ')
        return Table.from_pandas(df)

    @property
    def targets(self):
        return target_list_good_coord_format(self.target_list)

    @property
    def date_ranges_day_by_day(self):
        """ date range day by day

        Returns
        -------
        list
            list of date ranges

        """
        date_format = "%Y-%m-%d %H:%M:%S.%f"
        d2 = datetime.strptime(self.date_range[1].value, date_format)
        i = 0
        t = datetime.strptime(self.date_range[0].value, date_format)
        t_init = t
        date_ranges_day_by_day = []
        while t < d2:
            d = timedelta(days=1)
            t = t_init + d * i
            date_ranges_day_by_day.append(Time(t))
            i += 1
        return date_ranges_day_by_day

    @property
    def date_range_in_days(self):
        """ number of days in date range

        Returns
        -------
        int
            number of day between date start and date end

        """
        date_format = "%Y-%m-%d %H:%M:%S.%f"
        date_start = datetime.strptime(self.date_range[0].value, date_format)
        date_end = datetime.strptime(self.date_range[1].value, date_format)
        date_range_in_days = (date_end - date_start).days
        return date_range_in_days

    @property
    def months_obs(self):
        """ month of obs

        Returns
        -------
        int
            month number (between 0 and 11)

        """
        date_format = "%Y-%m-%d %H:%M:%S.%f"
        for i, t in enumerate(self.time_ranges):
            if (datetime.strptime(t[0].value, date_format) <=
                datetime.strptime(self.date_range[0].value, date_format) <= datetime.strptime(t[1].value, date_format))\
                    and (datetime.strptime(t[0].value, date_format) <=
                         datetime.strptime(self.date_range[1].value, date_format) <=
                         datetime.strptime(t[1].value, date_format)):
                return i
            if (datetime.strptime(t[0].value, date_format) <=
                datetime.strptime(self.date_range[0].value, date_format) <= datetime.strptime(t[1].value, date_format))\
                    and (datetime.strptime(t[1].value, date_format) <=
                         datetime.strptime(self.date_range[1].value, date_format)):
                if i < (len(self.time_ranges) - 1):
                    return i+1
                if i == (len(self.time_ranges) - 1):
                    return i

    @property
    def time_ranges(self):
        year = self.date_range[0].tt.datetime.strftime("%Y")
        time_ranges = [Time([year + '-01-01 12:00:00', year + '-01-31 12:00:00']),
                       Time([year + '-02-01 12:00:00', year + '-02-28 12:00:00']),
                       Time([year + '-03-01 15:00:00', year + '-03-31 15:00:00']),
                       Time([year + '-04-01 15:00:00', year + '-04-30 15:00:00']),
                       Time([year + '-05-01 15:00:00', year + '-05-31 15:00:00']),
                       Time([year + '-06-01 15:00:00', year + '-06-30 15:00:00']),
                       Time([year + '-07-01 12:00:00', year + '-07-31 12:00:00']),
                       Time([year + '-08-01 12:00:00', year + '-08-31 12:00:00']),
                       Time([year + '-09-01 12:00:00', year + '-09-30 12:00:00']),
                       Time([year + '-10-01 12:00:00', year + '-10-31 12:00:00']),
                       Time([year + '-11-01 12:00:00', year + '-11-30 12:00:00']),
                       Time([year + '-12-01 12:00:00', year + '-12-31 12:00:00'])]
        return time_ranges

    def load_parameters(self, date_range=None):
        self.observatory = charge_observatories(self.observatory_name)[0]
        if date_range is not None:
            self.date_range = Time(date_range)
        if self.Altitude_constraint:
            self.constraints.append(AltitudeConstraint(min=float(self.Altitude_constraint)*u.deg))
        if self.Moon_constraint:
            self.constraints.append(MoonSeparationConstraint(min=float(self.Moon_constraint)*u.deg))

    def observability_selection(self, day):
        day_fmt = Time(day.iso, out_subfmt='date').iso
        if os.path.exists(path_spock + '/SPOCK_files/sc_Ranking_months_' + str(self.observatory.name) + '_' +
                          str(day_fmt) + '_ndays_' + str(self.date_range_in_days) + '_' + str(
                len(self.targets)) + '.csv'):
            name_file = path_spock + '/SPOCK_files/sc_Ranking_months_' + str(self.observatory.name) + '_' +\
                        str(day_fmt) + '_ndays_' + str(self.date_range_in_days) + '_' + str(
                len(self.targets)) + '.csv'
            dataframe_ranking_months = pd.read_csv(name_file, delimiter=',')
            self.priority_table = Table.from_pandas(dataframe_ranking_months)
            # self.index_prio = np.argsort(self.priority_table['priority'])
            # self.priority_table_ranked = self.priority_table[self.index_prio]

        else:
            start_night_start = self.observatory.twilight_evening_nautical(day, which='nearest')  # * u.hour
            delta_midnight_start = np.linspace(0, self.observatory.twilight_morning_nautical(day, which='next').jd -
                                               self.observatory.twilight_evening_nautical(day, which='nearest').jd,
                                               100) * u.day  # Delta at the first day of schedule
            frame_start = AltAz(obstime=start_night_start + delta_midnight_start,
                                location=self.observatory.location)

            start_night_end = self.observatory.twilight_evening_nautical(day + self.date_range_in_days,
                                                                         which='nearest')  # * u.hour
            delta_midnight_end = np.linspace(0, self.observatory.twilight_morning_nautical(day +
                                                                                           self.date_range_in_days,
                                                                                           which='next').jd
                                             - self.observatory.twilight_evening_nautical(day +
                                                                                          self.date_range_in_days,
                                                                                          which='nearest').jd,
                                             100) * u.day  # Delta at the first day of schedule
            frame_end = AltAz(obstime=start_night_end + delta_midnight_end, location=self.observatory.location)

            target_alt = [[target.coord.transform_to(frame_start).alt,
                           target.coord.transform_to(frame_end).alt] for target in self.targets]  # This line takes time
            target_alt_start = np.asarray(target_alt)[:, 0, :]
            target_alt_end = np.asarray(target_alt)[:, 1, :]
            max_target_alt = list(map(max, target_alt_start))
            alt_set_start = list(map(first_elem_list, target_alt_start[:]))
            alt_rise_start = list(map(last_elem_list, target_alt_start[:]))
            alt_set_end = list(map(first_elem_list, target_alt_end[:]))
            alt_rise_end = list(map(last_elem_list, target_alt_end[:]))
            priority_value = [-0.5] * len(self.targets)  # For Mathilde: adapt this line
            set_or_rise = ['None'] * len(self.targets)
            df = pd.DataFrame({'priority': priority_value, 'set or rise': set_or_rise, 'alt set start': alt_set_start,
                               'alt rise start': alt_rise_start, 'alt set end': alt_set_end,
                               'alt rise end': alt_rise_end, 'max_alt': max_target_alt,
                               'Sp_ID': self.target_table_spc['Sp_ID']})
            self.priority_table = Table.from_pandas(df)
            month_opt = Table([[], [], [], [], []], names=['months', 'months_2nd', 'months_3rd',
                                                           'months_4th', 'months_5th'])
            [month_opt.add_row(month_option(target, self.reverse_df1)) for target in self.target_table_spc['Sp_ID']]
            # This line takes times
            idx_1rst_opt_monthobs = np.where((month_opt['months'] == self.months_obs))
            idx_2nd_opt_monthobs = np.where((month_opt['months_2nd'] == self.months_obs))
            idx_3rd_opt_monthobs = np.where((month_opt['months_3rd'] == self.months_obs))
            idx_4th_opt_monthobs = np.where((month_opt['months_4th'] == self.months_obs))
            idx_5th_opt_monthobs = np.where((month_opt['months_5th'] == self.months_obs))
            self.priority_table['priority'][idx_1rst_opt_monthobs] = \
                (self.priority_table['max_alt'][idx_1rst_opt_monthobs] - 30) * 10**4
            self.priority_table['priority'][idx_2nd_opt_monthobs] = \
                (self.priority_table['max_alt'][idx_2nd_opt_monthobs] - 30) * 10**3
            self.priority_table['priority'][idx_3rd_opt_monthobs] = \
                (self.priority_table['max_alt'][idx_3rd_opt_monthobs] - 30) * 10**2
            self.priority_table['priority'][idx_4th_opt_monthobs] = \
                (self.priority_table['max_alt'][idx_4th_opt_monthobs] - 30) * 10**1
            self.priority_table['priority'][idx_5th_opt_monthobs] = \
                (self.priority_table['max_alt'][idx_5th_opt_monthobs] - 30) * 10**0

            set_targets_index = (self.priority_table['alt set start'] > self.Altitude_constraint) & \
                                (self.priority_table['alt set end'] > self.Altitude_constraint)
            self.priority_table['set or rise'][set_targets_index] = 'set'

            rise_targets_index = (self.priority_table['alt rise start'] > self.Altitude_constraint) & \
                                 (self.priority_table['alt rise end'] > self.Altitude_constraint)
            self.priority_table['set or rise'][rise_targets_index] = 'rise'

            both_targets_index = (self.priority_table['alt rise start'] > self.Altitude_constraint) & \
                                 (self.priority_table['alt set start'] > self.Altitude_constraint)
            self.priority_table['set or rise'][both_targets_index] = 'both'

            not_observable_targets_index = (self.priority_table['max_alt'] < self.Altitude_constraint)
            self.priority_table['set or rise'][not_observable_targets_index] = 'not_observable'
            self.priority_table['priority'][not_observable_targets_index] = -0.5

            # self.index_prio = np.argsort(self.priority_table['priority'])
            # self.priority_table_ranked = self.priority_table[self.index_prio]

            dataframe_priority = self.priority_table.to_pandas()
            dataframe_priority.to_csv(path_spock + '/SPOCK_files/sc_Ranking_months_' + str(self.observatory.name) +
                                      '_' + str(day_fmt) + '_ndays_' + str(self.date_range_in_days) + '_' +
                                      str(len(self.targets)) + '.csv', sep=',', index=False)

    def exposure_time(self, day, i, telescope=None):
        """ calculation of the exposure time for a given target

        Parameters
        ----------
        day : date
            format 'yyyy-mm-dd'
        i : int
            index of target in target_list
        telescope: str
            name of the telescope

        Returns
        -------
        float
            exposure time

        """
        if telescope is not None:
            self.telescope = telescope
        if day is None:
            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' Not using moon phase in ETC')
        # moon_phase = round(moon_illumination(Time(day.iso, out_subfmt='date')), 2)
        if round(float(self.target_table_spc['SpT'][i])) <= 9:
            spt_type = 'M' + str(round(float(self.target_table_spc['SpT'][i])))
            if spt_type == 'M3':
                spt_type = 'M2'
        if round(float(self.target_table_spc['SpT'][i])) <= 2:
            spt_type = 'M2'
        elif (round(float(self.target_table_spc['SpT'][i])) == 12) or (round(float(self.target_table_spc['SpT'][i])) == 15)\
                or (int(float(self.target_table_spc['SpT'][i])) == 18):
            spt_type = 'M' + str(round(self.target_table_spc['SpT'][i]) - 10)
        elif round(float(self.target_table_spc['SpT'][i])) == 10:
            spt_type = 'M9'
        elif round(float(self.target_table_spc['SpT'][i])) == 11:
            spt_type = 'L2'
        elif round(float(self.target_table_spc['SpT'][i])) == 13:
            spt_type = 'L2'
        elif round(float(self.target_table_spc['SpT'][i])) == 14:
            spt_type = 'L5'
        elif round(float(self.target_table_spc['SpT'][i])) > 14:
            spt_type = 'L8'
        filt_ = str(self.target_table_spc['Filter_spc'][i])
        if (filt_ == 'z\'') or (filt_ == 'r\'') or (filt_ == 'i\'') or (filt_ == 'g\''):
            filt_ = filt_.replace('\'', '')
        filters = ['I+z', 'z', 'i', 'r']
        filt_idx = 0
        filt_ = filters[filt_idx]

        texp = 0

        while texp < 10:

            if filt_idx > 0:
                print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Change filter to avoid saturation!!')
                filt_idx += 1
                if filt_idx >= 3:
                    print(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' You have to defocus in we want to observe this target')
                    texp = 10.0001
                filt_ = filters[filt_idx]

            if self.telescope == 'Saint-Ex':
                print(day)
                moon_phase = round(moon_illumination(Time(day.iso, out_subfmt='date')), 2)
                if float(self.target_table_spc['J'][i]) != 0.:
                    a = (ETC.etc(mag_val=float(self.target_table_spc['J'][i]), mag_band='J', spt=spt_type,
                                 filt=filt_, airmass=1.1, moonphase=moon_phase, irtf=0.8, num_tel=1,
                                 seeing=1.0, gain=3.48, temp_ccd=-70, observatory_altitude=2780))
                else:
                    if (float(self.target_table_spc['J'][i]) == 0.) and (float(self.target_table_spc['V'][i]) != 0.):
                        a = (ETC.etc(mag_val=float(self.target_table_spc['V'][i]), mag_band='V', spt=spt_type,
                                     filt=filt_, airmass=1.1, moonphase=moon_phase, irtf=0.8, num_tel=1,
                                     seeing=1.0, gain=3.48, temp_ccd=-70, observatory_altitude=2780))
                    else:
                        sys.exit('ERROR: You must precise Vmag or Jmag for this target')
                texp = a.exp_time_calculator(ADUpeak=30000)[0]

            elif self.telescope == 'TN_Oukaimeden':
                if float(self.target_table_spc['J'][i]) != 0.:
                    a = (ETC.etc(mag_val=float(self.target_table_spc['J'][i]), mag_band='J', spt=spt_type,
                                 filt=filt_, airmass=1.1, moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
                else:
                    if (float(self.target_table_spc['J'][i]) == 0.) and (float(self.target_table_spc['V'][i]) != 0.):
                        a = (ETC.etc(mag_val=float(self.target_table_spc['J'][i]), mag_band='J', spt=spt_type,
                                     filt=filt_, airmass=1.1, moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
                    else:
                        sys.exit('ERROR: You must precise Vmag or Jmag for this target')
                texp = a.exp_time_calculator(ADUpeak=50000)[0]
                print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK +
                      ' Don\'t forget to  calculate exposure time for TRAPPIST observations!!')

            elif self.telescope == 'TS_La_Silla':
                if float(self.target_table_spc['J'][i]) != 0.:
                    a = (ETC.etc(mag_val=self.target_table_spc['J'][i], mag_band='J', spt=spt_type, filt=filt_,
                                 airmass=1.1, moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.4, gain=1.1))
                else:
                    if (float(self.target_table_spc['J'][i]) == 0.) and (float(self.target_table_spc['V'][i]) != 0.):
                        a = (ETC.etc(mag_val=self.target_table_spc['V'][i], mag_band='V', spt=spt_type, filt=filt_,
                                     airmass=1.1, moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.4, gain=1.1))
                    else:
                        sys.exit('ERROR: You must precise Vmag or Jmag for ' + str(self.target_table_spc['Sp_ID'][i]))
                texp = a.exp_time_calculator(ADUpeak=50000)[0]
                print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Don\'t forget to  '
                                                               'calculate exposure time for TRAPPIST observations!!')

            elif self.telescope == 'Artemis':
                if float(self.target_table_spc['J'][i]) != 0.:
                    a = (ETC.etc(mag_val=float(self.target_table_spc['J'][i]), mag_band='J', spt=spt_type,
                                 filt=filt_, airmass=1.1, moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
                else:
                    if (float(self.target_table_spc['J'][i]) == 0.) and (float(self.target_table_spc['V'][i]) != 0.):
                        a = (ETC.etc(mag_val=float(self.target_table_spc['V'][i]), mag_band='V', spt=spt_type,
                                     filt=filt_, airmass=1.1, moonphase=0.5, irtf=0.8, num_tel=1, seeing=1.0, gain=1.1))
                    else:
                        sys.exit('ERROR: You must precise Vmag or Jmag for this target')
                texp = a.exp_time_calculator(ADUpeak=45000)[0]

            elif self.telescope == 'Io' or self.telescope == 'Europa' \
                    or self.telescope == 'Ganymede' or self.telescope == 'Callisto':
                if float(self.target_table_spc['J'][i]) != 0.:
                    a = (ETC.etc(mag_val=float(self.target_table_spc['J'][i]), mag_band='J', spt=spt_type,
                                 filt=filt_, airmass=1.1, moonphase=0.5, irtf=0.8, num_tel=1, seeing=0.7, gain=1.1))
                else:
                    if (float(self.target_table_spc['J'][i]) == 0.) and (float(self.target_table_spc['V'][i]) != 0.):
                        a = (ETC.etc(mag_val=float(self.target_table_spc['V'][i]), mag_band='V', spt=spt_type,
                                     filt=filt_, airmass=1.1, moonphase=0.5, irtf=0.8, num_tel=1, seeing=0.7, gain=1.1))
                    else:
                        sys.exit('ERROR: You must precise Vmag or Jmag for this target')
                texp = a.exp_time_calculator(ADUpeak=45000)[0]

            filt_idx += 1

        return texp

    def is_moon_and_visibility_constraint(self, day):
        """

        Parameters
        ----------
        day

        Returns
        -------

        """
        dt_1day = Time('2018-01-02 00:00:00', scale='utc')-Time('2018-01-01 00:00:00', scale='utc')

        start_between_civil_nautical = Time((Time(
            self.observatory.twilight_evening_nautical(day,
                                                       which='next')).value +
                                             Time(self.observatory.twilight_evening_civil(
                                                 day,
                                                 which='next')).value) / 2,
                                            format='jd')

        end_between_nautical_civil = Time((Time(
            self.observatory.twilight_morning_nautical(day+dt_1day,
                                                       which='nearest')).value +
                                           Time(self.observatory.twilight_morning_civil(
                                               day+dt_1day,
                                               which='nearest')).value) / 2,
                                          format='jd')

        self.observability_table_day = observability_table(self.constraints, self.observatory, self.targets,
                                                           time_range=Time([start_between_civil_nautical.iso,
                                                                            end_between_nautical_civil.iso]))

        self.observability_table_day['fraction of time observable'] = \
            self.observability_table_day['fraction of time observable'] * self.night_duration(day).value * 24

        # is_visible_mid_night = is_observable(self.constraints, self.observatory, self.targets,
        #                                     times=Time(start_between_civil_nautical +
        #                                                self.night_duration(day).value/2))

        # idx_not_visible_mid_night = np.where((is_visible_mid_night is False))

        # self.observability_table_day['ever observable'][idx_not_visible_mid_night[0]] = False

        return self.observability_table_day

    def night_duration(self, day):
        """

        Parameters
        ----------
        day

        Returns
        -------

        """
        start_between_civil_nautical = Time((Time(
            self.observatory.twilight_evening_nautical(day,
                                                       which='next')).value +
                                             Time(self.observatory.twilight_evening_civil(
                                                 day,
                                                 which='next')).value) / 2,
                                            format='jd')

        end_between_nautical_civil = Time((Time(
            self.observatory.twilight_morning_nautical(day+1,
                                                       which='nearest')).value +
                                           Time(self.observatory.twilight_morning_civil(
                                               day+1,
                                               which='nearest')).value) / 2,
                                          format='jd')

        dura = end_between_nautical_civil - start_between_civil_nautical
        # dura = Time(Time(self.observatory.twilight_morning_nautical(day + dt_1day ,which='nearest')).jd - \
        #       Time(self.observatory.twilight_evening_nautical(day ,which='next')).jd,format='jd')
        return dura

    def make_schedule(self):

        self.reverse_df1 = reverse_observability(self.observatory, self.targets, self.constraints, self.time_ranges)

        for t in tqdm(range(0, self.date_range_in_days), desc="Scheduling "):
            print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' day is : ', Time(self.date_range[0] + t).iso)
            self.day_of_night = self.date_ranges_day_by_day[t]
            self.observability_selection(self.day_of_night)
            self.night_block = None
            self.moon_and_visibility_constraint_table = self.is_moon_and_visibility_constraint(self.day_of_night)
            for i in range(len(self.target_table_spc)):
                df_prediction = prediction(self.target_table_spc['Sp_ID'][i], self.target_table_spc['RA'][i],
                                           self.target_table_spc['DEC'][i],
                                           self.target_table_spc['T0'][i], self.target_table_spc['P'][i],
                                           self.target_table_spc['W'][i],
                                           self.day_of_night, self.ntr)
                if any(df_prediction[self.observatory_name]):
                    list_obs = ['SSO', 'TS_La_Silla', 'TN_Oukaimeden', 'Saint-Ex', 'SNO']
                    list_obs.remove(self.observatory_name)
                    df_prediction = df_prediction.drop(columns=list_obs)
                    df_prediction = df_prediction.loc[df_prediction[self.observatory_name]]
                    df_prediction['priority'] = [self.priority_table['priority'][i]]*len(df_prediction)
                    self.list_next_transit_window.append(df_prediction)

            self.next_transit_window = pd.concat(self.list_next_transit_window)
            self.next_transit_window = self.next_transit_window.sort_values(by=['priority'],ascending=False)
            self.next_transit_window = self.next_transit_window.reset_index()
            self.next_transit_window = self.next_transit_window.rename(columns={"index": "Sp_ID"})
            # delete line already used
            for j in range(len(self.next_transit_window)):
                if (self.next_transit_window['mid transit'][j] > self.date_range[0]) and \
                        (self.next_transit_window['mid transit'][j] < self.date_range[1]):
                    if self.night_block is not None:
                        new_block = self.transit_follow_up(input_name=self.next_transit_window['Sp_ID'][j])
                        #replace by  any
                        # for k in range(0,len(self.night_block)):
                        #     if new_block is not None:
                        #         if (new_block['end time (UTC)'][k] < self.night_block['start time (UTC)'][k]) or\
                        #                 (new_block['start time (UTC)'][k] > self.night_block['end time (UTC)'][k]):
                        #             self.night_block = pd.concat([new_block.to_pandas(),self.night_block.to_pandas()])

                    else:
                        self.night_block = self.transit_follow_up(input_name=self.next_transit_window['Sp_ID'][j])
                        if self.night_block is not None:
                            self.locking_observations()

            # self.index_prio = np.argsort(self.priority_table['priority'])
            # self.priority_table_ranked = self.priority_table[self.index_prio]
        print()

    def transit_follow_up(self, input_name):
        """
        Function to add a  night block with a transit of a follow candidate in the plans
        Parameters
        ----------
        input_name: name of the follow up candidate to schedule

        Returns
        -------
        Create a night block with transiting candidate observed for: transit duration * (1 + 1.5) centre on T0_predicted

        """
        self.observatory = charge_observatories(self.observatory_name)[0]

        constraints_follow_up = [AltitudeConstraint(min=self.Altitude_constraint * u.deg),
                                 MoonSeparationConstraint(min=self.Moon_constraint * u.deg),
                                 TimeConstraint(Time(self.day_of_night),Time(self.day_of_night+1))]

        if input_name is not None:
            i = int(np.where((self.target_table_spc['Sp_ID'] == input_name))[0])
            blocks = []
            epoch = Time(self.target_table_spc['T0'][i], format='jd')
            period = self.target_table_spc['P'][i] * u.day
            duration = self.target_table_spc['W'][i] * u.day
            oot_time = duration.value * 1.5 * u.day
            T0_err_transit = self.target_table_spc['T0_err'][i]
            P_err_transit = self.target_table_spc['P_err'][i]
            W_err_transit = self.target_table_spc['W_err'][i]
            target_transit = EclipsingSystem(primary_eclipse_time=epoch, orbital_period=period, duration=duration,
                                             name=self.target_table_spc['Sp_ID'][i])
            print(Fore.GREEN + Fore.GREEN + Fore.GREEN + 'INFO: ' + Fore.BLACK + ' ' + Fore.BLACK +
                  Fore.BLACK + str(self.target_table_spc['Sp_ID'][i]) + ' next transit: ')
            print(Time(target_transit.next_primary_eclipse_time(self.day_of_night, n_eclipses=1)).iso)
            timing_to_obs_jd = Time(target_transit.next_primary_eclipse_time(self.day_of_night, n_eclipses=1)).jd
            try:
                ing_egr = target_transit.next_primary_ingress_egress_time(self.day_of_night, n_eclipses=1)
                observable = is_event_observable(constraints_follow_up, self.observatory, self.targets[i],
                                                 times_ingress_egress=ing_egr)
            except ValueError:
                observable = False
                print(Fore.RED + 'ERROR: ' + Fore.BLACK + ' No transit of ', self.target_table_spc['Sp_ID'][i],
                      ' on the period chosen')
                pass

            if np.any(observable):
                err_T0_neg = timing_to_obs_jd[0] - (np.round((timing_to_obs_jd[0] - epoch.jd) / period.value, 1) *
                                                    (period.value - P_err_transit) + (epoch.jd - T0_err_transit))
                err_T0_pos = (np.round((timing_to_obs_jd[0] - epoch.jd) / period.value, 1) *
                              (period.value + P_err_transit) + (epoch.jd + T0_err_transit)) - timing_to_obs_jd[0]
                start_transit = Time(ing_egr[0][0].value - err_T0_neg - oot_time.value - W_err_transit,
                                     format='jd')
                end_transit = Time(ing_egr[0][1].value + err_T0_pos + oot_time.value + W_err_transit,
                                   format='jd')
                dur_obs_transit_oot_target = (end_transit - start_transit).value * 1. * u.day
                
                if (end_transit > Time(
                        self.observatory.twilight_morning_nautical(self.day_of_night + 1, which='nearest'))) \
                        or (start_transit <
                            Time(self.observatory.twilight_evening_nautical(self.day_of_night, which='next'))):
                    # Not the full baseline (transit +OOT) is observable
                    if (Time(ing_egr[0][1]) < Time(
                            self.observatory.twilight_morning_nautical(self.day_of_night + 1, which='nearest'))) \
                            and (Time(ing_egr[0][0]) >
                                 Time(self.observatory.twilight_evening_nautical(self.day_of_night, which='next'))):
                        # Transit is observable partially
                        constraints_transit_target = [AltitudeConstraint(min=self.Altitude_constraint * u.deg),
                                                      TimeConstraint(start_transit, end_transit),
                                                      MoonSeparationConstraint(min=self.Moon_constraint * u.deg),
                                                      AtNightConstraint(max_solar_altitude=-12 * u.deg)]
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ',
                              str(self.target_table_spc['Sp_ID'][i]), ' : ', start_transit.iso)
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit of ',
                              str(self.target_table_spc['Sp_ID'][i]), ' : ', end_transit.iso)
                        print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Time out of transit not optimal.')

                        a = ObservingBlock(self.targets[i], dur_obs_transit_oot_target, -1,
                                           constraints=constraints_transit_target,
                                           configuration={"filt": str(self.target_table_spc['Filter_spc'][i]),
                                                          "texp":  str(self.exposure_time(self.day_of_night, i))})
                        blocks.append(a)
                        transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
                        seq_schedule_ss1 = Schedule(self.day_of_night, self.day_of_night + 1)
                        sequen_scheduler_ss1 = SPECULOOSScheduler(constraints=constraints_transit_target,
                                                                  observer=self.observatory,
                                                                  transitioner=transitioner)
                        sequen_scheduler_ss1(blocks, seq_schedule_ss1)

                        if len(seq_schedule_ss1.to_table()['target']) == 0:
                            print(Fore.RED + 'ERROR:  ' + Fore.BLACK +
                                  ' The moon is too closed for the transit to be observed')
                            night_block = None

                    else:  # transit is full but not the baseline
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ',
                              str(self.target_table_spc['Sp_ID'][i]), ' : ', Time(ing_egr[0][0]).iso)
                        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' end_transit of ',
                              str(self.target_table_spc['Sp_ID'][i]), ' : ', Time(ing_egr[0][1]).iso)
                        print(Fore.YELLOW + 'WARNING: ' + Fore.BLACK + ' Transit not full.')
                        constraints_transit_target = [AltitudeConstraint(min=self.Altitude_constraint * u.deg),
                                                      MoonSeparationConstraint(min=self.Moon_constraint * u.deg),
                                                      TimeConstraint(Time(ing_egr[0][0]), Time(ing_egr[0][1]))]
                        a = ObservingBlock(self.targets[i], dur_obs_transit_oot_target, -1,
                                           constraints=constraints_transit_target,
                                           configuration={"filt": str(self.target_table_spc['Filter_spc'][i]),
                                                          "texp":  str(self.exposure_time(self.day_of_night, i))})
                        blocks.append(a)
                        transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
                        seq_schedule_ss1 = Schedule(self.day_of_night, self.day_of_night + 1)
                        sequen_scheduler_ss1 = SPECULOOSScheduler(constraints=constraints_transit_target,
                                                                  observer=self.observatory,
                                                                  transitioner=transitioner)
                        sequen_scheduler_ss1(blocks, seq_schedule_ss1)

                        if len(seq_schedule_ss1.to_table()['target']) == 0:
                            print(Fore.RED + 'ERROR:  ' + Fore.BLACK +
                                  ' The moon is too closed for the transit to be observed')
                            night_block = None

                else:  # full baseline is observable
                    constraints_transit_target = [AltitudeConstraint(min=self.Altitude_constraint * u.deg),
                                                  MoonSeparationConstraint(min=self.Moon_constraint * u.deg),
                                                  TimeConstraint(start_transit, end_transit)]
                    a = ObservingBlock(self.targets[i], dur_obs_transit_oot_target, -1,
                                       constraints=constraints_transit_target,
                                       configuration={"filt": str(self.target_table_spc['Filter_spc'][i]),
                                                      "texp": str(self.exposure_time(self.day_of_night, i))})
                    blocks.append(a)
                    transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
                    seq_schedule_ss1 = Schedule(self.day_of_night, self.day_of_night + 1)
                    sequen_scheduler_ss1 = SPECULOOSScheduler(constraints=constraints_transit_target,
                                                              observer=self.observatory,
                                                              transitioner=transitioner)
                    sequen_scheduler_ss1(blocks, seq_schedule_ss1)

                    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' start_transit of ',
                          str(self.target_table_spc['Sp_ID'][i]), ' : ', Time(ing_egr[0][0]).iso)
                    print(Fore.GREEN + Fore.GREEN + 'INFO: ' + Fore.BLACK + ' ' + Fore.BLACK + ' end_transit of ',
                          str(self.target_table_spc['Sp_ID'][i]), ' : ', Time(ing_egr[0][1]).iso)
                    print(Fore.GREEN + Fore.GREEN + 'INFO: ' + Fore.BLACK + ' ' + Fore.BLACK +
                          ' Transit is expected to be full.')
                    if len(seq_schedule_ss1.to_table()['target']) == 0:
                        print(Fore.RED + 'ERROR:  ' + Fore.BLACK +
                              ' The moon is too closed for the transit to be observed')
                        night_block = None

                if len(seq_schedule_ss1.to_table()['target']) != 0:
                    night_block = seq_schedule_ss1.to_table()

            else:
                night_block = None
                print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' no transit of ', self.target_table_spc['Sp_ID'][i], 
                      ' this day or does not respect constraint Moon separation')
            return night_block

    def locking_observations(self):
        night_block = self.night_block.to_pandas()

        night_block.to_csv(path_spock + '/DATABASE/' + self.telescope + '/' +
                           'Locked_obs/' + 'lock_night_block_' + self.telescope + '_' +
                           Time(self.day_of_night.iso, out_subfmt='date').iso + '.txt', index=None)

        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + 'Observation block saved as ' +
              path_spock + '/DATABASE/' + self.telescope + '/' +
              'Locked_obs/' + 'lock_night_block_' + self.telescope + '_' +
              Time(self.day_of_night.iso, out_subfmt='date').iso + '.txt')


def prediction(name, ra, dec, timing, period, duration, start_date, ntr):

    start_date = Time(start_date)

    constraints_predictions = [AltitudeConstraint(min=25 * u.deg), AtNightConstraint()]

    target_transit = EclipsingSystem(primary_eclipse_time=Time(timing, format='jd'),
                                     orbital_period=period*u.day, duration=duration*u.day,
                                     name=name)

    mid_transit_timing = Time(target_transit.next_primary_eclipse_time(start_date, n_eclipses=ntr)).iso
    mid_transit_timing_jd = Time(target_transit.next_primary_eclipse_time(start_date, n_eclipses=ntr)).jd

    ing_egr = target_transit.next_primary_ingress_egress_time(start_date, n_eclipses=ntr)

    target = [FixedTarget(coord=SkyCoord(ra=ra * u.degree, dec=dec * u.degree), name=name)]

    observable_sso = is_event_observable(constraints_predictions, charge_observatories('SSO')[0], target,
                                         times_ingress_egress=ing_egr)
    observable_sno = is_event_observable(constraints_predictions, charge_observatories('SNO')[0], target,
                                         times_ingress_egress=ing_egr)
    observable_saintex = is_event_observable(constraints_predictions, charge_observatories('Saint-Ex')[0], target,
                                             times_ingress_egress=ing_egr)
    observable_tn_oukaimeden = is_event_observable(constraints_predictions, charge_observatories('TN_Oukaimeden')[0],
                                                   target, times_ingress_egress=ing_egr)
    observable_ts_la_silla = is_event_observable(constraints_predictions, charge_observatories('TS_La_Silla')[0],
                                                 target, times_ingress_egress=ing_egr)

    df = pd.DataFrame(data=ing_egr.sort().iso, index=[name]*len(ing_egr), columns=["Ingress", "Egress"])

    df['mid transit'] = mid_transit_timing

    df['mid transit JD'] = np.round(mid_transit_timing_jd, 3)

    df['SSO'] = observable_sso[0]

    df['SNO'] = observable_sno[0]

    df['Saint-Ex'] = observable_saintex[0]

    df['TN_Oukaimeden'] = observable_tn_oukaimeden[0]

    df['TS_La_Silla'] = observable_ts_la_silla[0]

    return df
