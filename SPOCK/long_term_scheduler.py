#code for long term scheduling for SPECULOOS
import numpy as np
import matplotlib
from datetime import date , datetime , timedelta
from astropy.table import Table, Column
from astropy.table import join
from astropy.table import unique
from astropy import units as u
from astropy.coordinates import SkyCoord, get_sun, AltAz, EarthLocation
from astropy.utils import iers
iers.IERS_A_URL = 'http://toshi.nofs.navy.mil/ser7/finals2000A.all'
from astroplan import download_IERS_A
download_IERS_A()
from astroplan import FixedTarget, AltitudeConstraint, MoonSeparationConstraint,AtNightConstraint,observability_table, is_observable, months_observable,time_grid_from_range,LocalTimeConstraint, is_always_observable
#from eScheduler.spe_schedule import SPECULOOSScheduler, PriorityScheduler, Schedule, ObservingBlock, Transitioner, DateElsa
from astroplan import TimeConstraint
from astroplan import Observer
from astropy.time import Time
from astropy.table import Table, Column, hstack, vstack, unique, join
import pandas as pd
import functools
from functools import reduce
import os
import urllib.request
import requests
import re
import yaml
# from astroplan import download_IERS_A
# download_IERS_A()

# definition of the different fonctions
def Diff_list(li1, li2):
    """
    Inform on the difference between two lists

    Parameters
    ----------
    li1: list numero 1
    li2: list numero 2

    Returns
    -------
    Elements than are in list 1 but not in list 2
    """

    return (list(set(li1) - set(li2)))

def compare_target_lists(path_target_list):
    """
    Compare the target list from the given folder to the one on STARGATE and Cambridge server
    If different trigger a warning a tell how many targets are actually different from the referenced target list
    Parameters
    ----------
    path_target_list: path on your computer toward the target list, by default take the one on the Cambridge server

    Returns
    -------
    An idication about if the target list is the referenced one or not

    """
    TargetURL="http://www.mrao.cam.ac.uk/SPECULOOS/target_list_gaia.csv"
    user, password = 'educrot', '58JMSGgdmzTB'
    resp = requests.get(TargetURL, auth=(user, password))
    open('target_list_gaia.csv', 'wb').write(resp.content)
    content=resp.text.replace("\n", "")
    df_target_list_gaia = pd.read_csv('target_list_gaia.csv', delimiter=',',skipinitialspace=True,error_bad_lines=False)
    df_user = pd.read_csv(path_target_list, delimiter=' ',skipinitialspace=True,error_bad_lines=False)
    if Diff_list(df_user['Name'],df_target_list_gaia['spc']):
        print('WARNING ! Tragets in User\'s list but not in Cambridge server\'s list: ',Diff_list(df_user['Name'],df_target_list_gaia['spc'] ))
    if Diff_list(df_target_list_gaia['spc'],df_user['Name'] ):
        print('WARNING ! Tragets in Cambridge server\'s list but not in User\'s list: ',Diff_list(df_target_list_gaia['spc'],df_user['Name'] ))
    else:
        print('OK ! User\'s list is similar to the one on the Cambridge server')

def update_hours_target_list(path_target_list):
    """
    Update the targets hours of observation from Cambridge server

    Parameters
    ----------
    path_target_list: path on your computer toward the target list, by default take the one on the Cambridge server

    Returns
    -------
    Same target list but with updated number of hours of observation

    """
    TargetURL="http://www.mrao.cam.ac.uk/SPECULOOS/reports/SurveyTotal"
    user, password = 'educrot', '58JMSGgdmzTB'
    resp = requests.get(TargetURL, auth=(user, password))
    open('SurveyTotal.txt', 'wb').write(resp.content)
    content=resp.text.replace("\n", "")
    df_cambridge = pd.read_csv('SurveyTotal.txt', delimiter=' ',skipinitialspace=True,error_bad_lines=False)
    df_user = pd.read_csv(path_target_list, delimiter=' ',skipinitialspace=True,error_bad_lines=False)
    index_df2 = np.where([df_cambridge['Target']==target_df2 for target_df2 in df_user['Name']])
    index_df = np.where([df_user['Name']==target_df for target_df in df_cambridge['Target']])
    t = Table.from_pandas(df_cambridge)
    t2 = Table.from_pandas(df_user)
    for i in range(0,len(t['Hours'][index_df[0]])):
        t2['nb_hours_surved'][index_df2[0][i]]=max(float(t['Hours'][index_df[0][i]]),float(t2['nb_hours_surved'][index_df2[0][i]]))
    df2=pd.DataFrame(t2.to_pandas())
    df2.to_csv(path_target_list,sep=' ',index=False)
    print('OK ! Number of hours of observation have been updated from lastest version on Cambridge server')

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

    target_table_spc = Table.read(path_target_list, format='ascii')
    ra_spc = {}
    dec_spc = {}
    targets = []
    for i in range(0,len(target_table_spc)):
        if target_table_spc['DEC1'][i]<0:
            ra_spc[i]=15*(float(target_table_spc['RA1'][i])+(float(target_table_spc['RA2'][i])+float(target_table_spc['RA3'][i])/60)/60)
            dec_spc[i]=float(target_table_spc['DEC1'][i])-(float(target_table_spc['DEC2'][i])+float(target_table_spc['DEC3'][i])/60)/60
            targets.append(FixedTarget(coord=SkyCoord(ra=15*(target_table_spc['RA1'][i]+(target_table_spc['RA2'][i]+target_table_spc['RA3'][i]/60)/60)*u.deg, dec=(target_table_spc['DEC1'][i]-(target_table_spc['DEC2'][i]+target_table_spc['DEC3'][i]/60)/60)*u.deg), name=target_table_spc['Name'][i]))
        else:
            ra_spc[i]=15*(float(target_table_spc['RA1'][i])+(float(target_table_spc['RA2'][i])+float(target_table_spc['RA3'][i])/60)/60)
            dec_spc[i]=float(target_table_spc['DEC1'][i])-(float(target_table_spc['DEC2'][i])+float(target_table_spc['DEC3'][i])/60)/60
            targets.append(FixedTarget(coord=SkyCoord(ra=15*(target_table_spc['RA1'][i]+(target_table_spc['RA2'][i]+target_table_spc['RA3'][i]/60)/60)*u.deg, dec=(target_table_spc['DEC1'][i]+(target_table_spc['DEC2'][i]+target_table_spc['DEC3'][i]/60)/60)*u.deg), name=target_table_spc['Name'][i]))
    return(targets)


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
        observatories.append(Observer(location=location_saintex, name="saintex", timezone="UTC"))

    if 'TS La Silla' in str(Name):
        location_TSlasilla = EarthLocation.from_geodetic(-70.73000000000002*u.deg, -29.25666666666666*u.deg, 2346.9999999988418*u.m)
        observatories.append(Observer(location=location_TSlasilla, name="TSlasilla", timezone="UTC"))

    if 'TN Oukaimeden' in str(Name):
        location_TNOuka = EarthLocation.from_geodetic(31.20516*u.deg, -7.862263*u.deg, 2751*u.m)
        observatories.append(Observer(location=location_TNOuka, name="TNOuka", timezone="UTC"))

    if 'Munich' in str(Name):
        location_munich= EarthLocation.from_geodetic(48.2*u.deg, -11.6*u.deg, 600*u.m)
        observatories.append(Observer(location=location_munich, name="Munich", timezone="UTC"))

    return observatories

def _generate_24hr_grid(t0, start, end, N, for_deriv=False):
    """
    Generate a nearly linearly spaced grid of time durations.
    The midpoints of these grid points will span times from ``t0``+``start``
    to ``t0``+``end``, including the end points, which is useful when taking
    numerical derivatives.
    Parameters
    ----------
    t0 : `~astropy.time.Time`
        Time queried for, grid will be built from or up to this time.
    start : float
        Number of days before/after ``t0`` to start the grid.
    end : float
        Number of days before/after ``t0`` to end the grid.
    N : int
        Number of grid points to generate
    for_deriv : bool
        Generate time series for taking numerical derivative (modify
        bounds)?
    Returns
    -------
    `~astropy.time.Time`
    """

    if for_deriv:
        time_grid = np.concatenate([[start - 1/(N-1)],
                                    np.linspace(start, end, N)[1:-1],
                                    [end + 1/(N-1)]])*u.day
    else:
        time_grid = np.linspace(start, end, N)*u.day

    # broadcast so grid is first index, and remaining shape of t0
    # falls in later indices. e.g. if t0 is shape (10), time_grid
    # will be shape (N, 10). If t0 is shape (5, 2), time_grid is (N, 5, 2)
    while time_grid.ndim <= t0.ndim:
        time_grid = time_grid[:, np.newaxis]
    # we want to avoid 1D grids since we always want to broadcast against targets
    if time_grid.ndim == 1:
        time_grid = time_grid[:, np.newaxis]
    return t0 + time_grid

#altaz: function defined in Astroplan
def altaz(self, time, target=None, obswl=None, grid_times_targets=False):
    """
    Get an `~astropy.coordinates.AltAz` frame or coordinate.
    If ``target`` is None, generates an altitude/azimuth frame. Otherwise,
    calculates the transformation to that frame for the requested ``target``.
    Parameters
    ----------
    time : `~astropy.time.Time` or other (see below)
      The time at which the observation is taking place. Will be used as
      the ``obstime`` attribute in the resulting frame or coordinate. This
      will be passed in as the first argument to the `~astropy.time.Time`
      initializer, so it can be anything that `~astropy.time.Time` will
      accept (including a `~astropy.time.Time` object)
    target : `~astroplan.FixedTarget`, `~astropy.coordinates.SkyCoord`, or list (optional)
      Celestial object(s) of interest. If ``target`` is `None`, returns
      the `~astropy.coordinates.AltAz` frame without coordinates.
    obswl : `~astropy.units.Quantity` (optional)
      Wavelength of the observation used in the calculation.
    grid_times_targets: bool (optional)
      If True, the target object will have extra dimensions packed
      onto the end, so that calculations with M targets and N times
      will return an (M, N) shaped result. Otherwise, we rely on
      broadcasting the shapes together using standard numpy
      rules. Useful for grid searches for rise/set times etc.
    Returns
    -------
    `~astropy.coordinates.AltAz`
      If ``target`` is `None`, returns `~astropy.coordinates.AltAz` frame.
      If ``target`` is not `None`, returns the ``target`` transformed to
      the `~astropy.coordinates.AltAz` frame.
    Examples
    --------
    Create an instance of the `~astropy.coordinates.AltAz` frame for an
    observer at Apache Point Observatory at a particular time:
    >>> from astroplan import Observer
    >>> from astropy.time import Time
    >>> from astropy.coordinates import SkyCoord
    >>> apo = Observer.at_site("APO")
    >>> time = Time('2001-02-03 04:05:06')
    >>> target = SkyCoord(0*u.deg, 0*u.deg)
    >>> altaz_frame = apo.altaz(time)
    Now transform the target's coordinates to the alt/az frame:
    >>> target_altaz = target.transform_to(altaz_frame) # doctest: +SKIP
    Alternatively, construct an alt/az frame and transform the target to
    that frame all in one step:
    >>> target_altaz = apo.altaz(time, target) # doctest: +SKIP
    """
    if target is not None:
        time, target = self._preprocess_inputs(time, target, grid_times_targets)

        altaz_frame = AltAz(location=self.location, obstime=time,
                      pressure=self.pressure, obswl=obswl,
                      temperature=self.temperature,
                      relative_humidity=self.relative_humidity)
    if target is None:
      # Return just the frame
      return altaz_frame
    else:
        return target.transform_to(altaz_frame)

def type_spec(x):
    if x=='M6':
        y=210
    if x=='M7':
        y=160
    if x=='M8':
        y=130
    if x=='M9':
        y=100
    if x=='M10':
        y=90
    else:
        y=90
    return y


def Observability(j,time_range,observatory,targets,constraints):
    """
    Give a table with the observability score for each target of targets
    regarding the constraints given and for all ranges of time_range
    Parameters
    ----------
        j : list [start end], element of time range (that is to say, month of the year)
        time_range : List of astropy.Time range with start and end times
        observatory : observatory chosen
        targets : target list on the FixedTarget() format from astroplan
        constraints : general constraints for a target to be shceduled

    Returns
    -------
        Observability table: 12 columns (target name and each element of time_range, e.g months),
        rows= nb of targets
    """
    targets_observable=[]
    observable=np.argwhere(is_observable(constraints,observatory,targets,time_range=time_range))

    ##WARNING: Need to be replace by a np.where, but I can't find how
    [targets_observable.append(targets[int(obs)]) for obs in observable]

    table_observability=observability_table(constraints, observatory, targets_observable, time_range=time_range, time_grid_resolution = 0.5*u.hour) #calcul un edeuxieme fois is observable
    table_observability.remove_column('always observable')
    table_observability.remove_column('ever observable')
    table_observability.rename_column('fraction of time observable', 'Month' + str(j))

    return table_observability

def reverse_Observability(observatory,targets,constraints,time_ranges):
    """
    Reverse observability table, rows become columns
    Parameters
    ----------
        observatory : observatory chosen
        time_range : List of astropy.Time range with start and end times
        targets : target list on the FixedTarget() format from astroplan
        constraints : general constraints for a target to be shceduled

    Returns
    -------
        observability table with no NaN value (0 instead) inversed with targets as columns
        and elements of time_ranges (months) as rows
    """
    start_fmt = Time(time_ranges[0][0].iso , out_subfmt = 'date').iso
    end_fmt =  Time(time_ranges[len(time_ranges)-1][1].iso , out_subfmt = 'date').iso

    if os.path.exists('reverse_Obs_' + str(observatory.name) + '_' +  start_fmt + '_' + end_fmt + '_'  + str(len(targets)) + '.csv'):
        name_file = 'reverse_Obs_' + str(observatory.name) + '_' +  start_fmt + '_' + end_fmt + '_'  + str(len(targets)) +  '.csv'
        reverse_df1 = pd.read_csv(name_file, delimiter = ',')
        return reverse_df1

    else:
        #time_ranges = [Time(['2019-01-01 12:00:00', '2019-01-31 12:00:00']),Time(['2019-02-01 12:00:00', '2019-02-28 12:00:00']),Time(['2019-03-01 15:00:00', '2019-03-31 15:00:00']),Time(['2019-04-01 15:00:00', '2019-04-30 15:00:00']),Time(['2019-05-01 15:00:00', '2019-05-31 15:00:00']),Time(['2019-06-01 15:00:00', '2019-06-30 15:00:00']),Time(['2019-07-01 12:00:00', '2019-07-31 12:00:00']),Time(['2019-08-01 12:00:00', '2019-08-31 12:00:00']),Time(['2019-09-01 12:00:00', '2019-09-30 12:00:00']),Time(['2019-10-01 12:00:00', '2019-10-31 12:00:00']),Time(['2019-11-01 12:00:00', '2019-11-30 12:00:00']),Time(['2019-12-01 12:00:00', '2019-12-31 12:00:00'])]
        tables_observability=list(map((lambda x: Observability(x,time_ranges[x],observatory,targets,constraints)), range(0,len(time_ranges))))

        df=list(map((lambda x: pd.DataFrame(tables_observability[x].to_pandas())),range(0,len(time_ranges))))
        a=reduce(functools.partial(pd.merge,how='outer', on='target name'), df)
        df=a.replace(to_replace=float('NaN'), value=0.0)
        df1 = df.set_index('target name')
        reverse_df1=df1.T
        reverse_df1.to_csv('reverse_Obs_' + str(observatory.name) + '_' +  start_fmt + '_' + end_fmt + '_'  + str(len(targets)) + '.csv', sep= ',')
        return reverse_df1



def month_option(target_name,reverse_df1):
    """
        create a list of the best month for oservation for each target

    Parameters
    ----------
        target_name: name of the target
        reverse_df1: observability table with no NaN value (0 instead) inversed with targets as columns
        and elements of time_ranges (months) as rows

    Returns
    -------
        month: a list with the best month to observe the target
        month_2nd_option: same but for the second best month
        months_3rd_option: same but for the third best month
        etc

    Remarks
    -------
        the 2nd, 3rd etc choices are here in case the target list is not long enough to give
        solutions for all the telescopes each months, allows to avoid blancks in observations
    """
    months=[]
    months_2nd_option=[]
    months_3rd_option=[]
    months_4th_option=[]
    months_5th_option=[]
    try:
        months.append(reverse_df1[str(target_name)].idxmax())
        months_2nd_option.append(reverse_df1[str(target_name)].nlargest(2,keep='first').index[1])
        months_3rd_option.append(reverse_df1[str(target_name)].nlargest(3,keep='first').index[2])
        months_4th_option.append(reverse_df1[str(target_name)].nlargest(4,keep='first').index[3])
        months_5th_option.append(reverse_df1[str(target_name)].nlargest(5,keep='first').index[4])

    except KeyError:
        months.append('no')
        months_2nd_option.append('no')
        months_3rd_option.append('no')
        months_4th_option.append('no')
        months_5th_option.append('no')
    return {'months':months, 'months_2nd_option':months_2nd_option ,'months_3rd_option':months_3rd_option,'months_4th_option':months_4th_option,'months_5th_option':months_5th_option }



def Ranking_month(priority,i,date_range_in_days,observatory,month_obs,months,months_2nd_option,months_3rd_option,months_4th_option,months_5th_option,date_range,targets,name_spc,priospe_spc,reverse_df1,nb_hours_observed,nb_hours_threshold,Kmag_spc,Sptype_spc):
    """
        Rank all targets with respect to their altitude (whether they are up or not),
        also if their best month score is in agreement with the desired month of observation,
        if their Kmag is inferior to 10.5, if the priospe (equivalent to JWST SNR) is high,
        and if the numer of hours oserved is lower than the threshold (50 hours for the present strategy)

    Parameters
    ----------
        priority: astopy table, gives the observation priority, whether the target is
        more a set target (up a sunset for the whole month) a rise target (up at sunrise for the whole month)
        or both, the altitude at sunset at the beginning of the month,
        the altitude at sunrise at the beginning of the month, the altitude at sunset at the end of the month,
        the altitude at sunrise at the end of the month
        i: int, id target
        date_range_in_days : nb of days in the date range
        observatory: the observatory
        month obs: int, the number of the month of observation, among 0 to 11 (January to December)
        months: the best score of each target for each month
        months_2nd_option: the second best score of each target for each month
        months_3rd_option: the third best score of each target for each month
        months_4th_option: the fourth best score of each target for each month
        months_5th_option: the fifth best score of each target for each month
        time_ranges: list of astropy.Time range with start and end times
        targets : target list on the FixedTarget() format from astroplan
        name_spc: all the name of the targets in the target list
        priospe_spc: SNR JWST * 100, for the targets in the target list
        reverse_df1: inverse of observaility table
        nb_hours_observed: number of hours observed for each target of the target list
        nb_hours_threshold: number of hours to complete for each target of the target list
        Kmag_spc: Kmagnitude of each target of the target list

    Returns
    -------
        priority: the same priority as given in function arguments with updated first columns
        with values:
        -0.5 for non observable targets, <0 for targets where the numer of hours observed > nb hours nb_hours_threshold
        >0 for the others, the highest value belonging to the highest priority target

    """

    # if Sptype_spc[i]=='M6':
    #     factor_spctype=1
    # else:
    #     factor_spctype=2

    if priospe_spc[i]>=100:
        factor_priospe_spc=4
    else:
        factor_priospe_spc=1


    if (int(nb_hours_observed[i])>0) and (int(nb_hours_observed[i])<int(nb_hours_threshold[i])):
        factor_on_going=10**(1+1/(nb_hours_threshold[i]-nb_hours_observed[i]))
    if (int(nb_hours_observed[i])==0):
        factor_on_going=10**(1/(nb_hours_threshold[i]-nb_hours_observed[i]))
    if (int(nb_hours_observed[i])>=int(nb_hours_threshold[i])):
        factor_on_going=-1

    if months[i]=='Month' + str(month_obs) :
        #times = _generate_24hr_grid(date_range[0], 0, date_range_in_days, date_range_in_days)
        times = _generate_24hr_grid(date_range[0], 0, 1, 1)
        priority.add_row(((priospe_spc[i])**factor_priospe_spc*(2**(10*(reverse_df1[name_spc[i]][month_obs])))
        *max(observatory.altaz(times, targets[i], grid_times_targets=True).alt.value)*factor_on_going,name_spc[i],'None',
        observatory.altaz(observatory.twilight_evening_nautical(date_range[0],which='nearest'),targets[i]).alt.value,
        observatory.altaz(observatory.twilight_morning_nautical(date_range[0],which='next'),targets[i]).alt.value,
        observatory.altaz(observatory.twilight_evening_nautical(date_range[1],which='nearest'),targets[i]).alt.value,
        observatory.altaz(observatory.twilight_morning_nautical(date_range[1],which='next'),targets[i]).alt.value))

    else:
        if months_2nd_option[i]=='Month' + str(month_obs) or months_3rd_option[i]=='Month' + str(month_obs) or months_4th_option[i]=='Month' + str(month_obs) or months_5th_option[i]=='Month' + str(month_obs) :
            #times = _generate_24hr_grid(date_range[0], 0, date_range_in_days,date_range_in_days)
            times = _generate_24hr_grid(date_range[0], 0, 1, 1)
            priority.add_row(((priospe_spc[i])**factor_priospe_spc*(2**(10*(reverse_df1[name_spc[i]][month_obs])))
            *max(observatory.altaz(times, targets[i], grid_times_targets=True).alt.value)*factor_on_going,name_spc[i],'None',
            observatory.altaz(observatory.twilight_evening_nautical(date_range[0],which='nearest'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_morning_nautical(date_range[0],which='next'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_evening_nautical(date_range[1],which='nearest'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_morning_nautical(date_range[1],which='next'),targets[i]).alt.value))

        else:
            priority.add_row((-0.5,name_spc[i],'None',
            observatory.altaz(observatory.twilight_evening_nautical(date_range[0],which='nearest'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_morning_nautical(date_range[0],which='next'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_evening_nautical(date_range[1],which='nearest'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_morning_nautical(date_range[1],which='next'),targets[i]).alt.value))

    return priority



class scheduling(object):
    '''This classe gather all the function tools for scheduling the chosen targets'''

    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day
    constraints_ok=False
    def __init__(self, time_ranges,observatory,priority,targets,month_obs,nb_tot_telescopes,tel,filt_spc,tesxpspe_spc,nb_hours_observed, nb_hours_threshold,nb_hours,constraints_ok,idx_first_target,idx_second_target):
        self.time_ranges = time_ranges
        self.nb_hours_threshold = nb_hours_threshold
        self.observatory = observatory
        self.priority = priority
        self.targets = targets
        self.month_obs = month_obs
        self.nb_tot_telescopes=nb_tot_telescopes
        self.tel = tel
        self.filt_spc=filt_spc
        self.tesxpspe_spc=tesxpspe_spc
        self.nb_hours_observed=nb_hours_observed
        self.nb_hours=nb_hours
        self.idx_first_target=idx_first_target
        self.idx_second_target=idx_second_target

        self.constraints_ok=constraints_ok

    def idx_targets(self):
        """
            Give the index of the first and second targets as well as the
            corresponding row in the priority table

        Parameters
        ----------

        Returns
        -------
            idx_first_target: index first target in target list
            first_target: row number idx_first_target in priority table
            idx_second_target: index second target in target list
            second_target: row number idx_second_target in priority table

        """
        index_prio=np.argsort(self.priority['priority'])
        idx=1
        first_target=self.priority[index_prio[-idx]]
        print('------------- targets PRIO ----------------',self.priority,self.targets)
        self.idx_first_target=index_prio[-idx]
        idx_2nd=2
        dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day
        if first_target['set or rise'] == 'rise':
            for i in range(2,len(index_prio)):
                if self.priority['set or rise'][index_prio[-i]]=='set':
                    print('>>>>>>>>>>>>>>2nd suppose to be set >>>>>>>>>>>>>>>>>>>>>JAFFICHE ICICIIIIIII ',self.targets[index_prio[-i]])
                    print(self.observatory.target_rise_time(self.time_ranges[self.month_obs][0],self.targets[index_prio[-i]],which='nearest',horizon=24*u.deg),self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0],which='nearest'))
                    print(self.observatory.target_set_time(self.time_ranges[self.month_obs][0],self.targets[index_prio[-i]],which='next',horizon=24*u.deg),self.observatory.target_rise_time(self.time_ranges[self.month_obs][0],self.targets[self.idx_first_target],which='nearest',horizon=24*u.deg))
                    if (self.observatory.target_rise_time(self.time_ranges[self.month_obs][0],self.targets[index_prio[-i]],which='nearest',horizon=24*u.deg) < self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0],which='nearest')) and \
                    (self.observatory.target_set_time(self.time_ranges[self.month_obs][0],self.targets[index_prio[-i]],which='next',horizon=24*u.deg) > self.observatory.target_rise_time(self.time_ranges[self.month_obs][0],self.targets[self.idx_first_target],which='nearest',horizon=24*u.deg)):
                        print('>>>>>>>>>>>>>>2nd suppose to be set >>>>>>>>>>>>second step ',self.targets[index_prio[-i]])
                        self.idx_second_target=index_prio[-i]
                        second_target=self.priority[self.idx_second_target]
                        idx_2nd=i
                        break
                if self.priority['set or rise'][index_prio[-i]]=='both':
                    idx_2nd=i
        if first_target['set or rise'] == 'set':
            for i in range(2,len(index_prio)):
                if self.priority['set or rise'][index_prio[-i]]=='rise':
                    print('>>>>>>>>>>>>>>>2nd suppose to be rise>>>>>>>>>>>>>>>>>>>>JAFFICHE ICICIIIIIII ',self.targets[index_prio[-i]])
                    print(self.observatory.target_set_time(self.time_ranges[self.month_obs][0],self.targets[index_prio[-i]],which='next',horizon=24*u.deg),self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt,which='nearest'))
                    print(self.observatory.target_rise_time(self.time_ranges[self.month_obs][0],self.targets[index_prio[-i]],which='nearest',horizon=24*u.deg),self.observatory.target_set_time(self.time_ranges[self.month_obs][0],self.targets[self.idx_first_target],which='nearest',horizon=24*u.deg))

                    if (self.observatory.target_set_time(self.time_ranges[self.month_obs][0],self.targets[index_prio[-i]],which='next',horizon=24*u.deg) > self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt,which='nearest')) and \
                    (self.observatory.target_rise_time(self.time_ranges[self.month_obs][0],self.targets[index_prio[-i]],which='nearest',horizon=24*u.deg) < self.observatory.target_set_time(self.time_ranges[self.month_obs][0],self.targets[self.idx_first_target],which='nearest',horizon=24*u.deg)):
                        print('>>>>>>>>>>>>>>2nd suppose to be rise >>>>>>>>>>>>second step ',self.targets[index_prio[-i]])
                        self.idx_second_target=index_prio[-i]
                        second_target=self.priority[self.idx_second_target]
                        idx_2nd=i
                        break
                    if self.priority['set or rise'][index_prio[-i]]=='both':
                        idx_2nd=i
        if first_target['set or rise'] == 'both':
            second_target=first_target
            self.idx_second_target=self.idx_first_target
        print('IN idx_targets','idx_first_target',self.idx_first_target,'idx_second_target',self.idx_second_target)
        # return first_target,second_target

    def is_constraints_met_first_target(self,t,idx_set_targets_sorted,idx_rise_targets_sorted,index_prio):
        """
            Useful when the moon constrain (< 30°) or the hours constraint (nb_hours_observed>nb_hours_threshold)
            are not fullfilled anymore and the first target needs to be changed

        Parameters
        ----------
            t: int, days in the month
            idx_first_target: int, index associated to the first target (with regards to the targets list)
            idx_set_targets_sorted: list of int, index of all the set target ranked y priority
            idx_rise_targets_sorted: list of int, index of all the rise target ranked y priority
            index_prio: list of int, index of all the rise target ranked y priority

        Returns
        -------
            is_moon_constraint_met_first_target: boolean, say if moon constraint is fullfilled on first target
            hours_constraint: boolean, say if hour constraint is fullfilled on first target
            idx_first_target: int, index for the first target (if constraints where already fullfilled the Value
            is the same as input, else it is a new value for a new target)

        """
        dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')
        first_target=self.priority[self.idx_first_target]
        #dt_nautical_civil_evening_2=(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next')-self.observatory.twilight_evening_civil(self.time_ranges[self.month_obs][0]+dt*t,which='next'))/2
        #dt_civil_nautical_morning_2=(self.observatory.twilight_morning_civil(self.time_ranges[self.month_obs][0]+dt*(t+1),which='next')-self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='next'))/2
        night_duration=((Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest')))-(Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))))# print(night_duration.iso)
        moon_idx_set_target=0
        moon_idx_rise_target=0
        print('sunset/sunrise in 1rst target :',Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next')).iso,Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest')).iso)

        constraints_all = [AltitudeConstraint(min=24*u.deg), MoonSeparationConstraint(min=35*u.deg), \
        TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))), \
        (Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest') )))]
        is_moon_constraint_met_first_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_first_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
        hours_constraint_first=self.nb_hours_ok(self.idx_first_target,self.nb_hours_observed)
        idx_safe=1
        idx_safe_1rst=1
        idx_safe_2nd=1
        print('1ere target is:',first_target)
        print('constraints 1ere:',is_moon_constraint_met_first_target, hours_constraint_first)
        while not (is_moon_constraint_met_first_target & hours_constraint_first):

            before_change_first_target=self.priority[self.idx_first_target]
            before_change_idx_first_target=self.idx_first_target

            print('------lune sur 1rst ou heures dépassées-------')
            print('prio est :',self.priority['priority'][self.idx_first_target])

            if first_target['set or rise'] == 'set':
                print('set')
                moon_idx_set_target+=1
                if ((moon_idx_set_target+self.nb_tot_telescopes)>=len(idx_set_targets_sorted)):
                    idx_safe_1rst+=1
                    self.idx_first_target=index_prio[-idx_safe_1rst]
                    first_target=self.priority[self.idx_first_target]
                    is_moon_constraint_met_first_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_first_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                    print('safety solu',self.targets[self.idx_first_target],is_moon_constraint_met_first_target,hours_constraint_first)
                else:
                    self.idx_first_target=idx_set_targets_sorted[-(self.nb_tot_telescopes+moon_idx_set_target)]
                    if self.priority['priority'][self.idx_first_target]!=float('-inf'):
                        first_target=self.priority[idx_set_targets_sorted[-(self.nb_tot_telescopes+moon_idx_set_target)]]
                        is_moon_constraint_met_first_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_first_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                        print('before safety',self.targets[self.idx_first_target],is_moon_constraint_met_first_target,hours_constraint_first)
                    else:
                        is_moon_constraint_met_first_target=False

            if before_change_first_target['set or rise'] == 'rise':
                print('rise')
                moon_idx_rise_target+=1
                if ((moon_idx_rise_target+self.nb_tot_telescopes)>=len(idx_rise_targets_sorted)):
                    idx_safe_2nd+=1
                    self.idx_first_target=index_prio[-idx_safe_2nd]
                    first_target=self.priority[self.idx_first_target]
                    # print(self.targets[self.idx_first_target])
                    is_moon_constraint_met_first_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_first_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                    print('safety solu',self.targets[self.idx_first_target],is_moon_constraint_met_first_target)

                else:
                    self.idx_first_target=idx_rise_targets_sorted[-(self.nb_tot_telescopes+moon_idx_rise_target)]
                    if self.priority['priority'][self.idx_first_target]!=float('-inf'):
                        first_target=self.priority[idx_rise_targets_sorted[-(self.nb_tot_telescopes+moon_idx_rise_target)]]
                        is_moon_constraint_met_first_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_first_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                        print('before safety',is_moon_constraint_met_first_target,hours_constraint_first)
                    else:
                        is_moon_constraint_met_first_target=False


            if (before_change_first_target['set or rise'] != 'rise') and (before_change_first_target['set or rise'] != 'set'):
                print('other')
                print(first_target['set or rise'])
                idx_safe+=1
                self.idx_first_target=index_prio[-idx_safe]
                first_target=self.priority[self.idx_first_target]
                is_moon_constraint_met_first_target = is_observable(constraints_all,self.observatory, self.targets[self.idx_first_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t)))
                print('safety solu',self.targets[self.idx_first_target],is_moon_constraint_met_first_target)
            hours_constraint_first=self.nb_hours_ok(self.idx_first_target,self.nb_hours_observed)
            print('end of loop constraints',)

        variable=self.update_hours_observed_first(t)
        self.nb_hours_observed=variable[0]
        self.nb_hours=variable[1]
        print(self.targets[self.idx_first_target],'nb_hours_first',self.nb_hours[self.idx_first_target])
        if is_moon_constraint_met_first_target:
            is_moon_constraint_met_first_target=True
        else:
            is_moon_constraint_met_first_target=False
        return is_moon_constraint_met_first_target,hours_constraint_first,self.idx_first_target


    def is_constraints_met_second_target(self,t,idx_set_targets_sorted,idx_rise_targets_sorted,index_prio):
        """
            Useful when the moon constrain (< 30°) or the hours constraint (nb_hours_observed>nb_hours_threshold)
            are not fullfilled anymore and the second target needs to be changed

        Parameters
        ----------
            t: int, days in the month
            idx_first_target: int, index associated to the first target (with regards to the targets list)
            idx_second_target: int, index associated to the second target (with regards to the targets list)
            idx_set_targets_sorted: list of int, index of all the set target ranked y priority
            idx_rise_targets_sorted: list of int, index of all the rise target ranked y priority
            index_prio: list of int, index of all the rise target ranked y priority

        Returns
        -------
            is_moon_constraint_met_second_target: boolean, say if moon constraint is fullfilled on second target
            hours_constraint: boolean, say if hour constraint is fullfilled on second target
            idx_second_target: int, index for the second target (if constraints where already fullfilled the Value
            is the same as input, else it is a new value for a new target)

        """
        dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')
        first_target=self.priority[self.idx_first_target]
        second_target=self.priority[self.idx_second_target]
        #dt_nautical_civil_evening_2=(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next')-self.observatory.twilight_evening_civil(self.time_ranges[self.month_obs][0]+dt*t,which='next'))/2
        #dt_civil_nautical_morning_2=(self.observatory.twilight_morning_civil(self.time_ranges[self.month_obs][0]+dt*(t+1),which='next')-self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='next'))/2
        night_duration=((Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest') ))-(Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))))# print(night_duration.iso)
        moon_idx_set_target=0
        moon_idx_rise_target=0
        print('sunset/sunrise in 2nd target :',Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next')).iso,Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest')).iso)

        constraints_all = [AltitudeConstraint(min=24*u.deg), MoonSeparationConstraint(min=30*u.deg), \
        TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))), \
        (Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest') )))]
        is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
        hours_constraint=self.nb_hours_ok(self.idx_second_target,self.nb_hours_observed)
        idx_safe=1
        idx_safe_1rst=1
        idx_safe_2nd=1
        print('2nd target is:',second_target)
        print('constraints 2nd:',is_moon_constraint_met_second_target, hours_constraint)
        if self.idx_first_target==self.idx_second_target:
            print(' 2nd=1ere ')
        else:
            while not (is_moon_constraint_met_second_target & hours_constraint):

                before_change_second_target=self.priority[self.idx_second_target]
                before_change_idx_second_target=self.idx_second_target
                print('------lune sur 2nd ou heures dépassées------')
                print('prio est :',self.priority['priority'][self.idx_second_target])

                if before_change_second_target['set or rise'] == 'set':
                    print('set')
                    moon_idx_set_target+=1
                    if ((moon_idx_set_target+self.nb_tot_telescopes)>=len(idx_set_targets_sorted)):
                        idx_safe_1rst+=1
                        self.idx_second_target=index_prio[-idx_safe_1rst]
                        second_target=self.priority[self.idx_second_target]
                        is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                        print('safety solu if 1rst=set',self.targets[self.idx_second_target],is_moon_constraint_met_second_target,hours_constraint)
                        print('before safety if 1rst=set prio est :',self.priority['priority'][self.idx_second_target])
                    else:
                        self.idx_second_target=idx_set_targets_sorted[-(self.nb_tot_telescopes+moon_idx_set_target)]
                        if self.priority['priority'][self.idx_second_target]!=float('-inf'):
                            second_target=self.priority[idx_set_targets_sorted[-(self.nb_tot_telescopes+moon_idx_set_target)]]
                            is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                            print('before safety if 1rst=set',self.targets[self.idx_second_target],is_moon_constraint_met_second_target,hours_constraint)
                            print('before safety if 1rst=set prio est :',self.priority['priority'][self.idx_second_target])
                        else:
                            is_moon_constraint_met_second_target=False

                if before_change_second_target['set or rise'] == 'rise':
                    print('rise')
                    moon_idx_rise_target+=1
                    if ((moon_idx_rise_target+self.nb_tot_telescopes)>=len(idx_rise_targets_sorted)):
                        idx_safe_2nd+=1
                        self.idx_second_target=index_prio[-idx_safe_2nd]
                        second_target=self.priority[self.idx_second_target]
                        # print(self.targets[self.idx_second_target])
                        is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                        print('safety solu if 1rst=rise',self.targets[self.idx_second_target],is_moon_constraint_met_second_target)

                    else:
                        self.idx_second_target=idx_rise_targets_sorted[-(self.nb_tot_telescopes+moon_idx_rise_target)]
                        if self.priority['priority'][self.idx_second_target]!=float('-inf'):
                            second_target=self.priority[idx_rise_targets_sorted[-(self.nb_tot_telescopes+moon_idx_rise_target)]]
                            is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                            print('before safety if 1rst=rise prio est',self.targets[self.idx_second_target],is_moon_constraint_met_second_target)
                        else:
                            is_moon_constraint_met_second_target=False

                if (before_change_second_target['set or rise'] != 'rise') and (before_change_second_target['set or rise'] != 'set'):
                    print('other')
                    # print(second_target['set or rise'])
                    if first_target['set or rise'] == 'rise':
                        moon_idx_set_target+=1
                        print('moon_idx_set_target+self.nb_tot_telescopes',moon_idx_set_target+self.nb_tot_telescopes,'len(idx_set_targets_sorted)',len(idx_set_targets_sorted))
                        if ((moon_idx_set_target+self.nb_tot_telescopes)>=len(idx_set_targets_sorted)):
                            idx_safe_1rst+=1
                            self.idx_second_target=index_prio[-idx_safe_1rst]
                            second_target=self.priority[self.idx_second_target]
                            is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                            print('safety solu if 1rst=rise',self.targets[self.idx_second_target],is_moon_constraint_met_second_target)
                        else:
                            self.idx_second_target=idx_set_targets_sorted[-(self.nb_tot_telescopes+moon_idx_set_target)]
                            if self.priority['priority'][self.idx_second_target]!=float('-inf'):
                                second_target=self.priority[idx_set_targets_sorted[-(self.nb_tot_telescopes+moon_idx_set_target)]]
                                is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                                print('before safety if 1rst=rise',self.targets[self.idx_second_target],is_moon_constraint_met_second_target)
                            else:
                                is_moon_constraint_met_second_target=False

                    if first_target['set or rise'] == 'set':
                        moon_idx_rise_target+=1
                        if ((moon_idx_rise_target+self.nb_tot_telescopes)>=len(idx_rise_targets_sorted)):
                            idx_safe_2nd+=1
                            self.idx_second_target=index_prio[-idx_safe_2nd]
                            second_target=self.priority[self.idx_second_target]
                            # print(self.targets[self.idx_second_target])
                            is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                            print('safety solu if 1rst=set',self.targets[self.idx_second_target],is_moon_constraint_met_second_target)

                        else:
                            self.idx_second_target=idx_rise_targets_sorted[-(self.nb_tot_telescopes+moon_idx_rise_target)]
                            if self.priority['priority'][self.idx_second_target]!=float('-inf'):
                                second_target=self.priority[idx_rise_targets_sorted[-(self.nb_tot_telescopes+moon_idx_rise_target)]]
                                is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
                                print('before safety if 1rst=set',self.targets[self.idx_second_target],is_moon_constraint_met_second_target)
                            else:
                                is_moon_constraint_met_second_target=False

                    if first_target['set or rise'] == 'both':
                        #idx_safe+=1
                        self.idx_second_target=self.idx_first_target
                        second_target=self.idx_first_target
                        is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next')))
                        print('safety solu if 1rst=both, 2nd=1rst',self.targets[self.idx_second_target],self.priority['priority'][self.idx_second_target],is_moon_constraint_met_second_target)

                hours_constraint=self.nb_hours_ok(self.idx_second_target,self.nb_hours_observed)

        variable=self.update_hours_observed_second(t)
        self.nb_hours_observed=variable[0]
        self.nb_hours=variable[1]
        print(self.targets[self.idx_second_target],'nb_hours_second',self.nb_hours[self.idx_second_target])

        if is_moon_constraint_met_second_target:
            is_moon_constraint_met_second_target=True
        else:
            is_moon_constraint_met_second_target=False
        return is_moon_constraint_met_second_target,hours_constraint,self.idx_second_target

    # def combined_saintex_artemis():
    #

    def nb_hours_ok(self,idx_target,nb_hours_observed):
        """
            Check if number of hours is ok

        Parameters
        ----------
            idx_target: int, index of the target you want to check
            nb_hours_observed: list of all the number of hours oserved for each targets of target listnb_hours_observed
        Returns
        -------
            is_hours_constraint_met_target: boolean, say the hour constraint is ok or not

        """
        is_hours_constraint_met_target = True
        a=(1-nb_hours_observed[idx_target]/(self.nb_hours_threshold[idx_target]+20))
        if a < 1E-5:
            is_hours_constraint_met_target = False
        return is_hours_constraint_met_target

    def moon_ok_2nd(self):
        """
            Check if moon > 30° is ok on second target

        Parameters
        ----------
            idx_second_target: int, index of the second target
        Returns
        -------
            is_moon_constraint_met_target: boolean, say the moon constraint is ok or not on the second target

        """
        dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')
        #dt_nautical_civil_evening_2=(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next')-self.observatory.twilight_evening_civil(self.time_ranges[self.month_obs][0]+dt*t,which='next'))/2
        #dt_civil_nautical_morning_2=(self.observatory.twilight_morning_civil(self.time_ranges[self.month_obs][0]+dt*(t+1),which='next')-self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='next'))/2
        night_duration=((Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest')))-(Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))))# print(night_duration.iso)

        constraints_all = [AltitudeConstraint(min=24*u.deg), MoonSeparationConstraint(min=30*u.deg), \
        TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))), \
        (Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest') )))]

        is_moon_constraint_met_second_target = is_observable(constraints_all, self.observatory, self.targets[self.idx_second_target], times=Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + (night_duration/(2*u.day))*u.day))
        return is_moon_constraint_met_second_target

    def schedule_blocks(self,t):
        """
            schedule the lock thanks to astroplan tools

        Parameters
        ----------
            idx_first_target: int, index of the first target
            idx_second_target: int, index of the second target
        Returns
        -------
            SS1_night_blocks: astropy.table with name, start time, end time, duration ,
            coordinates (RE and DEC) and configuration (filter and exposure time) for each target of the night

        """
        first_target=self.priority[self.idx_first_target]
        second_target=self.priority[self.idx_second_target]
        dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')

        #dt_nautical_civil_evening_2=(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next')-self.observatory.twilight_evening_civil(self.time_ranges[self.month_obs][0]+dt*t,which='next'))/2
        #dt_civil_nautical_morning_2=(self.observatory.twilight_morning_civil(self.time_ranges[self.month_obs][0]+dt*(t+1),which='next')-self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='next'))/2
        night_duration=((Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest') ))-(Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))))# print(night_duration.iso)
        # print('sunset/sunrise in schedule blocks',Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest')).iso,Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next')).iso)
        dur_obs_both_target=(night_duration/(2*u.day))*2*u.day
        # dur_obs_set_target=(night_duration/(2*u.day))*u.day+1*(1*u.hour-t/15*u.hour)
        # dur_obs_rise_target=(night_duration/(2*u.day))*u.day-1*(1*u.hour-t/15*u.hour)
        aa=len_time_range=self.time_ranges[self.month_obs][1]-self.time_ranges[self.month_obs][0]
        dur_obs_set_target=(night_duration/(2*u.day))*u.day+1*((aa.value/30)*u.hour-t/(aa.value/2)*u.hour) #30 cause change of 2hour approximatively in one month (of course dep on the target)
        dur_obs_rise_target=(night_duration/(2*u.day))*u.day-1*((aa.value/30)*u.hour-t/(aa.value/2)*u.hour)

        #Constraints on first and second targets
        constraints_set_target=[AltitudeConstraint(min=24*u.deg), MoonSeparationConstraint(min=30*u.deg),
        TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))),(Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + dur_obs_set_target)))]
        print('INFO in SCHEDULE BLOCK 1rst target',Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next')).iso,(Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + dur_obs_set_target).iso))
        print('INFO in SCHEDULE BLOCK 2nd target',(Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + dur_obs_set_target).iso),Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest')).iso)

        constraints_rise_target=[AltitudeConstraint(min=24*u.deg), MoonSeparationConstraint(min=30*u.deg),
        TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next') + dur_obs_set_target)), (Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest') )))] #AtNightConstraint.twilight_astronomical()

        constraints_all = [AltitudeConstraint(min=24*u.deg), MoonSeparationConstraint(min=30*u.deg), \
        TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))), \
        (Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest') )))]

        blocks=[]
        if first_target['set or rise'] == 'set':
            print('case1')
            a=ObservingBlock(self.targets[self.idx_first_target],dur_obs_set_target,-1,constraints= constraints_set_target,configuration={'filt=' + str(self.filt_spc[self.idx_first_target]),'texp=' + str(self.tesxpspe_spc[self.idx_first_target])})
            blocks.append(a)
            b=ObservingBlock(self.targets[self.idx_second_target],dur_obs_rise_target,-1,constraints=constraints_rise_target,configuration={'filt=' + str(self.filt_spc[self.idx_second_target]),'texp=' + str(self.tesxpspe_spc[self.idx_second_target])})
            blocks.append(b)

        if first_target['set or rise'] == 'both':
            print('case2')
            a=ObservingBlock(self.targets[self.idx_first_target],dur_obs_both_target,-1,constraints= constraints_all,configuration={'filt=' + str(self.filt_spc[self.idx_first_target]),'texp=' + str(self.tesxpspe_spc[self.idx_first_target])})
            blocks.append(a)

        if first_target['set or rise'] == 'rise':
            print('case3')
            b=ObservingBlock(self.targets[self.idx_second_target],dur_obs_set_target,-1,constraints=constraints_set_target,configuration={'filt=' + str(self.filt_spc[self.idx_second_target]),'texp=' + str(self.tesxpspe_spc[self.idx_second_target])})
            blocks.append(b)
            a=ObservingBlock(self.targets[self.idx_first_target],dur_obs_rise_target,-1,constraints=constraints_rise_target,configuration={'filt=' + str(self.filt_spc[self.idx_first_target]),'texp=' + str(self.tesxpspe_spc[self.idx_first_target])})
            blocks.append(a)

        transitioner = Transitioner(slew_rate= 11*u.deg/u.second)
        seq_schedule_SS1=Schedule(self.time_ranges[self.month_obs][0]+dt*t,self.time_ranges[self.month_obs][0]+dt*(t+1))
        sequen_scheduler_SS1=SPECULOOSScheduler(constraints=constraints_all, observer=self.observatory,transitioner=transitioner)
        sequen_scheduler_SS1(blocks,seq_schedule_SS1)

        #make night blocks
        SS1_night_blocks=seq_schedule_SS1.to_table()
        name_all=SS1_night_blocks['target']
        name=[]
        for i,nam in enumerate(name_all):
            name.append(nam)
        return SS1_night_blocks,name

    def update_hours_observed_first(self,t):
        """
            update number of hours observed for the corresponding first target

        Parameters
        ----------
            t: int, day of the month
            idx_first_target: int, index of the first target
        Returns
        -------
            self.nb_hours_observed
            self.nb_hours
            self.nb_hours[idx_first_target]

        """
        dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')

        first_target=self.priority[self.idx_first_target]
        night_duration=((Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest') ))-(Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))))# print(night_duration.iso)

        dur_obs_both_target=(night_duration/(2*u.day))*2*u.day
        # dur_obs_set_target=(night_duration/(2*u.day))*u.day+1*(1*u.hour-t/15*u.hour)
        # dur_obs_rise_target=(night_duration/(2*u.day))*u.day-1*(1*u.hour-t/15*u.hour)
        aa=len_time_range=self.time_ranges[self.month_obs][1]-self.time_ranges[self.month_obs][0]
        dur_obs_set_target=(night_duration/(2*u.day))*u.day+1*((aa.value/30)*u.hour-t/(aa.value/2)*u.hour)
        dur_obs_rise_target=(night_duration/(2*u.day))*u.day-1*((aa.value/30)*u.hour-t/(aa.value/2)*u.hour)


        if first_target['set or rise'] == 'set':
            nb_hours__1rst_old=self.nb_hours_observed[self.idx_first_target]
            self.nb_hours_observed[self.idx_first_target]+=dur_obs_set_target.value*24 #in hours
            self.nb_hours[self.idx_first_target]=[nb_hours__1rst_old,self.nb_hours_observed[self.idx_first_target]]
        if first_target['set or rise'] == 'both':
            nb_hours__1rst_old=self.nb_hours_observed[self.idx_first_target]
            self.nb_hours_observed[self.idx_first_target]+=dur_obs_both_target.value*24 #in hours
            self.nb_hours[self.idx_first_target]=[nb_hours__1rst_old,self.nb_hours_observed[self.idx_first_target]]
        if first_target['set or rise'] == 'rise':
            nb_hours__1rst_old=self.nb_hours_observed[self.idx_first_target]
            self.nb_hours_observed[self.idx_first_target]+=dur_obs_rise_target.value*24 #in hours
            self.nb_hours[self.idx_first_target]=[nb_hours__1rst_old,self.nb_hours_observed[self.idx_first_target]]

        return self.nb_hours_observed,self.nb_hours,self.nb_hours[self.idx_first_target]

    def update_hours_observed_second(self,t):
        """
            update number of hours observed for the corresponding second target

        Parameters
        ----------
            t: int, day of the month
            idx_second_target: int, index of the second target
        Returns
        -------
            self.nb_hours_observed
            self.nb_hours
            self.nb_hours[idx_second_target]

        """
        dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')
        second_target=self.priority[self.idx_second_target]
        dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')
        # dt_nautical_civil_evening_2=(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next')-self.observatory.twilight_evening_civil(self.time_ranges[self.month_obs][0]+dt*t,which='next'))/2
        # dt_civil_nautical_morning_2=(self.observatory.twilight_morning_civil(self.time_ranges[self.month_obs][0]+dt*(t+1),which='next')-self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='next'))/2
        night_duration=((Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt*(t+1),which='nearest') ))-(Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0]+dt*t,which='next'))))# print(night_duration.iso)

        dur_obs_both_target=(night_duration/(2*u.day))*2*u.day
        # dur_obs_set_target=(night_duration/(2*u.day))*u.day+1*(1*u.hour-t/15*u.hour)
        # dur_obs_rise_target=(night_duration/(2*u.day))*u.day-1*(1*u.hour-t/15*u.hour)
        aa=len_time_range=self.time_ranges[self.month_obs][1]-self.time_ranges[self.month_obs][0]
        dur_obs_set_target=(night_duration/(2*u.day))*u.day+1*((aa.value/30)*u.hour-t/(aa.value/2)*u.hour)
        dur_obs_rise_target=(night_duration/(2*u.day))*u.day-1*((aa.value/30)*u.hour-t/(aa.value/2)*u.hour)


        if second_target['set or rise'] == 'rise':
            nb_hours__2nd_old=self.nb_hours_observed[self.idx_second_target]
            self.nb_hours_observed[self.idx_second_target]+=dur_obs_rise_target.value*24 #in hours
            self.nb_hours[self.idx_second_target]=[nb_hours__2nd_old,self.nb_hours_observed[self.idx_second_target]]
        if second_target['set or rise'] == 'set':
            nb_hours__2nd_old=self.nb_hours_observed[self.idx_second_target]
            self.nb_hours_observed[self.idx_second_target]+=dur_obs_set_target.value*24 #in hours
            self.nb_hours[self.idx_second_target]=[nb_hours__2nd_old,self.nb_hours_observed[self.idx_second_target]]

        return self.nb_hours_observed,self.nb_hours,self.nb_hours[self.idx_second_target]

    def update_month_plan(self,SS1_night_blocks,file_txt):

        '''Modify Month_plan but also return
        the colmuns "name" with the name of
        all the target that where observed
        in the night block'''

        SS1_night_blocks.add_index('target')
        try:
            index_to_delete=SS1_night_blocks.loc['TransitionBlock'].index
            SS1_night_blocks.remove_row(index_to_delete)
        except KeyError:
            print('no transition block')
        name_all=SS1_night_blocks['target']
        Start_all=SS1_night_blocks['start time (UTC)']
        Finish_all=SS1_night_blocks['end time (UTC)']
        Duration_all=SS1_night_blocks['duration (minutes)']
        with open(file_txt, 'a') as file_plan:
            for i,nam in enumerate(name_all):
                if i == 1:
                    if nam==name_all[i-1]:
                        Start_all[i]=Start_all[i-1]
                        Duration_all[i]=(Time(Finish_all[i])-Time(Start_all[i])).value*24*60
                        SS1_night_blocks.remove_row(0)

                idx_target=np.where((nam==self.priority['target name']))[0]
                if nam!='TransitionBlock':
                    file_plan.write(nam + ' ' + '\"' + Start_all[i] + '\"' + ' ' + '\"' + Finish_all[i] + '\"' + ' ' + self.tel + ' '  + \
                    str(np.round(self.nb_hours[idx_target[0]][0],3)) + '/' + str(self.nb_hours_threshold[idx_target[0]]) + ' ' + \
                    str(np.round(self.nb_hours[idx_target[0]][1],3)) + '/' + str(self.nb_hours_threshold[idx_target[0]]) + '\n')

    def make_night_block(self,t,SS1_night_blocks):

        dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')

        SS1_night_blocks.add_index('target')
        try:
            index_to_delete=SS1_night_blocks.loc['TransitionBlock'].index
            SS1_night_blocks.remove_row(index_to_delete)
        except KeyError:
            print('no transition block')
        print('avant astuce',SS1_night_blocks)
        panda_table=SS1_night_blocks.to_pandas()
        panda_table.to_csv(os.path.join('night_blocks_' + self.tel + '_' +  Time(self.time_ranges[self.month_obs][0]+dt*t).tt.datetime.strftime("%Y-%m-%d ") + '.txt'),sep=' ',index=False)

    def update_priority(self,name):
        for nam in name:
            print('nam',nam)
            scheduled_targets_idx=(self.priority['target name']==nam)
            self.priority['priority'][scheduled_targets_idx]=float('-inf')
        return self.priority


def run_all(a,idx_set_targets_sorted,idx_rise_targets_sorted,index_prio):

    SS1_night_blocks_all=[]
    df=[]
    df2=[]
    name=[]
    month_obs=int(a.month_obs)
    time_ranges=a.time_ranges
    delta_t=Time(time_ranges[month_obs][1])-Time(time_ranges[month_obs][0])
    for x in range(0,int(delta_t.sec/60/60/24)):
        print('le jour est :',x)
        # b=a.idx_targets()
        df.append(a.is_constraints_met_first_target(x,idx_set_targets_sorted,idx_rise_targets_sorted,index_prio))
        # idx_first_target=df[x][2] #A AJIOUTER OU SUPPRIMER ?
        df2.append(a.is_constraints_met_second_target(x,idx_set_targets_sorted,idx_rise_targets_sorted,index_prio))
        # self.idx_second_target=df2[x][2] #A AJOUTER OU SUPPRIMER ?
        SS1_night_blocks_all.append(a.schedule_blocks(x)[0])
        print('targets idx are',a.idx_first_target,a.idx_second_target)
        a.update_month_plan(SS1_night_blocks_all[x],'month_plan_test.txt')
        a.make_night_block(x,SS1_night_blocks_all[x])
        name.append(a.schedule_blocks(x)[1])
        # a.priority=a.update_priority(name)
    NAMES=np.concatenate(name)
    a.priority=a.update_priority(NAMES)

    print('name targets scheduled',name)
    return df,df2,SS1_night_blocks_all,a.priority,name





class make_schedule:

    def __init__(self):
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
        self.target_list = None
        self.observatory = None #observatory
        self.date_range = None #date_range
        self.strategy = None
        self.duration_segments = None #duration_segments
        self.nb_segments =  None #nb_segments
        self.constraints = None
        self.priority = None
        self.priority_by_day = []
        self.index_prio = None
        self.index_prio_by_day = []
        self.priority_ranked = None
        self.priority_ranked_by_day = []
        self.reverse_df1 = None
        self.targets = None
        self.idx_first_target = None
        self.idx_first_target_by_day = []
        self.idx_second_target = None
        self.idx_second_target_by_day = []
        self.first_target = None
        self.first_target_by_day = []
        self.second_target = None
        self.second_target_by_day = []
        self.moon_and_visibility_contraint_table = None
        self.time_ranges = [Time(['2019-01-01 12:00:00', '2019-01-31 12:00:00']),Time(['2019-02-01 12:00:00', '2019-02-28 12:00:00']),Time(['2019-03-01 15:00:00', '2019-03-31 15:00:00']),Time(['2019-04-01 15:00:00', '2019-04-30 15:00:00']),Time(['2019-05-01 15:00:00', '2019-05-31 15:00:00']),Time(['2019-06-01 15:00:00', '2019-06-30 15:00:00']),Time(['2019-07-01 12:00:00', '2019-07-31 12:00:00']),Time(['2019-08-01 12:00:00', '2019-08-31 12:00:00']),Time(['2019-09-01 12:00:00', '2019-09-30 12:00:00']),Time(['2019-10-01 12:00:00', '2019-10-31 12:00:00']),Time(['2019-11-01 12:00:00', '2019-11-30 12:00:00']),Time(['2019-12-01 12:00:00', '2019-12-31 12:00:00'])]

    @property
    def idx_rise_targets_sorted(self):
        idx_rise_targets=(self.priority_ranked['set or rise']=='rise')
        idx_rise_targets_sorted=self.index_prio[idx_rise_targets]
        return idx_rise_targets_sorted
    @property
    def idx_set_targets_sorted(self):
        idx_set_targets=(self.priority_ranked['set or rise']=='set')
        idx_set_targets_sorted=self.index_prio[idx_set_targets]
        return idx_set_targets_sorted

    @property
    def target_table_spc(self):
         target_table_spc = Table.read(self.target_list, format='ascii')
         return target_table_spc

    @property
    def months_obs(self):
        date_format = "%Y-%m-%d %H:%M:%S.%f"
        for i,t in enumerate(self.time_ranges):
            if (datetime.strptime(t[0].value,date_format) <= datetime.strptime(self.date_range[0].value,date_format) <= datetime.strptime(t[1].value,date_format)) and\
                (datetime.strptime(t[0].value,date_format) <= datetime.strptime(self.date_range[1].value,date_format) <= datetime.strptime(t[1].value,date_format)):
                return  i
            if (datetime.strptime(t[0].value,date_format) <= datetime.strptime(self.date_range[0].value,date_format) <= datetime.strptime(t[1].value,date_format)) and\
                (datetime.strptime(t[1].value,date_format) <= datetime.strptime(self.date_range[1].value,date_format)) :
                if i < (len(self.time_ranges) - 1):
                    return i+1
                if i == (len(self.time_ranges) - 1):
                    return i

    @property
    def date_range_in_days(self):
        date_format = "%Y-%m-%d %H:%M:%S.%f"
        date_start = datetime.strptime(self.date_range[0].value, date_format)
        date_end = datetime.strptime(self.date_range[1].value, date_format)
        date_range_in_days = (date_end - date_start).days
        return date_range_in_days

    @property
    def nb_hours_threshold(self):
        nb_hours_threshold = [50]* len(self.target_table_spc)
        return nb_hours_threshold

    @property
    def date_ranges_day_by_day(self):
        date_format = "%Y-%m-%d %H:%M:%S.%f"
        d2 = datetime.strptime(self.date_range[1].value, date_format)
        i = 1
        t = datetime.strptime(self.date_range[0].value, date_format)
        u = t
        date_ranges_day_by_day = []
        while t < d2:
            d = timedelta(days=1)
            t = u + d * i
            i += 1
            date_ranges_day_by_day.append(Time(t))
        return date_ranges_day_by_day


    def load_parameters(self,input_file):
        with open(input_file, "r") as f:
            Inputs = yaml.load(f)
            self.target_list = Inputs['target_list']
            self.observatory = charge_observatories(Inputs['observatories'][0])[0]
            self.date_range = Time(Inputs['date_range'])#,Time(Inputs['date_range'][1])]
            self.strategy = Inputs['strategy']
            self.duration_segments = Inputs['duration_segments']
            self.nb_segments = Inputs['nb_segments']
            self.constraints = [AtNightConstraint.twilight_astronomical()]


    def schedule(self, Altitude_constraint = None, Moon_constraint = None):
        import time
        start = time.time()
        if Altitude_constraint:
            self.constraints.append(AltitudeConstraint(min=float(Altitude_constraint)*u.deg))
        if Moon_constraint:
            self.constraints.append(MoonSeparationConstraint(min=float(Moon_constraint)*u.deg))

        self.targets = target_list_good_coord_format(self.target_list)

        if str(self.strategy) == 'alternative':
            self.is_moon_and_visibility_contraint()
            print()

        if str(self.strategy) == 'continuous':

            self.reverse_df1=reverse_Observability(self.observatory,self.targets,self.constraints,self.time_ranges)
            end = time.time()
            print('reverse_df1',end - start)


            for t in range(1,self.date_range_in_days):

                product = self.table_priority_prio(self.date_ranges_day_by_day[t])
                end = time.time()
                print('index_prio', end - start)
                self.index_prio = product[0]
                self.priority = product[1]
                end = time.time()
                print('priority', end - start)
                self.idx_targets()
                end = time.time()
                print('idx_targets', end - start)

                self.moon_and_visibility_contraint_table = self.is_moon_and_visibility_contraint(self.date_ranges_day_by_day[t])

                print('day number:',t)

                self.priority_by_day.append(self.priority)
                self.index_prio_by_day.append(self.index_prio)
                self.priority_ranked_by_day.append(self.priority_ranked)

                if self.is_constraints_met_first_target(t):
                    self.first_target = self.priority[self.idx_first_target]

                    self.first_target_by_day.append(self.first_target)
                    self.idx_first_target_by_day.append(self.idx_first_target)


                if not self.is_constraints_met_first_target(t):
                    print('WARNING: impossible to find the first target that respects the constraints')

                if self.idx_second_target is None:
                    print('WARNING: no second target')

                if self.idx_second_target != None and self.is_constraints_met_second_target(t):

                    self.second_target = self.priority[self.idx_second_target]
                    self.second_target_by_day.append(self.second_target )
                    self.idx_second_target_by_day.append(self.idx_second_target)



        if str(self.strategy) == 'segmented':
            print()

    def idx_targets(self):
        """
            Give the index of the first and second targets as well as the
            corresponding row in the priority table

        Parameters
        ----------

        Returns
        -------
            idx_first_target: index first target in target list
            first_target: row number idx_first_target in priority table
            idx_second_target: index second target in target list
            second_target: row number idx_second_target in priority table

        """

        self.idx_first_target=self.index_prio[-1]
        self.first_target=self.priority[self.idx_first_target]
        dt_1day=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day
        for i in range(2,len(self.index_prio)):
            if self.first_target['set or rise'] == 'rise':
                if self.priority['set or rise'][self.index_prio[-i]]=='set':
                    if (self.observatory.target_rise_time(self.date_range[0],self.targets[self.index_prio[-i]],which='nearest',horizon=24*u.deg) \
                        < self.observatory.twilight_evening_nautical(self.date_range[0],which='nearest')) and \
                        (self.observatory.target_set_time(self.date_range[0],self.targets[self.index_prio[-i]],which='next',horizon=24*u.deg) \
                        > self.observatory.target_rise_time(self.date_range[0],self.targets[self.idx_first_target],which='nearest',horizon=24*u.deg)):

                        self.idx_second_target=self.index_prio[-i]
                        self.second_target=self.priority[self.idx_second_target]
                        break
                    else:
                        self.second_target = None
                        self.idx_second_target = None

                if self.priority['set or rise'][self.index_prio[-i]]=='both':
                    self.second_target = None
                    self.idx_second_target = None
                    print('still searching for the second target that is no both')

            if self.first_target['set or rise'] == 'set':
                if self.priority['set or rise'][self.index_prio[-i]]=='rise':
                    print('ici','index prio',self.index_prio[-i])
                    if (self.observatory.target_set_time(self.date_range[0],self.targets[self.index_prio[-i]],which='next',horizon=24*u.deg) \
                        > self.observatory.twilight_morning_nautical(self.date_range[0]+dt_1day,which='nearest')) and \
                        (self.observatory.target_rise_time(self.date_range[0],self.targets[self.index_prio[-i]],which='nearest',horizon=24*u.deg) \
                        < self.observatory.target_set_time(self.date_range[0],self.targets[self.idx_first_target],which='nearest',horizon=24*u.deg)):
                        print('tu as passé l\'étape 1')
                        self.idx_second_target=self.index_prio[-i]
                        self.second_target=self.priority[self.idx_second_target]
                        break
                    else:
                        self.second_target = None
                        self.idx_second_target = None

                if self.priority['set or rise'][self.index_prio[-i]]=='both':
                    self.second_target = None
                    self.idx_second_target = None
                    print('still searching for the second target that is no both')

            if self.first_target['set or rise'] == 'both':
                self.second_target = None
                self.idx_second_target = None
                break

            if self.idx_second_target is None:
                print('no second target available')

        #print('IN idx_targets','idx_first_target',self.idx_first_target,'idx_second_target',self.idx_second_target)
        #return self.first_target,self.second_target
        # idx_set_targets=(self.priority_ranked['set or rise']=='set')
        # idx_rise_targets=(self.priority_ranked['set or rise']=='set')
        # self.idx_set_targets_sorted=self.index_prio[idx_set_targets]
        # self.idx_rise_targets_sorted=self.index_prio[idx_rise_targets]

    def table_priority_prio(self,day):

        self.targets = target_list_good_coord_format(self.target_list)
        self.priority=Table(names=('priority','target name','set or rise','alt set start','alt rise start','alt set end','alt rise end'),dtype=('f4','S11','S4','f4','f4','f4','f4'))
        self.Ranking_month(day)

        #self.priority = self.priority[0]
        set_targets_index = (self.priority['alt set start']>20) & (self.priority['alt set end']>20)
        self.priority['set or rise'][set_targets_index] = 'set'
        rise_targets_index = (self.priority['alt rise start']>20) & (self.priority['alt rise end']>20)
        self.priority['set or rise'][rise_targets_index]='rise'
        both_targets_index=(self.priority['alt rise start']>20) & (self.priority['alt set start']>20) & (self.priority['alt rise end']>20) & (self.priority['alt set end']>20)
        self.priority['set or rise'][both_targets_index] = 'both'
        self.priority['priority'][both_targets_index] = self.priority['priority'][both_targets_index]*10
        priority_non_observable_idx = (self.priority['set or rise']=='None')
        self.priority['priority'][priority_non_observable_idx] = 0.5
        self.index_prio = np.argsort(self.priority['priority'])
        self.priority_ranked = self.priority[self.index_prio]

        return self.index_prio, self.priority, self.priority_ranked

    def Ranking_month(self,day):
        """
            Rank all targets with respect to their altitude (whether they are up or not),
            also if their best month score is in agreement with the desired month of observation,
            if their Kmag is inferior to 10.5, if the priospe (equivalent to JWST SNR) is high,
            and if the numer of hours oserved is lower than the threshold (50 hours for the present strategy)

        Parameters
        ----------
            priority: astopy table, gives the observation priority, whether the target is
            more a set target (up a sunset for the whole month) a rise target (up at sunrise for the whole month)
            or both, the altitude at sunset at the beginning of the month,
            the altitude at sunrise at the beginning of the month, the altitude at sunset at the end of the month,
            the altitude at sunrise at the end of the month
            i: int, id target
            date_range_in_days : nb of days in the date range
            observatory: the observatory
            month obs: int, the number of the month of observation, among 0 to 11 (January to December)
            months: the best score of each target for each month
            months_2nd_option: the second best score of each target for each month
            months_3rd_option: the third best score of each target for each month
            months_4th_option: the fourth best score of each target for each month
            months_5th_option: the fifth best score of each target for each month
            time_ranges: list of astropy.Time range with start and end times
            targets : target list on the FixedTarget() format from astroplan
            name_spc: all the name of the targets in the target list
            priospe_spc: SNR JWST * 100, for the targets in the target list
            reverse_df1: inverse of observaility table
            nb_hours_observed: number of hours observed for each target of the target list
            nb_hours_threshold: number of hours to complete for each target of the target list
            Kmag_spc: Kmagnitude of each target of the target list

        Returns
        -------
            priority: the same priority as given in function arguments with updated first columns
            with values:
            -0.5 for non observable targets, <0 for targets where the numer of hours observed > nb hours nb_hours_threshold
            >0 for the others, the highest value belonging to the highest priority target

        """

        # if Sptype_spc[i]=='M6':
        #     factor_spctype=1
        # else:
        #     factor_spctype=2
        b = []
        nb_days = 1
        dt_1day = Time('2018-01-02 00:00:00', scale='tcg') - Time('2018-01-01 00:00:00', scale='tcg')
        day_fmt = Time(day.iso , out_subfmt = 'date').iso

        if os.path.exists('Ranking_months_' + str(self.observatory.name) + '_' + str(day_fmt) + '_' + str(len(self.targets)) + '.csv'):
            name_file = 'Ranking_months_' + str(self.observatory.name) + '_' + str(day_fmt) + '_' + str(len(self.targets)) + '.csv'
            dataframe_ranking_months = pd.read_csv(name_file, delimiter=',')
            self.priority = Table.from_pandas(dataframe_ranking_months)

        else:
            for i in range(len(self.target_table_spc['Name'])):

                print(self.target_table_spc['Name'][i])
                target = self.targets[i]
                b.append(month_option(self.target_table_spc['Name'][i], self.reverse_df1))
                month_option_var = vstack(Table(b))
                months = month_option_var['months']
                print(months)
                months_2nd_option = month_option_var['months_2nd_option']
                months_3rd_option = month_option_var['months_3rd_option']
                months_4th_option = month_option_var['months_4th_option']
                months_5th_option = month_option_var['months_5th_option']

                if self.target_table_spc['prio_spc'][i] >= 100:
                    factor_priospe_spc = 4
                else:
                    factor_priospe_spc = 1

                if (int(self.target_table_spc['nb_hours_surved'][i]) > 0) and (int(self.target_table_spc['nb_hours_surved'][i]) < int(self.nb_hours_threshold[i])):
                    factor_on_going = 10 ** (1 + 1 / (self.nb_hours_threshold[i] - self.target_table_spc['nb_hours_surved'][i]))

                if (int(self.target_table_spc['nb_hours_surved'][i]) == 0):
                    factor_on_going = 10 ** (1 / (self.nb_hours_threshold[i] - self.target_table_spc['nb_hours_surved'][i]))

                if (int(self.target_table_spc['nb_hours_surved'][i]) >= int(self.nb_hours_threshold[i])):
                    factor_on_going = -1

                if months[i] == 'Month' + str(self.months_obs):
                    # times = _generate_24hr_grid(day, 0, date_range_in_days, date_range_in_days)
                    times = _generate_24hr_grid(day, 0, nb_days, 1)
                    self.priority.add_row(
                        ((self.target_table_spc['prio_spc'][i]) ** factor_priospe_spc * (2 ** (10 * (self.reverse_df1[self.target_table_spc['Name'][i]][self.months_obs])))
                         * max(self.observatory.altaz(times, target, grid_times_targets=True).alt.value) * factor_on_going,
                         self.target_table_spc['Name'][i], 'None',
                         self.observatory.altaz(self.observatory.twilight_evening_nautical(day, which='nearest'),
                                           target).alt.value,
                         self.observatory.altaz(self.observatory.twilight_morning_nautical(day, which='next'),
                                           target).alt.value,
                         self.observatory.altaz(self.observatory.twilight_evening_nautical(day + nb_days * dt_1day, which='nearest'),
                                           target).alt.value,
                         self.observatory.altaz(self.observatory.twilight_morning_nautical(day + nb_days * dt_1day, which='next'),
                                           target).alt.value))

                else:
                    if months_2nd_option[i] == str(self.months_obs) or months_3rd_option[i] ==  str(self.months_obs) or \
                            months_4th_option[i] == str(self.months_obs) or months_5th_option[i] ==  str(
                            self.months_obs):
                        # times = _generate_24hr_grid(date_range[0], 0, date_range_in_days,date_range_in_days)
                        times = _generate_24hr_grid(day, 0, 1, 1)
                        self.priority.add_row(
                            ((self.target_table_spc['prio_spc'][i]) ** factor_priospe_spc * (2 ** (10 * (self.reverse_df1[self.target_table_spc['Name'][i]][self.months_obs])))
                             * max(self.observatory.altaz(times, target, grid_times_targets=True).alt.value) * factor_on_going,
                             self.target_table_spc['Name'][i], 'None',
                             self.observatory.altaz(self.observatory.twilight_evening_nautical(day, which='nearest'),
                                               target).alt.value,
                             self.observatory.altaz(self.observatory.twilight_morning_nautical(day, which='next'),
                                               target).alt.value,
                             self.observatory.altaz(self.observatory.twilight_evening_nautical(day + nb_days * dt_1day, which='nearest'),
                                               target).alt.value,
                             self.observatory.altaz(self.observatory.twilight_morning_nautical(day + nb_days * dt_1day, which='next'),
                                               target).alt.value))

                    else:
                        self.priority.add_row((-0.5, self.target_table_spc['Name'][i], 'None',
                                          self.observatory.altaz(self.observatory.twilight_evening_nautical(day, which='nearest'),
                                                            target).alt.value,
                                          self.observatory.altaz(self.observatory.twilight_morning_nautical(day, which='next'),
                                                            target).alt.value,
                                          self.observatory.altaz(self.observatory.twilight_evening_nautical(day + nb_days * dt_1day, which='nearest'),
                                                            target).alt.value,
                                          self.observatory.altaz(self.observatory.twilight_morning_nautical(day + nb_days * dt_1day, which='next'),
                                                            target).alt.value))

            dataframe_priority = self.priority.to_pandas()
            dataframe_priority.to_csv('Ranking_months_' + str(self.observatory.name) + '_' + str(day_fmt) + '_' + str(len(self.targets)) + '.csv', sep= ',',index=False)

    def run_all(a,idx_set_targets_sorted,idx_rise_targets_sorted,index_prio):

        SS1_night_blocks_all=[]
        df=[]
        df2=[]
        name=[]
        delta_t=Time(self.date_range[1])-Time(self.date_range[month_obs][0])
        for x in range(0,int(delta_t.sec/60/60/24)):
            print('le jour est :',x)
            # b=a.idx_targets()
            df.append(a.is_constraints_met_first_target(x,idx_set_targets_sorted,idx_rise_targets_sorted,index_prio))
            # idx_first_target=df[x][2] #A AJIOUTER OU SUPPRIMER ?
            df2.append(a.is_constraints_met_second_target(x,idx_set_targets_sorted,idx_rise_targets_sorted,index_prio))
            # self.idx_second_target=df2[x][2] #A AJOUTER OU SUPPRIMER ?
            SS1_night_blocks_all.append(a.schedule_blocks(x)[0])
            print('targets idx are',a.idx_first_target,a.idx_second_target)
            a.update_month_plan(SS1_night_blocks_all[x],'month_plan_test.txt')
            a.make_night_block(x,SS1_night_blocks_all[x])
            name.append(a.schedule_blocks(x)[1])
            # a.priority=a.update_priority(name)
        NAMES=np.concatenate(name)
        a.priority=a.update_priority(NAMES)

        print('name targets scheduled',name)
        return df,df2,SS1_night_blocks_all,a.priority,name

    def shift_hours_observation(self,idx_target):

        date_format = "%Y-%m-%d %H:%M:%S.%f"

        if self.priority['set or rise'][idx_target]=='set':

            set_time_begin = datetime.strptime(self.observatory.target_set_time(self.date_range[0],self.targets[idx_target],which='next',horizon=24*u.deg).iso, date_format)
            set_time_end = datetime.strptime(self.observatory.target_set_time(self.date_range[1],self.targets[idx_target],which='next',horizon=24*u.deg).iso, date_format)
            shift_hours_observation = (set_time_end - set_time_begin).hours

        if self.priority['set or rise'][idx_target]=='rise':

            rise_time_begin = datetime.strptime(self.observatory.target_rise_time(self.date_range[0],self.targets[idx_target],which='next',horizon=24*u.deg).iso, date_format)
            rise_time_end = datetime.strptime(self.observatory.target_rise_time(self.date_range[1],self.targets[idx_target],which='next',horizon=24*u.deg).iso, date_format)
            shift_hours_observation = (rise_time_end - rise_time_begin).hours

    def schedule_blocks(self,t):
        """
            schedule the lock thanks to astroplan tools

        Parameters
        ----------
            idx_first_target: int, index of the first target
            idx_second_target: int, index of the second target
        Returns
        -------
            SS1_night_blocks: astropy.table with name, start time, end time, duration ,
            coordinates (RE and DEC) and configuration (filter and exposure time) for each target of the night

        """

        first_target=self.priority[self.idx_first_target]
        second_target=self.priority[self.idx_second_target]
        dt_1day=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')

        night_duration=((Time(self.observatory.twilight_morning_nautical(self.date_range[0]+dt_1day*(t+1),which='nearest') ))\
            -(Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next'))))

        dur_obs_both_target=(night_duration/(2*u.day))*2*u.day

        dur_obs_set_target=(night_duration/(2*u.day))*u.day+1*((aa.value/30)*u.hour-t/(aa.value/2)*u.hour) #30 cause change of 2hour approximatively in one month (of course dep on the target)

        dur_obs_rise_target=(night_duration/(2*u.day))*u.day-1*((aa.value/30)*u.hour-t/(aa.value/2)*u.hour)

        #Constraints on first and second targets
        constraints_set_target=[AltitudeConstraint(min=24*u.deg), MoonSeparationConstraint(min=30*u.deg),
        TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next'))),(Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next') + dur_obs_set_target)))]
        print('INFO in SCHEDULE BLOCK 1rst target',Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next')).iso,(Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next') + dur_obs_set_target).iso))
        print('INFO in SCHEDULE BLOCK 2nd target',(Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next') + dur_obs_set_target).iso),Time(self.observatory.twilight_morning_nautical(self.date_range[0]+dt_1day*(t+1),which='nearest')).iso)

        constraints_rise_target=[AltitudeConstraint(min=24*u.deg), MoonSeparationConstraint(min=30*u.deg),
        TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next') + dur_obs_set_target)), (Time(self.observatory.twilight_morning_nautical(self.date_range[0]+dt_1day*(t+1),which='nearest') )))] #AtNightConstraint.twilight_astronomical()

        constraints_all = [AltitudeConstraint(min=24*u.deg), MoonSeparationConstraint(min=30*u.deg), \
        TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next'))), \
        (Time(self.observatory.twilight_morning_nautical(self.date_range[0]+dt_1day*(t+1),which='nearest') )))]

        blocks=[]
        if first_target['set or rise'] == 'set':
            print('case1')
            a=ObservingBlock(self.targets[self.idx_first_target],dur_obs_set_target,-1,constraints= constraints_set_target,configuration={'filt=' + str(self.filt_spc[self.idx_first_target]),'texp=' + str(self.tesxpspe_spc[self.idx_first_target])})
            blocks.append(a)
            b=ObservingBlock(self.targets[self.idx_second_target],dur_obs_rise_target,-1,constraints=constraints_rise_target,configuration={'filt=' + str(self.filt_spc[self.idx_second_target]),'texp=' + str(self.tesxpspe_spc[self.idx_second_target])})
            blocks.append(b)

        if first_target['set or rise'] == 'both':
            print('case2')
            a=ObservingBlock(self.targets[self.idx_first_target],dur_obs_both_target,-1,constraints= constraints_all,configuration={'filt=' + str(self.filt_spc[self.idx_first_target]),'texp=' + str(self.tesxpspe_spc[self.idx_first_target])})
            blocks.append(a)

        if first_target['set or rise'] == 'rise':
            print('case3')
            b=ObservingBlock(self.targets[self.idx_second_target],dur_obs_set_target,-1,constraints=constraints_set_target,configuration={'filt=' + str(self.filt_spc[self.idx_second_target]),'texp=' + str(self.tesxpspe_spc[self.idx_second_target])})
            blocks.append(b)
            a=ObservingBlock(self.targets[self.idx_first_target],dur_obs_rise_target,-1,constraints=constraints_rise_target,configuration={'filt=' + str(self.filt_spc[self.idx_first_target]),'texp=' + str(self.tesxpspe_spc[self.idx_first_target])})
            blocks.append(a)

        transitioner = Transitioner(slew_rate= 11*u.deg/u.second)
        seq_schedule_SS1=Schedule(self.date_range[0]+dt_1day*t,self.date_range[0]+dt_1day*(t+1))
        sequen_scheduler_SS1=SPECULOOSScheduler(constraints=constraints_all, observer=self.observatory,transitioner=transitioner)
        sequen_scheduler_SS1(blocks,seq_schedule_SS1)

        #make night blocks
        SS1_night_blocks=seq_schedule_SS1.to_table()
        name_all=SS1_night_blocks['target']
        name=[]
        for i,nam in enumerate(name_all):
            name.append(nam)
        return SS1_night_blocks,name

    def is_constraint_hours(self,idx_target):
        """
            Check if number of hours is ok

        Parameters
        ----------
            idx_target: int, index of the target you want to check

        Returns
        -------
            is_hours_constraint_met_target: boolean, say the hour constraint is ok or not

        """
        nb_hours_observed = self.target_table_spc['nb_hours_surved']
        #self.nb_hours_threshold = [50] * len(self.target_table_spc)
        is_hours_constraint_met_target = True
        a=(1-nb_hours_observed[idx_target]/(self.nb_hours_threshold[idx_target]+20))
        if a < 1E-5:
            is_hours_constraint_met_target = False
        return is_hours_constraint_met_target

    def __is_moon_and_visibility_contraint(self,t,idx_target):
        """
            Check if number of moon not too close and target is visible

        Parameters
        ----------
            idx_target: int, index of the target you want to check
        Returns
        -------
            is_hours_constraint_met_target: boolean, say the hour constraint is ok or not

        """
        if isinstance(idx_target, int):
            idx_target = [idx_target]

        dt_1day=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')
        night_duration = ((Time(self.observatory.twilight_morning_nautical(self.date_range[0]+dt_1day*(t+1),which='nearest')))\
            -(Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next'))))# print(night_duration.iso)

        is_moon_constraint_met = is_observable(self.constraints, self.observatory, self.targets, \
            times= Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next') + night_duration/2))
            #Time([self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next').iso,self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next').iso]))
            #Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next') + (night_duration/(2*u.day))*u.day))

        return is_moon_constraint_met

    def night_duration(self, day):
        '''

        :param day: day str format '%y%m%d HH:MM:SS.sss'
        :return:
        '''
        dt_1day=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')
        return ((Time(self.observatory.twilight_morning_nautical(day + dt_1day ,which='nearest')))\
            -(Time(self.observatory.twilight_evening_nautical(day ,which='next'))))

    def info_obs_possible(self,day):
        '''

        :param day: day str format '%y%m%d HH:MM:SS.sss'
        :return:
        '''
        duration_obs_possible = self._is_moon_and_visibility_contraint(day)['fraction of time observable']#np.subtract(np.asarray(self.set_time_targets(day).value) , np.asarray(self.rise_time_targets(day).value))
        #duration_obs_possible = duration_obs_possible*24
        duration_obs_possible[duration_obs_possible <0]= 0
        return duration_obs_possible

    def rise_time_targets(self,day):
        '''

        :param day: day str format '%y%m%d HH:MM:SS.sss'
        :return:
        '''
        rise_time = self.observatory.target_rise_time(day, self.targets, which='next', horizon=24 * u.deg)
        return rise_time

    def set_time_targets(self,day):
        '''

        :param day: day str format '%y%m%d HH:MM:SS.sss'
        :return:
        '''
        dt_1day = Time('2018-01-02 00:00:00', scale='tcg') - Time('2018-01-01 00:00:00', scale='tcg')
        set_time = self.observatory.target_set_time(day+1*dt_1day,self.targets, which='nearest', horizon=24 * u.deg)
        return set_time

    def idx_is_julien_criterion(self,day):
        '''

        :param day: day str format '%y%m%d HH:MM:SS.sss'
        :return:
        '''
        nb_hour_wanted = 3.5 #hours
        idx_is_julien_criterion = (self.info_obs_possible(day) > nb_hour_wanted )
        return idx_is_julien_criterion



    def is_moon_and_visibility_contraint(self,day):
        """
            Check if number of moon not too close and target is visible

        Parameters
        ----------
            idx_target: int, index of the target you want to check
        Returns
        -------
            is_hours_constraint_met_target: boolean, say the hour constraint is ok or not

        """


        dt_1day=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')
        nd = Time(self.night_duration(day))

        observability_table_day = observability_table(self.constraints, self.observatory, self.targets , \
            time_range = Time([Time(self.observatory.twilight_evening_nautical(Time(day),which='next')).iso,Time(self.observatory.twilight_morning_nautical(Time(day+dt_1day),which='nearest')).iso]) ) # + nd/2

        observability_table_day['fraction of time observable'] =  observability_table_day['fraction of time observable'] * \
                                                                  (self.observatory.twilight_morning_nautical(day + dt_1day,which='nearest') - self.observatory.twilight_evening_nautical( day, which='nearest')) * 24

        return observability_table_day

    def is_constraints_met_first_target(self,t):
        """
            Useful when the moon constrain (< 30°) or the hours constraint (nb_hours_observed>nb_hours_threshold)
            are not fullfilled anymore and the first target needs to be changed

        Parameters
        ----------
            t: int, days in the month
            idx_first_target: int, index associated to the first target (with regards to the targets list)
            idx_set_targets_sorted: list of int, index of all the set target ranked y priority
            idx_rise_targets_sorted: list of int, index of all the rise target ranked y priority
            index_prio: list of int, index of all the rise target ranked y priority

        Returns
        -------
            is_moon_constraint_met_first_target: boolean, say if moon constraint is fullfilled on first target
            hours_constraint: boolean, say if hour constraint is fullfilled on first target
            idx_first_target: int, index for the first target (if constraints where already fullfilled the Value
            is the same as input, else it is a new value for a new target)

        """
        dt_1day=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')
        night_duration = ((Time(self.observatory.twilight_morning_nautical(self.date_range[0]+dt_1day*(t+1),which='nearest')))\
            -(Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next'))))# print(night_duration.iso)


        # self.constraints.append(TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next'))), \
        #     (Time(self.observatory.twilight_morning_nautical(self.date_range[0]+dt_1day*(t+1),which='nearest') ))) )


        #see if moon constraint is met
        is_moon_constraint_met_first_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_first_target]
        hours_constraint_first = self.is_constraint_hours(self.idx_first_target)

        idx_safe=1
        idx_safe_1rst_set=1
        idx_safe_1rst_rise=1
        moon_idx_set_target=0
        moon_idx_rise_target=0

        while not (is_moon_constraint_met_first_target & hours_constraint_first):

            before_change_first_target=self.priority[self.idx_first_target]
            before_change_idx_first_target=self.idx_first_target

            if before_change_first_target['set or rise'] == 'set':
                moon_idx_set_target+=1
                if ( moon_idx_set_target >= len(self.idx_set_targets_sorted)):
                    idx_safe_1rst_set+=1
                    self.idx_first_target=self.index_prio[-idx_safe_1rst_set]
                    self.first_target=self.priority[self.idx_first_target]
                    is_moon_constraint_met_first_target = self.moon_and_visibility_contraint_table ['ever observable'][self.idx_first_target]
                    hours_constraint_first=self.is_constraint_hours(self.idx_first_target)
                else:
                    self.idx_first_target=self.idx_set_targets_sorted[-moon_idx_set_target]
                    self.first_target=self.priority[self.idx_first_target]
                    if self.priority['priority'][self.idx_first_target]!=float('-inf'):
                        is_moon_constraint_met_first_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_first_target]
                        hours_constraint_first=self.is_constraint_hours(self.idx_first_target)
                    else:
                        is_moon_constraint_met_first_target = False
                        hours_constraint_first = False

            if before_change_first_target['set or rise'] == 'rise':
                moon_idx_rise_target+=1
                if ((moon_idx_rise_target) >= len(self.idx_rise_targets_sorted)):
                    idx_safe_1rst_rise+=1
                    self.idx_first_target=self.index_prio[-idx_safe_1rst_rise]
                    self.first_target=self.priority[self.idx_first_target]
                    is_moon_constraint_met_first_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_first_target]
                    hours_constraint_first=self.is_constraint_hours(self.idx_first_target)
                else:
                    self.idx_first_target=self.idx_rise_targets_sorted[-(moon_idx_rise_target)]
                    self.first_target=self.priority[self.idx_first_target]
                    if self.priority['priority'][self.idx_first_target]!=float('-inf'):
                        is_moon_constraint_met_first_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_first_target]
                        hours_constraint_first=self.is_constraint_hours(self.idx_first_target)
                    else:
                        is_moon_constraint_met_first_target=False
                        hours_constraint_first = False


            if (before_change_first_target['set or rise'] != 'rise') and (before_change_first_target['set or rise'] != 'set'):
                idx_safe+=1
                self.idx_first_target=self.index_prio[-idx_safe]
                self.first_target=self.priority[self.idx_first_target]
                is_moon_constraint_met_first_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_first_target]
                hours_constraint_first=self.is_constraint_hours(self.idx_first_target)

        # variable=self.update_hours_observed_first(t)
        # self.nb_hours_observed=variable[0]
        # self.nb_hours=variable[1]

        if is_moon_constraint_met_first_target and hours_constraint_first:
            is_constraints_met_first_target=True
        else:
            is_constraints_met_first_target=False
        return is_constraints_met_first_target

    def is_constraints_met_second_target(self,t):
        """
            Useful when the moon constrain (< 30°) or the hours constraint (nb_hours_observed>nb_hours_threshold)
            are not fullfilled anymore and the second target needs to be changed

        Parameters
        ----------
            t: int, days in the month
            idx_first_target: int, index associated to the first target (with regards to the targets list)
            idx_second_target: int, index associated to the second target (with regards to the targets list)
            idx_set_targets_sorted: list of int, index of all the set target ranked y priority
            idx_rise_targets_sorted: list of int, index of all the rise target ranked y priority
            index_prio: list of int, index of all the rise target ranked y priority

        Returns
        -------
            is_moon_constraint_met_second_target: boolean, say if moon constraint is fullfilled on second target
            hours_constraint: boolean, say if hour constraint is fullfilled on second target
            idx_second_target: int, index for the second target (if constraints where already fullfilled the Value
            is the same as input, else it is a new value for a new target)

        """
        dt_1day=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg')
        night_duration = ((Time(self.observatory.twilight_morning_nautical(self.date_range[0]+dt_1day*(t+1),which='nearest')))\
            -(Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next'))))# print(night_duration.iso)


        # first_target=self.priority[self.idx_first_target]
        # second_target=self.priority[self.idx_second_target]


        # self.constraints.append(TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.date_range[0]+dt_1day*t,which='next'))), \
        #     (Time(self.observatory.twilight_morning_nautical(self.date_range[0]+dt_1day*(t+1),which='nearest') ))) )

        is_moon_constraint_met_second_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_second_target]
        hours_constraint_second=self.is_constraint_hours(self.idx_second_target)

        idx_safe=1
        idx_safe_2nd_set=1
        idx_safe_2nd_rise=1
        moon_idx_set_target=0
        moon_idx_rise_target=0

        if self.idx_first_target==self.idx_second_target:
            print(' 2nd=1ere ')
        if self.idx_second_target == None:
            print(' No second target for that night')
        else:
            while not (is_moon_constraint_met_second_target & hours_constraint_second):

                before_change_second_target=self.priority[self.idx_second_target]
                before_change_idx_second_target=self.idx_second_target

                if before_change_second_target['set or rise'] == 'set':
                    moon_idx_set_target += 1
                    if ((moon_idx_set_target) >= len(self.idx_set_targets_sorted)):
                        idx_safe_2nd_set += 1
                        self.idx_second_target = self.index_prio[-idx_safe_2nd_set]
                        self.second_target = self.priority[self.idx_second_target]
                        is_moon_constraint_met_second_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_second_target]
                        hours_constraint_second = self.is_constraint_hours(self.idx_second_target)
                    else:
                        self.idx_second_target=self.idx_set_targets_sorted[-(moon_idx_set_target)]
                        self.second_target=self.priority[self.idx_second_target]
                        if self.priority['priority'][self.idx_second_target]!=float('-inf'):
                            is_moon_constraint_met_second_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_second_target]
                            hours_constraint_second=self.is_constraint_hours(self.idx_second_target)
                        else:
                            is_moon_constraint_met_second_target=False
                            hours_constraint_second = False

                if before_change_second_target['set or rise'] == 'rise':
                    moon_idx_rise_target+=1
                    if ((moon_idx_rise_target) >=len (self.idx_rise_targets_sorted)):
                        idx_safe_2nd_rise+=1
                        self.idx_second_target = self.index_prio[-idx_safe_2nd_rise]
                        self.second_target = self.priority[self.idx_second_target]
                        is_moon_constraint_met_second_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_second_target]
                        hours_constraint_second = self.is_constraint_hours(self.idx_second_target)
                    else:
                        self.idx_second_target=idx_rise_targets_sorted[-(moon_idx_rise_target)]
                        self.second_target=self.priority[self.idx_second_target]
                        if self.priority['priority'][self.idx_second_target]!=float('-inf'):
                            is_moon_constraint_met_second_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_second_target]
                            hours_constraint_second=self.is_constraint_hours(self.idx_second_target)
                        else:
                            is_moon_constraint_met_second_target=False
                            hours_constraint_second = False

                if (before_change_second_target['set or rise'] != 'rise') and (before_change_second_target['set or rise'] != 'set'):
                    if first_target['set or rise'] == 'rise':
                        moon_idx_set_target+=1
                        if ((moon_idx_set_target) >= len(self.idx_set_targets_sorted)):
                            idx_safe_2nd_set += 1
                            self.idx_second_target = self.index_prio[-idx_safe_2nd_set]
                            self.second_target = self.priority[self.idx_second_target]
                            is_moon_constraint_met_second_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_second_target]
                            hours_constraint_second=self.is_constraint_hours(self.idx_second_target)
                        else:
                            self.idx_second_target=idx_set_targets_sorted[-(moon_idx_set_target)]
                            self.second_target=self.priority[self.idx_second_target]
                            if self.priority['priority'][self.idx_second_target]!=float('-inf'):
                                is_moon_constraint_met_second_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_second_target]
                                hours_constraint_second=self.is_constraint_hours(self.idx_second_target)
                            else:
                                is_moon_constraint_met_second_target=False
                                hours_constraint_second = False

                    if first_target['set or rise'] == 'set':
                        moon_idx_rise_target+=1
                        if ((moon_idx_rise_target) >= len(self.idx_rise_targets_sorted)):
                            idx_safe_2nd_rise += 1
                            self.idx_second_target = self.index_prio[-idx_safe_2nd_rise]
                            self.second_target=self.priority[self.idx_second_target]
                            is_moon_constraint_met_second_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_second_target]
                            hours_constraint_second=self.is_constraint_hours(self.idx_second_target)

                        else:
                            self.idx_second_target=idx_rise_targets_sorted[-(moon_idx_rise_target)]
                            self.second_target=self.priority[self.idx_second_target]
                            if self.priority['priority'][self.idx_second_target]!=float('-inf'):
                                is_moon_constraint_met_second_target = self.moon_and_visibility_contraint_table['ever observable'][self.idx_second_target]
                                hours_constraint_second=self.is_constraint_hours(self.idx_second_target)
                            else:
                                is_moon_constraint_met_second_target=False
                                hours_constraint_second = False

                    if first_target['set or rise'] == 'both':
                        self.idx_second_target = self.idx_first_target
                        self.second_target = self.first_target
                        is_moon_constraint_met_second_target = self.is_moon_and_visibility_contraint(t,self.idx_second_target)

                # hours_constraint=self.nb_hours_ok(self.idx_second_target,self.nb_hours_observed)

        # variable=self.update_hours_observed_second(t)
        # self.nb_hours_observed=variable[0]
        # self.nb_hours=variable[1]
        if is_moon_constraint_met_second_target and hours_constraint_second:
            is_constraints_met_second_target=True
        else:
            is_constraints_met_second_target=False
        return is_constraints_met_second_target

    def reference_table(self):

        a = np.zeros((len(self.date_ranges_day_by_day), len(self.target_table_spc['Name']))).astype("object")
        b = np.zeros((len(self.date_ranges_day_by_day), len(self.target_table_spc['Name']))).astype("object")
        c = np.zeros((len(self.date_ranges_day_by_day), len(self.target_table_spc['Name']))).astype("object")
        d = np.zeros((len(self.date_ranges_day_by_day), len(self.target_table_spc['Name']))).astype("object")
        idx_true = [0] * len(self.date_ranges_day_by_day)
        is_julien_criterion = [False] * len(self.date_ranges_day_by_day)

        for i, day in enumerate(self.date_ranges_day_by_day):
            a[i] = self._is_moon_and_visibility_contraint(day)['ever observable']
            b[i] = self._is_moon_and_visibility_contraint(day)['fraction of time observable']
            idx_true = len(np.where(a[i])[0])
            is_julien_criterion[i] = np.sum(self.idx_is_julien_criterion(day))
            c[i] = self.rise_time_targets(day)
            d[i] = self.set_time_targets(day)

        E = np.zeros((len(self.date_ranges_day_by_day), len(self.target_table_spc['Name']), 5)).astype("object")
        E[:, :, 0] = self.target_table_spc['Name']
        E[:, :, 1] = a
        E[:, :, 2] = b
        E[:, :, 3] = c
        E[:, :, 4] = d

        np.save('array' + str(self.observatory.name) + '.npy', E, allow_pickle=True)

        df = pd.DataFrame({'day': Time(self.date_ranges_day_by_day), 'nb_observable_target': idx_true,
                           'julien_criterion': is_julien_criterion})

        df.to_csv('nb_observable_target_prio_50_' + str(self.observatory.name) + '.csv', sep=',', index=False)

