#code for long term scheduling for SPECULOOS
import numpy as np
import matplotlib
from astropy.table import Table, Column
from astropy.table import join
from astropy.table import unique
from astropy import units as u
from astroplan import FixedTarget, AltitudeConstraint, MoonSeparationConstraint,AtNightConstraint,observability_table, is_observable, months_observable,time_grid_from_range,LocalTimeConstraint, is_always_observable
from eScheduler.spe_schedule import SPECULOOSScheduler, PriorityScheduler, Schedule, ObservingBlock, Transitioner, DateElsa
from astroplan import TimeConstraint
from astropy.time import Time
import pandas as pd
import functools
from functools import reduce
import os

# definition of the different fonctions

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

    table_observability=observability_table(constraints, observatory, targets_observable, time_range=time_range) #calcul un edeuxieme fois is observable
    table_observability.remove_column('always observable')
    table_observability.remove_column('ever observable')
    table_observability.rename_column('fraction of time observable', 'Month' + str(j))

    return table_observability

def reverse_Observability(observatory,time_ranges,targets,constraints):
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
    tables_observability=list(map((lambda x: Observability(x,time_ranges[x],observatory,targets,constraints)), range(0,len(time_ranges))))
    df=list(map((lambda x: pd.DataFrame(tables_observability[x].to_pandas())),range(0,len(time_ranges))))
    a=reduce(functools.partial(pd.merge,how='outer', on='target name'), df)
    df=a.replace(to_replace=float('NaN'), value=0.0)
    df1 = df.set_index('target name')
    reverse_df1=df1.T
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



def Ranking_month(priority,i,observatory,month_obs,months,months_2nd_option,months_3rd_option,months_4th_option,months_5th_option,time_ranges,targets,name_spc,priospe_spc,reverse_df1,nb_hours_observed,nb_hours_threshold,Kmag_spc,Sptype_spc):
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
        i: int, the day of the month
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
#CHANGE HERE
    if priospe_spc[i]>=100:
        factor_priospe_spc=4
    else:
        factor_priospe_spc=1

    if (int(nb_hours_observed[i])>0) and (int(nb_hours_observed[i])<int(nb_hours_threshold[i])):
        factor_on_going=10**(1+1/(nb_hours_threshold[i]-nb_hours_observed[i]))
    if (int(nb_hours_observed[i])==0):
        factor_on_going=10**(1/(nb_hours_threshold[i]-nb_hours_observed[i]))
    # if (nb_hours_observed[i]==nb_hours_threshold[i]):
    #     factor_on_going=2*10**(1/(nb_hours_observed[i]-nb_hours_threshold[i]))
    if (int(nb_hours_observed[i])>=int(nb_hours_threshold[i])):
        factor_on_going=-1
    # print('sunset/sunrise :',Time(observatory.twilight_evening_nautical(time_ranges[month_obs][0],which='nearest')).iso,Time(observatory.twilight_morning_nautical(time_ranges[month_obs][0],which='next')).iso)
    altaz_all_targets=[]
    if months[i]=='Month' + str(month_obs) :
        times = _generate_24hr_grid(time_ranges[month_obs][0], 0, 30, 30)
        priority.add_row(((priospe_spc[i])**factor_priospe_spc*(2**(10*(reverse_df1[name_spc[i]][month_obs])))
        *max(observatory.altaz(times, targets[i], grid_times_targets=True).alt.value)*factor_on_going,name_spc[i],'None',
        observatory.altaz(observatory.twilight_evening_nautical(time_ranges[month_obs][0],which='nearest'),targets[i]).alt.value,
        observatory.altaz(observatory.twilight_morning_nautical(time_ranges[month_obs][0],which='next'),targets[i]).alt.value,
        observatory.altaz(observatory.twilight_evening_nautical(time_ranges[month_obs][1],which='nearest'),targets[i]).alt.value,
        observatory.altaz(observatory.twilight_morning_nautical(time_ranges[month_obs][1],which='next'),targets[i]).alt.value))

    else:
        if months_2nd_option[i]=='Month' + str(month_obs) or months_3rd_option[i]=='Month' + str(month_obs) or months_4th_option[i]=='Month' + str(month_obs) or months_5th_option[i]=='Month' + str(month_obs) :
            times = _generate_24hr_grid(time_ranges[month_obs][0], 0, 30,30)
            priority.add_row(((priospe_spc[i])**factor_priospe_spc*(2**(10*(reverse_df1[name_spc[i]][month_obs])))
            *max(observatory.altaz(times, targets[i], grid_times_targets=True).alt.value)*factor_on_going,name_spc[i],'None',
            observatory.altaz(observatory.twilight_evening_nautical(time_ranges[month_obs][0],which='nearest'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_morning_nautical(time_ranges[month_obs][0],which='next'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_evening_nautical(time_ranges[month_obs][1],which='nearest'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_morning_nautical(time_ranges[month_obs][1],which='next'),targets[i]).alt.value))

        else:
            priority.add_row((-0.5,name_spc[i],'None',
            observatory.altaz(observatory.twilight_evening_nautical(time_ranges[month_obs][0],which='nearest'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_morning_nautical(time_ranges[month_obs][0],which='next'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_evening_nautical(time_ranges[month_obs][1],which='nearest'),targets[i]).alt.value,
            observatory.altaz(observatory.twilight_morning_nautical(time_ranges[month_obs][1],which='next'),targets[i]).alt.value))

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
        # print('sunset/sunrise in idx_targets()',Time(self.observatory.twilight_evening_nautical(self.time_ranges[self.month_obs][0],which='nearest')).iso,Time(self.observatory.twilight_morning_nautical(self.time_ranges[self.month_obs][0]+dt,which='nearest')).iso)
        index_prio=np.argsort(self.priority['priority'])
        idx=1
        first_target=self.priority[index_prio[-idx]]
        priority_prio=self.priority[index_prio]
        # targets_prio=self.targets[index_prio]
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
