#!/anaconda3/bin/python3.6

import os
import math
from datetime import datetime
import importlib
import pandas as pd
import numpy as np
from jdcal import gcal2jd, jd2gcal
import csv
import sys
import astropy
from astropy import units as u
import datetime
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_sun, AltAz, EarthLocation
from astropy.coordinates import ICRS, get_moon
from astroplan import Observer
import julian
import time
from astroplan.plots import plot_airmass,plot_sky
from astroplan import FixedTarget, AltitudeConstraint, MoonSeparationConstraint,AtNightConstraint,AirmassConstraint,observability_table, is_observable, months_observable,time_grid_from_range,LocalTimeConstraint, is_always_observable
from astroplan import TimeConstraint
#from astroplan.Observer import twilight_evening_astronomical
import matplotlib
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.table import join
from astropy.table import unique
import pylab
from astroplan.plots import plot_schedule_airmass,plot_airmass,plot_sky
from eScheduler.spe_schedule import SPECULOOSScheduler, PriorityScheduler, Schedule, ObservingBlock, Transitioner, DateElsa
from eScheduler.spe_schedule import SequentialScheduler, Transitioner, PriorityScheduler, SimpleScheduler
from astropy.time import TimeDelta
from progress.bar import Bar
from astroplan.periodic import EclipsingSystem
from astroplan.constraints import is_event_observable



dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

constraints = [AltitudeConstraint(min=24*u.deg), AtNightConstraint()] #MoonSeparationConstraint(min=30*u.deg)


def scheduled_table(t,telescope,SS1_night_blocks_exists):

    Path='/Users/elsaducrot/Documents/GitHub/Scheduler_global/Python'
    tel = telescope
    try:
        os.path.exists(os.path.join(Path,tel,'night_blocks_' + tel + '_' +  Time(t).tt.datetime.strftime("%Y-%m-%d ") + '.txt'))
        print(os.path.join(Path,tel,'night_blocks_' + tel + '_' +  Time(t).tt.datetime.strftime("%Y-%m-%d ") + '.txt'))
    except NameError:
        print('no input night_block for this day')
    if SS1_night_blocks_exists:
        SS1_night_blocks_old=SS1_night_blocks
        scheduled_table=SS1_night_blocks_old
    else:
        scheduled_table=Table.read(os.path.join(Path,tel,'night_blocks_' + tel + '_' +  Time(t).tt.datetime.strftime("%Y-%m-%d ") + '.txt'), format='ascii')
    return scheduled_table
    # print('scheduled_table ---->',scheduled_table)

    # tel='Europa'
    t_now=time.strftime("%Y-%m-%d 12:00:00")
    today = input('Today? (y/n)')
    if today=='y':
        t_now=time.strftime("%Y-%m-%d 12:00:00")
    else:
        t_now = input('Enter date (%Y-%m-%d 12:00:00) : ')
    #Date different than todays's date
    # t_now='2018-11-01 12:00:00.000'

    nb_jours=int(input('nb days to prepare:'))
    t0=Time(t_now)
    # print(t0)


def transit_follow_up(t_now,t,observatory,constraints,targets_transit,name_spc_transit,timings_transit,T0_err_transit,periods_transit,P_err_transit,durations_transit,W_err_transit,filt_spc_transit,tesxpspe_spc_transit):
########################### TRANSITS FOLLOW UP ##################################
    observing_time = Time(t_now)
    t=Time(t_now)
    nb_jours=1#input('nb jours for follow up:')
    observing_time_end =Time(t+dt*nb_jours)
    # print(observing_time.iso,observing_time_end.iso)
    constraints = [AltitudeConstraint(min=24*u.deg), AtNightConstraint()] #MoonSeparationConstraint(min=30*u.deg)

    for i,nam in enumerate(name_spc_transit):
        blocks=[]
        try:
            SS1_night_blocks
        except NameError:
            SS1_night_blocks_exists = False
        else:
            SS1_night_blocks_exists = True
            SS1_night_blocks_old=SS1_night_blocks
        epoch = Time(timings_transit[i], format='jd')
        period = periods_transit[i] * u.day
        duration = durations_transit[i] * u.day
        oot_time=1.*u.hour
        target_transit=EclipsingSystem(primary_eclipse_time=epoch, orbital_period=period, duration=duration,
                               name=str(nam))
        print('target_transit',target_transit.next_primary_eclipse_time(observing_time,n_eclipses=1))
        timing_to_obs_jd=Time(target_transit.next_primary_eclipse_time(observing_time,n_eclipses=1)).jd
        # print('timing_to_obs_jd',timing_to_obs_jd[0],(timing_to_obs_jd[0]-epoch.jd)/periods_transit[i]*P_err_transit[i])

        target = FixedTarget.from_name(str(nam))
        n_transits =1#int(((observing_time_end-observing_time)/periods[i]).value)  # This is the roughly number of transits per year
        try:
            ing_egr = target_transit.next_primary_ingress_egress_time(observing_time,
                                                                n_eclipses=n_transits)
            # print(ing_egr.iso)
        except ValueError:
            print('No transit of ',nam,' on the period chosen')

        observable = is_event_observable(constraints, observatory, target,times_ingress_egress=ing_egr)
        observable_half=is_event_observable(constraints, observatory, target,times_ingress_egress=ing_egr)
        if np.any(observable):
            start_transit=Time(ing_egr[0][0]-T0_err_transit[i]*u.day-W_err_transit[i]*u.day-oot_time-((timing_to_obs_jd[0]-epoch.jd)/periods_transit[i]*P_err_transit[i])*1*u.day)
            print('start_transit',start_transit.iso)
            end_transit=Time(ing_egr[0][1]+T0_err_transit[i]*u.day+W_err_transit[i]*u.day+oot_time+((timing_to_obs_jd[0]-epoch.jd)/periods_transit[i]*P_err_transit[i])*1*u.day)
            print('end_transit',end_transit.iso)
            dur_obs_transit_target=(end_transit-start_transit).value*1.*u.day
            constraints_transit_target=[AltitudeConstraint(min=24*u.deg),MoonSeparationConstraint(min=30*u.deg), \
            TimeConstraint(start_transit,end_transit)]
            idx_first_target=int(np.where((name_spc_transit==nam))[0])
            a=ObservingBlock(targets_transit[idx_first_target],dur_obs_transit_target,-1,constraints= constraints_transit_target,configuration={'filt=' + str(filt_spc_transit[idx_first_target]),'texp=' + str(tesxpspe_spc_transit[idx_first_target])})
            blocks.append(a)
            transitioner = Transitioner(slew_rate= 11*u.deg/u.second)
            seq_schedule_SS1=Schedule(t,t+dt)
            sequen_scheduler_SS1=SPECULOOSScheduler(constraints=constraints_transit_target, observer=observatory,transitioner=transitioner)
            sequen_scheduler_SS1(blocks,seq_schedule_SS1)
            print(start_transit.iso,end_transit.iso)
            if (len(seq_schedule_SS1.to_table()['target'])!=0):
                SS1_night_blocks=seq_schedule_SS1.to_table() #Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
                return SS1_night_blocks
        else:
            print('pas de transit ce jour')

def special_target(t,observatory,input_name,targets,name_spc,filt_spc,tesxpspe_spc):
############## Scheduling of the special target ##################################
    # dt_nautical_civil_evening_2=(observatory.twilight_evening_nautical(t,which='next')-observatory.twilight_evening_civil(t,which='next'))/2
    # dt_civil_nautical_morning_2=(observatory.twilight_morning_civil(t+dt,which='nearest')-observatory.twilight_morning_nautical(t+dt,which='nearest'))/2
    night_duration=((Time(observatory.twilight_morning_nautical(t+dt,which='nearest')))-(Time(observatory.twilight_evening_nautical(t,which='next'))))#observatory.twilight_morning_astronomical(time_ranges[month_obs][0]+dt*t,which='next')-observatory.twilight_evening_astronomical(time_ranges[month_obs][0]+dt*t,which='nearest') #night_duration
    dur_obs_both_target=(night_duration/(2*u.day))*2*u.day
    # print(dur_obs_both_target)
    constraints_special_target=[AltitudeConstraint(min=24*u.deg),MoonSeparationConstraint(min=30*u.deg), \
    TimeConstraint((Time(observatory.twilight_evening_nautical(t,which='next'))), \
    (Time(observatory.twilight_morning_nautical(t+dt,which='nearest'))))]
    # print(Time(observatory.twilight_evening_nautical(t,which='next')-dt_nautical_civil_evening_2).iso,Time(observatory.twilight_morning_civil(t+dt,which='nearest') - dt_civil_nautical_morning_2).iso)
    idx_first_target=int(np.where((name_spc==input_name))[0])
    # print(idx_first_target)
    # nb_hours_start_special=nb_hours_observed[idx_first_target]
    # nb_hours_observed[idx_first_target]+=dur_obs_both_target.value*24
    # nb_hours[idx_first_target]=[nb_hours_start_special,nb_hours_observed[idx_first_target]]
    # nb_hours_special_target=nb_hours[idx_first_target]
    idx_special_target=idx_first_target
    # print(targets[idx_first_target])
    blocks=[]
    a=ObservingBlock(targets[idx_first_target],dur_obs_both_target,-1,constraints= constraints_special_target,configuration={'filt=' + str(filt_spc[idx_first_target]),'texp=' + str(tesxpspe_spc[idx_first_target])})
    blocks.append(a)
    transitioner = Transitioner(slew_rate= 11*u.deg/u.second)
    # print(t.iso,(t+dt).iso)
    seq_schedule_SS1=Schedule(t,t+dt)
    sequen_scheduler_SS1=SPECULOOSScheduler(constraints=constraints_special_target, observer=observatory,transitioner=transitioner)
    sequen_scheduler_SS1(blocks,seq_schedule_SS1)
    SS1_night_blocks=seq_schedule_SS1.to_table() #Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
    return SS1_night_blocks

def nb_hours_special(t,observatory,input_name,name_spc,nb_hours_observed,nb_hours):
    idx_first_target=int(np.where((name_spc==input_name))[0])
    # dt_nautical_civil_evening_2=(observatory.twilight_evening_nautical(t,which='next')-observatory.twilight_evening_civil(t,which='next'))/2
    # dt_civil_nautical_morning_2=(observatory.twilight_morning_civil(t+dt,which='nearest')-observatory.twilight_morning_nautical(t+dt,which='nearest'))/2
    night_duration=((Time(observatory.twilight_morning_nautical(t+dt,which='nearest')))-(Time(observatory.twilight_evening_nautical(t,which='next'))))#observatory.twilight_morning_astronomical(time_ranges[month_obs][0]+dt*t,which='next')-observatory.twilight_evening_astronomical(time_ranges[month_obs][0]+dt*t,which='nearest') #night_duration
    dur_obs_both_target=(night_duration/(2*u.day))*2*u.day
    nb_hours_start_special=nb_hours_observed[idx_first_target]
    nb_hours_observed[idx_first_target]+=dur_obs_both_target.value*24
    nb_hours[idx_first_target]=[nb_hours_start_special,nb_hours_observed[idx_first_target]]
    nb_hours_special_target=nb_hours[idx_first_target]
    return nb_hours_special_target

def monitoring(t,observatory,input_name,name_spc_mon,targets_mon,filt_spc,tesxpspe_spc,time_monitoring,airmass_max,airmax_min):
    ############################# MONITORING ##############################################
    # dt_nautical_civil_evening_2=(observatory.twilight_evening_nautical(t,which='next')-observatory.twilight_evening_civil(t,which='next'))/2
    # dt_civil_nautical_morning_2=(observatory.twilight_morning_civil(t+dt,which='nearest')-observatory.twilight_morning_nautical(t+dt,which='nearest'))/2
    #Scheduling of the monitiring target
    idx_first_target=int(np.where((name_spc_mon==input_name))[0])
    # print('target moni',targets_mon[idx_first_target])
    dur_mon_target=float(time_monitoring[idx_first_target])*u.minute
    # print(dur_mon_target)
    constraints_monitoring_target=[AltitudeConstraint(min=24*u.deg),MoonSeparationConstraint(min=30*u.deg), AirmassConstraint(airmass_max, airmax_min), \
    TimeConstraint((Time(observatory.twilight_evening_nautical(t,which='next'))), \
    (Time(observatory.twilight_morning_nautical(t+dt,which='nearest'))))]
    blocks=[]
    a=ObservingBlock(targets_mon[idx_first_target],dur_mon_target,-1,constraints= constraints_monitoring_target,configuration={'filt=' + str(filt_spc[idx_first_target]),'texp=' + str(tesxpspe_spc[idx_first_target])})
    blocks.append(a)

    transitioner = Transitioner(slew_rate= 11*u.deg/u.second)
    seq_schedule_SS1=Schedule(t,t+dt)
    sequen_scheduler_SS1=SPECULOOSScheduler(constraints=constraints_monitoring_target, observer=observatory,transitioner=transitioner)
    sequen_scheduler_SS1(blocks,seq_schedule_SS1)
    SS1_night_blocks=seq_schedule_SS1.to_table() #Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
    return SS1_night_blocks

# #######################
def planification(t,observatory,SS1_night_blocks_exists,scheduled_table,SS1_night_blocks):
    # dt_nautical_civil_evening_2=(observatory.twilight_evening_nautical(t,which='next')-observatory.twilight_evening_civil(t,which='next'))/2
    # dt_civil_nautical_morning_2=(observatory.twilight_morning_civil(t+dt,which='nearest')-observatory.twilight_morning_nautical(t+dt,which='nearest'))/2

    name=scheduled_table['target']
    scheduled_date_start=scheduled_table['start time (UTC)']
    scheduled_date_end=scheduled_table['end time (UTC)']
    scheduled_duration=scheduled_table['duration (minutes)']
    scheduled_ra1=scheduled_table['ra (h)']
    scheduled_ra2=scheduled_table['ra (m)']
    scheduled_ra3=scheduled_table['ra (s)']
    scheduled_dec1=scheduled_table['dec (d)']
    scheduled_dec2=scheduled_table['dec (m)']
    scheduled_dec3=scheduled_table['dec (s)']
    scheduled_configuration=scheduled_table['configuration']
    pd.DataFrame(scheduled_table.to_pandas())


    name_special_targets=SS1_night_blocks['target']
    Start_special_targets=SS1_night_blocks['start time (UTC)']
    Finish_special_targets=SS1_night_blocks['end time (UTC)']
    duration_special_targets=SS1_night_blocks['duration (minutes)']
    special_targets_ra1=SS1_night_blocks['ra (h)']
    special_targets_ra2=SS1_night_blocks['ra (m)']
    special_targets_ra3=SS1_night_blocks['ra (s)']
    special_targets_dec1=SS1_night_blocks['dec (d)']
    special_targets_dec2=SS1_night_blocks['dec (m)']
    special_targets_dec3=SS1_night_blocks['dec (s)']
    special_targets_configuration=SS1_night_blocks['configuration']


    for i,nam in enumerate(name):
        # print(nam)
        end_before_cut=scheduled_date_end[i]
        start_before_cut=scheduled_date_start[i]
        if SS1_night_blocks_exists:
            print('Several transit this night')
            if (np.any(str(SS1_night_blocks['target'][0]).find('Trappist-1'))==0) and (np.any(str(SS1_night_blocks_old['target'][0]).find('Trappist-1'))==0):
                print('deux trappist1!')
                Start_special_targets[0]=min(SS1_night_blocks['start time (UTC)'][0],SS1_night_blocks_old['start time (UTC)'][0])
                Finish_special_targets[0]=max(SS1_night_blocks['end time (UTC)'][0],SS1_night_blocks_old['end time (UTC)'][0])
        if (Start_special_targets[0]<Time(observatory.twilight_evening_nautical(t,which='next'))):
            Start_special_targets[0]=Time(observatory.twilight_evening_nautical(t,which='next')).iso
        if (scheduled_date_start[i] <= Start_special_targets[0]) and (Start_special_targets[0] <= scheduled_date_end[i]):
            if (Finish_special_targets[0]<= scheduled_date_end[i]): #special target in the middle of a scheduled target block
                #first part of scheduled target, change end time and duration
                print('case1')
                end_before_cut=scheduled_date_end[i]
                scheduled_date_end[i]=Start_special_targets[0]
                scheduled_duration[i]=(Time(Start_special_targets[0])-Time(scheduled_date_start[i])).value*24*60
                #row added for special target
                scheduled_table.add_row((name_special_targets[0],Start_special_targets[0],Finish_special_targets[0],duration_special_targets[0],special_targets_ra1[0],special_targets_ra2[0], \
                special_targets_ra3[0],special_targets_dec1[0],special_targets_dec2[0],special_targets_dec3[0],special_targets_configuration[0]))
                #row added for second part of scheduled target, change start time and duration
                scheduled_table.add_row(('_2',Finish_special_targets[0],end_before_cut,(Time(end_before_cut)-Time(Finish_special_targets[0])).value*24*60,scheduled_ra1[i],scheduled_ra2[i],scheduled_ra3[i], \
                scheduled_dec1[i],scheduled_dec2[i],scheduled_dec3[i],scheduled_configuration[i]))
                # print('ici',scheduled_table)
            if (Finish_special_targets[0]>= end_before_cut):
                print('case2')
                #Change end time of scheduled target
                scheduled_date_end[i]=Start_special_targets[0]
                scheduled_duration[i]=(Time(Start_special_targets[0])-Time(scheduled_date_start[i])).value*24*60
                #row added for special target
                scheduled_table.add_row((name_special_targets[0],Start_special_targets[0],Finish_special_targets[0],duration_special_targets[0],special_targets_ra1[0],special_targets_ra2[0], \
                special_targets_ra3[0],special_targets_dec1[0],special_targets_dec2[0],special_targets_dec3[0],special_targets_configuration[0]))
        if (start_before_cut >= Start_special_targets[0]):
            if (Finish_special_targets[0] <= end_before_cut):
                print('case3')
                if (Finish_special_targets[0] <= start_before_cut):
                    print('no change')
                else:
                    #Change start time of scheduled target and duration
                    scheduled_date_start[i]=Finish_special_targets[0]
                    scheduled_duration[i]=(Time(scheduled_date_end[i])-Time(Finish_special_targets[0])).value*24*60
                    #delete and replace scheduled target because change of start time and duration
                    scheduled_table.add_index('target')
                    index_to_delete=scheduled_table.loc[nam].index
                    scheduled_table.remove_row(index_to_delete)
                    scheduled_table.add_row((str(nam),Finish_special_targets[0],end_before_cut,scheduled_duration[i],scheduled_ra1[i],scheduled_ra2[i],scheduled_ra3[i], \
                    scheduled_dec1[i],scheduled_dec2[i],scheduled_dec3[i],scheduled_configuration[i]))
                    #row added for special target
                    scheduled_table.add_row((name_special_targets[0],Start_special_targets[0],Finish_special_targets[0],duration_special_targets[0],special_targets_ra1[0],special_targets_ra2[0], \
                    special_targets_ra3[0],special_targets_dec1[0],special_targets_dec2[0],special_targets_dec3[0],special_targets_configuration[0]))
            if (Finish_special_targets[0] >= end_before_cut):
                # print('case4')
                #delete scheduled target because special target overlap entirely
                scheduled_table.add_index('target')
                index_to_delete=scheduled_table.loc[nam].index
                scheduled_table.remove_row(index_to_delete)
                #row added for special target
                scheduled_table.add_row((name_special_targets[0],Start_special_targets[0],Finish_special_targets[0],duration_special_targets[0],special_targets_ra1[0],special_targets_ra2[0], \
                special_targets_ra3[0],special_targets_dec1[0],special_targets_dec2[0],special_targets_dec3[0],special_targets_configuration[0]))

    # idx_Delete=np.where(scheduled_table['duration (minutes)']<5)[0]
    # print('IDX DELETE',idx_Delete)
    # if len(idx_Delete)!=0:
    #     scheduled_table.remove_row(int(idx_Delete))
    scheduled_table_unique=unique(scheduled_table, keys='target')
    index_prio=np.argsort(scheduled_table_unique['start time (UTC)'])
    scheduled_table_sorted=scheduled_table_unique[index_prio]

    return scheduled_table_sorted

def update_month_plan(t,tel,input_name,name_spc,nb_hours_special_target,nb_hours_threshold,plan_file,scheduled_table_sorted):
    ##### MODIFY MONTH PLAN #######################################################

    Path='/Users/elsaducrot/Documents/GitHub/Scheduler_global/Python/'
    idx_special_target=int(np.where((name_spc==input_name))[0])
    name_final_targets=scheduled_table_sorted['target']
    Start_final_targets=scheduled_table_sorted['start time (UTC)']
    Finish_final_targets=scheduled_table_sorted['end time (UTC)']
    duration_final_targets=scheduled_table_sorted['duration (minutes)']
    target_table_month_plan = Table.read(plan_file+'.txt', format='ascii')
    hours_start=target_table_month_plan['nb_hours_start']
    hours_finish=target_table_month_plan['nb_hours_end']
    start_times=target_table_month_plan['Start']
    old_scheduled_table=Table.read(os.path.join(Path,tel,'night_blocks_' + tel + '_' +  Time(t).tt.datetime.strftime("%Y-%m-%d ") + '.txt'), format='ascii')
    old_scheduled_name=old_scheduled_table['target']
    old_scheduled_date_start=old_scheduled_table['start time (UTC)']
    old_scheduled_date_finish=old_scheduled_table['end time (UTC)']
    try :
        hours_start_special=str(np.round(nb_hours_special_target[0],3)) + '/' + str(nb_hours_threshold[idx_special_target]) #+ '/' + str(nb_hours_threshold[idx_target[0]]) + ' ' + \
        hours_finish_special=str(np.round(nb_hours_special_target[1],3)) + '/' + str(nb_hours_threshold[idx_special_target]) + '\n'
    except NameError:
        print('not a special target')

    for i,old in enumerate(old_scheduled_date_start):
        # idx_to_delete=(start_times==old_scheduled_date_start[0])
        for j in range(0,len(name_final_targets)):
            if (duration_final_targets[int(j)]>20):
                if (name_final_targets[int(j)]==input_name):
                    target_table_month_plan.add_row([[str(name_final_targets[int(j)])], [str(Start_final_targets[int(j)])], [str(Finish_final_targets[int(j)])], [tel + '_s'],[hours_start_special],[hours_finish_special]])
                else:
                    target_table_month_plan.add_row([[str(name_final_targets[int(j)])], [str(Start_final_targets[int(j)])], [str(Finish_final_targets[int(j)])], [tel + '_s'],[str(hours_start[int(j)])],[str(hours_finish[int(j)])]])

        # target_table_month_plan_unique=unique(target_table_month_plan, keys='Start')
        index_prio=np.argsort(target_table_month_plan['Start'])
        target_table_month_plan_sorted=target_table_month_plan[index_prio]
        target_table_month_plan_sorted.add_index(['Name','Start'])
        # print(old_scheduled_date_start[int(i)],old_scheduled_name[int(i)])
        try:
            ind_delete=(target_table_month_plan_sorted['Start']==old_scheduled_date_start[int(i)]) & (target_table_month_plan_sorted['Name']==old_scheduled_name[int(i)]) & (target_table_month_plan_sorted['Finish']==old_scheduled_date_finish[int(i)])
            index_to_delete=np.where(ind_delete)[0]
            target_table_month_plan_sorted.remove_row(int(index_to_delete))
        except TypeError:
            print('this block has been modified, no need to delete')
        # print(ind_delete,target_table_month_plan_sorted[ind_delete])
        # # print(target_table_month_plan_sorted.loc['2018-11-15 23:45:20.025'])
        # # index_to_delete=target_table_month_plan_sorted.loc(old_scheduled_date_start[int(i)])#.index
        # index_to_delete=np.where(ind_delete)[0]
        # print(np.where(ind_delete)[0])
        # target_table_month_plan_sorted.remove_row(int(index_to_delete))

    index_prio_tel=np.argsort(target_table_month_plan_sorted['Telescope'])
    target_table_month_plan_sorted=unique(target_table_month_plan_sorted[index_prio_tel],keys=['Name','Start','Finish','Telescope'])
    panda_table=target_table_month_plan_sorted.to_pandas()
    panda_table.to_csv(str(plan_file) + '.txt',sep=' ',index=False)
    return pd.DataFrame(panda_table)

    ###############################################################################
