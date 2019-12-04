#!/anaconda3/bin/python3.6

import os
import pandas as pd
import numpy as np
from SPOCK.upload_night_plans import upload_np_calli, upload_np_gany, upload_np_io, upload_np_euro,upload_np_artemis
from SPOCK.make_night_plans import make_np
import subprocess
import sys
from astropy import units as u
from astropy.time import Time
from astroplan.plots import dark_style_sheet, plot_airmass
from astropy.coordinates import SkyCoord, get_sun, AltAz, EarthLocation
from astroplan import Observer
import yaml
import shutil
from astroplan import FixedTarget, AltitudeConstraint, MoonSeparationConstraint,AtNightConstraint,AirmassConstraint,observability_table, is_observable, months_observable,time_grid_from_range,LocalTimeConstraint, is_always_observable
from astroplan import TimeConstraint
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.table import unique
from astroplan.plots import plot_schedule_airmass,plot_airmass,plot_sky
from eScheduler.spe_schedule import SPECULOOSScheduler, PriorityScheduler, Schedule, ObservingBlock, Transitioner, DateElsa
from eScheduler.spe_schedule import SequentialScheduler, Transitioner, PriorityScheduler, SimpleScheduler
from astroplan.periodic import EclipsingSystem
from astroplan.constraints import is_event_observable


dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

constraints = [AltitudeConstraint(min=24*u.deg), AtNightConstraint()] #MoonSeparationConstraint(min=30*u.deg)


# def scheduled_table(t,telescope,SS1_night_blocks_exists):
#
#     Path='/Users/elsaducrot/Documents/GitHub/Scheduler_global/Python'
#     tel = telescope
#     try:
#         os.path.exists(os.path.join(Path,tel,'night_blocks_' + tel + '_' +  Time(t).tt.datetime.strftime("%Y-%m-%d ") + '.txt'))
#         print(os.path.join(Path,tel,'night_blocks_' + tel + '_' +  Time(t).tt.datetime.strftime("%Y-%m-%d ") + '.txt'))
#     except NameError:
#         print('no input night_block for this day')
#     if SS1_night_blocks_exists:
#         SS1_night_blocks_old=SS1_night_blocks
#         scheduled_table=SS1_night_blocks_old
#     else:
#         scheduled_table=Table.read(os.path.join(Path,tel,'night_blocks_' + tel + '_' +  Time(t).tt.datetime.strftime("%Y-%m-%d ") + '.txt'), format='ascii')
#     return scheduled_table
#     # print('scheduled_table ---->',scheduled_table)
#
#     # tel='Europa'
#     t_now=time.strftime("%Y-%m-%d 12:00:00")
#     today = input('Today? (y/n)')
#     if today=='y':
#         t_now=time.strftime("%Y-%m-%d 12:00:00")
#     else:
#         t_now = input('Enter date (%Y-%m-%d 12:00:00) : ')
#     #Date different than todays's date
#     # t_now='2018-11-01 12:00:00.000'
#
#     nb_jours=int(input('nb days to prepare:'))
#     t0=Time(t_now)
#     # print(t0)
#
#
# def transit_follow_up(t_now,t,observatory,constraints,targets_transit,name_spc_transit,timings_transit,T0_err_transit,periods_transit,P_err_transit,durations_transit,W_err_transit,filt_spc_transit,tesxpspe_spc_transit):
# ########################### TRANSITS FOLLOW UP ##################################
#     observing_time = Time(t_now)
#     t=Time(t_now)
#     nb_jours=1#input('nb jours for follow up:')
#     observing_time_end =Time(t+dt*nb_jours)
#     # print(observing_time.iso,observing_time_end.iso)
#     constraints = [AltitudeConstraint(min=24*u.deg), AtNightConstraint()] #MoonSeparationConstraint(min=30*u.deg)
#
#     for i,nam in enumerate(name_spc_transit):
#         blocks=[]
#         try:
#             SS1_night_blocks
#         except NameError:
#             SS1_night_blocks_exists = False
#         else:
#             SS1_night_blocks_exists = True
#             SS1_night_blocks_old=SS1_night_blocks
#         epoch = Time(timings_transit[i], format='jd')
#         period = periods_transit[i] * u.day
#         duration = durations_transit[i] * u.day
#         oot_time=1.*u.hour
#         target_transit=EclipsingSystem(primary_eclipse_time=epoch, orbital_period=period, duration=duration,
#                                name=str(nam))
#         print('target_transit',target_transit.next_primary_eclipse_time(observing_time,n_eclipses=1))
#         timing_to_obs_jd=Time(target_transit.next_primary_eclipse_time(observing_time,n_eclipses=1)).jd
#         # print('timing_to_obs_jd',timing_to_obs_jd[0],(timing_to_obs_jd[0]-epoch.jd)/periods_transit[i]*P_err_transit[i])
#
#         target = FixedTarget.from_name(str(nam))
#         n_transits =1#int(((observing_time_end-observing_time)/periods[i]).value)  # This is the roughly number of transits per year
#         try:
#             ing_egr = target_transit.next_primary_ingress_egress_time(observing_time,
#                                                                 n_eclipses=n_transits)
#             # print(ing_egr.iso)
#         except ValueError:
#             print('No transit of ',nam,' on the period chosen')
#
#         observable = is_event_observable(constraints, observatory, target,times_ingress_egress=ing_egr)
#         observable_half=is_event_observable(constraints, observatory, target,times_ingress_egress=ing_egr)
#         if np.any(observable):
#             start_transit=Time(ing_egr[0][0]-T0_err_transit[i]*u.day-W_err_transit[i]*u.day-oot_time-((timing_to_obs_jd[0]-epoch.jd)/periods_transit[i]*P_err_transit[i])*1*u.day)
#             print('start_transit',start_transit.iso)
#             end_transit=Time(ing_egr[0][1]+T0_err_transit[i]*u.day+W_err_transit[i]*u.day+oot_time+((timing_to_obs_jd[0]-epoch.jd)/periods_transit[i]*P_err_transit[i])*1*u.day)
#             print('end_transit',end_transit.iso)
#             dur_obs_transit_target=(end_transit-start_transit).value*1.*u.day
#             constraints_transit_target=[AltitudeConstraint(min=24*u.deg),MoonSeparationConstraint(min=30*u.deg), \
#             TimeConstraint(start_transit,end_transit)]
#             idx_first_target=int(np.where((name_spc_transit==nam))[0])
#             a=ObservingBlock(targets_transit[idx_first_target],dur_obs_transit_target,-1,constraints= constraints_transit_target,configuration={'filt=' + str(filt_spc_transit[idx_first_target]),'texp=' + str(tesxpspe_spc_transit[idx_first_target])})
#             blocks.append(a)
#             transitioner = Transitioner(slew_rate= 11*u.deg/u.second)
#             seq_schedule_SS1=Schedule(t,t+dt)
#             sequen_scheduler_SS1=SPECULOOSScheduler(constraints=constraints_transit_target, observer=observatory,transitioner=transitioner)
#             sequen_scheduler_SS1(blocks,seq_schedule_SS1)
#             print(start_transit.iso,end_transit.iso)
#             if (len(seq_schedule_SS1.to_table()['target'])!=0):
#                 SS1_night_blocks=seq_schedule_SS1.to_table() #Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
#                 return SS1_night_blocks
#         else:
#             print('pas de transit ce jour')
#
# def special_target(t,observatory,input_name,targets,name_spc,filt_spc,tesxpspe_spc):
# ############## Scheduling of the special target ##################################
#     # dt_nautical_civil_evening_2=(observatory.twilight_evening_nautical(t,which='next')-observatory.twilight_evening_civil(t,which='next'))/2
#     # dt_civil_nautical_morning_2=(observatory.twilight_morning_civil(t+dt,which='nearest')-observatory.twilight_morning_nautical(t+dt,which='nearest'))/2
#     night_duration=((Time(observatory.twilight_morning_nautical(t+dt,which='nearest')))-(Time(observatory.twilight_evening_nautical(t,which='next'))))#observatory.twilight_morning_astronomical(time_ranges[month_obs][0]+dt*t,which='next')-observatory.twilight_evening_astronomical(time_ranges[month_obs][0]+dt*t,which='nearest') #night_duration
#     dur_obs_both_target=(night_duration/(2*u.day))*2*u.day
#     # print(dur_obs_both_target)
#     constraints_special_target=[AltitudeConstraint(min=24*u.deg),MoonSeparationConstraint(min=30*u.deg), \
#     TimeConstraint((Time(observatory.twilight_evening_nautical(t,which='next'))), \
#     (Time(observatory.twilight_morning_nautical(t+dt,which='nearest'))))]
#     # print(Time(observatory.twilight_evening_nautical(t,which='next')-dt_nautical_civil_evening_2).iso,Time(observatory.twilight_morning_civil(t+dt,which='nearest') - dt_civil_nautical_morning_2).iso)
#     idx_first_target=int(np.where((name_spc==input_name))[0])
#     # print(idx_first_target)
#     # nb_hours_start_special=nb_hours_observed[idx_first_target]
#     # nb_hours_observed[idx_first_target]+=dur_obs_both_target.value*24
#     # nb_hours[idx_first_target]=[nb_hours_start_special,nb_hours_observed[idx_first_target]]
#     # nb_hours_special_target=nb_hours[idx_first_target]
#     idx_special_target=idx_first_target
#     # print(targets[idx_first_target])
#     blocks=[]
#     a=ObservingBlock(targets[idx_first_target],dur_obs_both_target,-1,constraints= constraints_special_target,configuration={'filt=' + str(filt_spc[idx_first_target]),'texp=' + str(tesxpspe_spc[idx_first_target])})
#     blocks.append(a)
#     transitioner = Transitioner(slew_rate= 11*u.deg/u.second)
#     # print(t.iso,(t+dt).iso)
#     seq_schedule_SS1=Schedule(t,t+dt)
#     sequen_scheduler_SS1=SPECULOOSScheduler(constraints=constraints_special_target, observer=observatory,transitioner=transitioner)
#     sequen_scheduler_SS1(blocks,seq_schedule_SS1)
#     SS1_night_blocks=seq_schedule_SS1.to_table() #Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
#     return SS1_night_blocks
#
# def nb_hours_special(t,observatory,input_name,name_spc,nb_hours_observed,nb_hours):
#     idx_first_target=int(np.where((name_spc==input_name))[0])
#     # dt_nautical_civil_evening_2=(observatory.twilight_evening_nautical(t,which='next')-observatory.twilight_evening_civil(t,which='next'))/2
#     # dt_civil_nautical_morning_2=(observatory.twilight_morning_civil(t+dt,which='nearest')-observatory.twilight_morning_nautical(t+dt,which='nearest'))/2
#     night_duration=((Time(observatory.twilight_morning_nautical(t+dt,which='nearest')))-(Time(observatory.twilight_evening_nautical(t,which='next'))))#observatory.twilight_morning_astronomical(time_ranges[month_obs][0]+dt*t,which='next')-observatory.twilight_evening_astronomical(time_ranges[month_obs][0]+dt*t,which='nearest') #night_duration
#     dur_obs_both_target=(night_duration/(2*u.day))*2*u.day
#     nb_hours_start_special=nb_hours_observed[idx_first_target]
#     nb_hours_observed[idx_first_target]+=dur_obs_both_target.value*24
#     nb_hours[idx_first_target]=[nb_hours_start_special,nb_hours_observed[idx_first_target]]
#     nb_hours_special_target=nb_hours[idx_first_target]
#     return nb_hours_special_target
#
# def monitoring(t,observatory,input_name,name_spc_mon,targets_mon,filt_spc,tesxpspe_spc,time_monitoring,airmass_max,airmax_min):
#     ############################# MONITORING ##############################################
#     # dt_nautical_civil_evening_2=(observatory.twilight_evening_nautical(t,which='next')-observatory.twilight_evening_civil(t,which='next'))/2
#     # dt_civil_nautical_morning_2=(observatory.twilight_morning_civil(t+dt,which='nearest')-observatory.twilight_morning_nautical(t+dt,which='nearest'))/2
#     #Scheduling of the monitiring target
#     idx_first_target=int(np.where((name_spc_mon==input_name))[0])
#     # print('target moni',targets_mon[idx_first_target])
#     dur_mon_target=float(time_monitoring[idx_first_target])*u.minute
#     # print(dur_mon_target)
#     constraints_monitoring_target=[AltitudeConstraint(min=24*u.deg),MoonSeparationConstraint(min=30*u.deg), AirmassConstraint(airmass_max, airmax_min), \
#     TimeConstraint((Time(observatory.twilight_evening_nautical(t,which='next'))), \
#     (Time(observatory.twilight_morning_nautical(t+dt,which='nearest'))))]
#     blocks=[]
#     a=ObservingBlock(targets_mon[idx_first_target],dur_mon_target,-1,constraints= constraints_monitoring_target,configuration={'filt=' + str(filt_spc[idx_first_target]),'texp=' + str(tesxpspe_spc[idx_first_target])})
#     blocks.append(a)
#
#     transitioner = Transitioner(slew_rate= 11*u.deg/u.second)
#     seq_schedule_SS1=Schedule(t,t+dt)
#     sequen_scheduler_SS1=SPECULOOSScheduler(constraints=constraints_monitoring_target, observer=observatory,transitioner=transitioner)
#     sequen_scheduler_SS1(blocks,seq_schedule_SS1)
#     SS1_night_blocks=seq_schedule_SS1.to_table() #Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
#     return SS1_night_blocks
#
# # #######################
# def planification(t,observatory,SS1_night_blocks_exists,scheduled_table,SS1_night_blocks):
#     # dt_nautical_civil_evening_2=(observatory.twilight_evening_nautical(t,which='next')-observatory.twilight_evening_civil(t,which='next'))/2
#     # dt_civil_nautical_morning_2=(observatory.twilight_morning_civil(t+dt,which='nearest')-observatory.twilight_morning_nautical(t+dt,which='nearest'))/2
#
#     name=scheduled_table['target']
#     scheduled_date_start=scheduled_table['start time (UTC)']
#     scheduled_date_end=scheduled_table['end time (UTC)']
#     scheduled_duration=scheduled_table['duration (minutes)']
#     scheduled_ra1=scheduled_table['ra (h)']
#     scheduled_ra2=scheduled_table['ra (m)']
#     scheduled_ra3=scheduled_table['ra (s)']
#     scheduled_dec1=scheduled_table['dec (d)']
#     scheduled_dec2=scheduled_table['dec (m)']
#     scheduled_dec3=scheduled_table['dec (s)']
#     scheduled_configuration=scheduled_table['configuration']
#     pd.DataFrame(scheduled_table.to_pandas())
#
#
#     name_special_targets=SS1_night_blocks['target']
#     Start_special_targets=SS1_night_blocks['start time (UTC)']
#     Finish_special_targets=SS1_night_blocks['end time (UTC)']
#     duration_special_targets=SS1_night_blocks['duration (minutes)']
#     special_targets_ra1=SS1_night_blocks['ra (h)']
#     special_targets_ra2=SS1_night_blocks['ra (m)']
#     special_targets_ra3=SS1_night_blocks['ra (s)']
#     special_targets_dec1=SS1_night_blocks['dec (d)']
#     special_targets_dec2=SS1_night_blocks['dec (m)']
#     special_targets_dec3=SS1_night_blocks['dec (s)']
#     special_targets_configuration=SS1_night_blocks['configuration']
#
#
#     for i,nam in enumerate(name):
#         # print(nam)
#         end_before_cut=scheduled_date_end[i]
#         start_before_cut=scheduled_date_start[i]
#         if SS1_night_blocks_exists:
#             print('Several transit this night')
#             if (np.any(str(SS1_night_blocks['target'][0]).find('Trappist-1'))==0) and (np.any(str(SS1_night_blocks_old['target'][0]).find('Trappist-1'))==0):
#                 print('deux trappist1!')
#                 Start_special_targets[0]=min(SS1_night_blocks['start time (UTC)'][0],SS1_night_blocks_old['start time (UTC)'][0])
#                 Finish_special_targets[0]=max(SS1_night_blocks['end time (UTC)'][0],SS1_night_blocks_old['end time (UTC)'][0])
#         if (Start_special_targets[0]<Time(observatory.twilight_evening_nautical(t,which='next'))):
#             Start_special_targets[0]=Time(observatory.twilight_evening_nautical(t,which='next')).iso
#         if (scheduled_date_start[i] <= Start_special_targets[0]) and (Start_special_targets[0] <= scheduled_date_end[i]):
#             if (Finish_special_targets[0]<= scheduled_date_end[i]): #special target in the middle of a scheduled target block
#                 #first part of scheduled target, change end time and duration
#                 print('case1')
#                 end_before_cut=scheduled_date_end[i]
#                 scheduled_date_end[i]=Start_special_targets[0]
#                 scheduled_duration[i]=(Time(Start_special_targets[0])-Time(scheduled_date_start[i])).value*24*60
#                 #row added for special target
#                 scheduled_table.add_row((name_special_targets[0],Start_special_targets[0],Finish_special_targets[0],duration_special_targets[0],special_targets_ra1[0],special_targets_ra2[0], \
#                 special_targets_ra3[0],special_targets_dec1[0],special_targets_dec2[0],special_targets_dec3[0],special_targets_configuration[0]))
#                 #row added for second part of scheduled target, change start time and duration
#                 scheduled_table.add_row(('_2',Finish_special_targets[0],end_before_cut,(Time(end_before_cut)-Time(Finish_special_targets[0])).value*24*60,scheduled_ra1[i],scheduled_ra2[i],scheduled_ra3[i], \
#                 scheduled_dec1[i],scheduled_dec2[i],scheduled_dec3[i],scheduled_configuration[i]))
#                 # print('ici',scheduled_table)
#             if (Finish_special_targets[0]>= end_before_cut):
#                 print('case2')
#                 #Change end time of scheduled target
#                 scheduled_date_end[i]=Start_special_targets[0]
#                 scheduled_duration[i]=(Time(Start_special_targets[0])-Time(scheduled_date_start[i])).value*24*60
#                 #row added for special target
#                 scheduled_table.add_row((name_special_targets[0],Start_special_targets[0],Finish_special_targets[0],duration_special_targets[0],special_targets_ra1[0],special_targets_ra2[0], \
#                 special_targets_ra3[0],special_targets_dec1[0],special_targets_dec2[0],special_targets_dec3[0],special_targets_configuration[0]))
#         if (start_before_cut >= Start_special_targets[0]):
#             if (Finish_special_targets[0] <= end_before_cut):
#                 print('case3')
#                 if (Finish_special_targets[0] <= start_before_cut):
#                     print('no change')
#                 else:
#                     #Change start time of scheduled target and duration
#                     scheduled_date_start[i]=Finish_special_targets[0]
#                     scheduled_duration[i]=(Time(scheduled_date_end[i])-Time(Finish_special_targets[0])).value*24*60
#                     #delete and replace scheduled target because change of start time and duration
#                     scheduled_table.add_index('target')
#                     index_to_delete=scheduled_table.loc[nam].index
#                     scheduled_table.remove_row(index_to_delete)
#                     scheduled_table.add_row((str(nam),Finish_special_targets[0],end_before_cut,scheduled_duration[i],scheduled_ra1[i],scheduled_ra2[i],scheduled_ra3[i], \
#                     scheduled_dec1[i],scheduled_dec2[i],scheduled_dec3[i],scheduled_configuration[i]))
#                     #row added for special target
#                     scheduled_table.add_row((name_special_targets[0],Start_special_targets[0],Finish_special_targets[0],duration_special_targets[0],special_targets_ra1[0],special_targets_ra2[0], \
#                     special_targets_ra3[0],special_targets_dec1[0],special_targets_dec2[0],special_targets_dec3[0],special_targets_configuration[0]))
#             if (Finish_special_targets[0] >= end_before_cut):
#                 # print('case4')
#                 #delete scheduled target because special target overlap entirely
#                 scheduled_table.add_index('target')
#                 index_to_delete=scheduled_table.loc[nam].index
#                 scheduled_table.remove_row(index_to_delete)
#                 #row added for special target
#                 scheduled_table.add_row((name_special_targets[0],Start_special_targets[0],Finish_special_targets[0],duration_special_targets[0],special_targets_ra1[0],special_targets_ra2[0], \
#                 special_targets_ra3[0],special_targets_dec1[0],special_targets_dec2[0],special_targets_dec3[0],special_targets_configuration[0]))
#
#     # idx_Delete=np.where(scheduled_table['duration (minutes)']<5)[0]
#     # print('IDX DELETE',idx_Delete)
#     # if len(idx_Delete)!=0:
#     #     scheduled_table.remove_row(int(idx_Delete))
#     scheduled_table_unique=unique(scheduled_table, keys='target')
#     index_prio=np.argsort(scheduled_table_unique['start time (UTC)'])
#     scheduled_table_sorted=scheduled_table_unique[index_prio]
#
#     return scheduled_table_sorted
#
# def update_month_plan(t,tel,input_name,name_spc,nb_hours_special_target,nb_hours_threshold,plan_file,scheduled_table_sorted):
#     ##### MODIFY MONTH PLAN #######################################################
#
#     Path='/Users/elsaducrot/Documents/GitHub/Scheduler_global/Python/'
#     idx_special_target=int(np.where((name_spc==input_name))[0])
#     name_final_targets=scheduled_table_sorted['target']
#     Start_final_targets=scheduled_table_sorted['start time (UTC)']
#     Finish_final_targets=scheduled_table_sorted['end time (UTC)']
#     duration_final_targets=scheduled_table_sorted['duration (minutes)']
#     target_table_month_plan = Table.read(plan_file+'.txt', format='ascii')
#     hours_start=target_table_month_plan['nb_hours_start']
#     hours_finish=target_table_month_plan['nb_hours_end']
#     start_times=target_table_month_plan['Start']
#     old_scheduled_table=Table.read(os.path.join(Path,tel,'night_blocks_' + tel + '_' +  Time(t).tt.datetime.strftime("%Y-%m-%d ") + '.txt'), format='ascii')
#     old_scheduled_name=old_scheduled_table['target']
#     old_scheduled_date_start=old_scheduled_table['start time (UTC)']
#     old_scheduled_date_finish=old_scheduled_table['end time (UTC)']
#     try :
#         hours_start_special=str(np.round(nb_hours_special_target[0],3)) + '/' + str(nb_hours_threshold[idx_special_target]) #+ '/' + str(nb_hours_threshold[idx_target[0]]) + ' ' + \
#         hours_finish_special=str(np.round(nb_hours_special_target[1],3)) + '/' + str(nb_hours_threshold[idx_special_target]) + '\n'
#     except NameError:
#         print('not a special target')
#
#     for i,old in enumerate(old_scheduled_date_start):
#         # idx_to_delete=(start_times==old_scheduled_date_start[0])
#         for j in range(0,len(name_final_targets)):
#             if (duration_final_targets[int(j)]>20):
#                 if (name_final_targets[int(j)]==input_name):
#                     target_table_month_plan.add_row([[str(name_final_targets[int(j)])], [str(Start_final_targets[int(j)])], [str(Finish_final_targets[int(j)])], [tel + '_s'],[hours_start_special],[hours_finish_special]])
#                 else:
#                     target_table_month_plan.add_row([[str(name_final_targets[int(j)])], [str(Start_final_targets[int(j)])], [str(Finish_final_targets[int(j)])], [tel + '_s'],[str(hours_start[int(j)])],[str(hours_finish[int(j)])]])
#
#         # target_table_month_plan_unique=unique(target_table_month_plan, keys='Start')
#         index_prio=np.argsort(target_table_month_plan['Start'])
#         target_table_month_plan_sorted=target_table_month_plan[index_prio]
#         target_table_month_plan_sorted.add_index(['Name','Start'])
#         # print(old_scheduled_date_start[int(i)],old_scheduled_name[int(i)])
#         try:
#             ind_delete=(target_table_month_plan_sorted['Start']==old_scheduled_date_start[int(i)]) & (target_table_month_plan_sorted['Name']==old_scheduled_name[int(i)]) & (target_table_month_plan_sorted['Finish']==old_scheduled_date_finish[int(i)])
#             index_to_delete=np.where(ind_delete)[0]
#             target_table_month_plan_sorted.remove_row(int(index_to_delete))
#         except TypeError:
#             print('this block has been modified, no need to delete')
#         # print(ind_delete,target_table_month_plan_sorted[ind_delete])
#         # # print(target_table_month_plan_sorted.loc['2018-11-15 23:45:20.025'])
#         # # index_to_delete=target_table_month_plan_sorted.loc(old_scheduled_date_start[int(i)])#.index
#         # index_to_delete=np.where(ind_delete)[0]
#         # print(np.where(ind_delete)[0])
#         # target_table_month_plan_sorted.remove_row(int(index_to_delete))
#
#     index_prio_tel=np.argsort(target_table_month_plan_sorted['Telescope'])
#     target_table_month_plan_sorted=unique(target_table_month_plan_sorted[index_prio_tel],keys=['Name','Start','Finish','Telescope'])
#     panda_table=target_table_month_plan_sorted.to_pandas()
#     panda_table.to_csv(str(plan_file) + '.txt',sep=' ',index=False)
#     return pd.DataFrame(panda_table)
#
#     ###############################################################################
#

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

    if 'TS_La_Silla' in str(Name):
        location_TSlasilla = EarthLocation.from_geodetic(-70.73000000000002*u.deg, -29.25666666666666*u.deg, 2346.9999999988418*u.m)
        observatories.append(Observer(location=location_TSlasilla, name="TSlasilla", timezone="UTC"))

    if 'TN_Oukaimeden' in str(Name):
        location_TNOuka = EarthLocation.from_geodetic(31.20516*u.deg, -7.862263*u.deg, 2751*u.m)
        observatories.append(Observer(location=location_TNOuka, name="TNOuka", timezone="UTC"))

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

class schedules:
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

    def load_parameters(self, input_file_short_term, nb_observatory):
        with open(input_file_short_term, "r") as f:
            Inputs = yaml.load(f, Loader=yaml.FullLoader)
            self.target_list = Inputs['target_list']
            self.use = Inputs['use']
            df = pd.DataFrame.from_dict(Inputs['observatories'])
            self.observatory = charge_observatories(df[nb_observatory]['name'])[0]
            self.telescopes = df[nb_observatory]['telescopes']
            self.telescope = self.telescopes[0]
            self.day_of_night = Time(Inputs['day_of_night'])  # ,Time(Inputs['date_range'][1])]
            self.start_end_range = Time(Inputs['start_end_range'])
            if self.start_end_range[1] <= self.start_end_range[0]:
                sys.exit('ERROR: end date inferior to start date')
            self.constraints = [AtNightConstraint.twilight_nautical()]
            df = pd.read_csv(self.target_list, delimiter=' ')
            self.target_table_spc = Table.from_pandas(df)
            self.targets = target_list_good_coord_format(self.target_list)

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
        constraints_monitoring_target = [AltitudeConstraint(min=24 * u.deg), MoonSeparationConstraint(min=30 * u.deg),\
                                         AirmassConstraint(max = airmass_max,boolean_constraint=True), \
                                         TimeConstraint((Time(self.observatory.twilight_evening_nautical(self.day_of_night, which='next'))), \
                                                        (Time(self.observatory.twilight_morning_nautical(self.day_of_night+1, which='nearest'))))]
        blocks = []
        a = ObservingBlock(self.targets[idx_first_target], dur_mon_target, -1, constraints=constraints_monitoring_target,\
                           configuration={'filt=' + str(self.target_table_spc['Filter'][idx_first_target]),'texp=' + str(self.target_table_spc['texp_spc'][idx_first_target])})
        blocks.append(a)
        transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
        seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night+1)
        sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_monitoring_target, observer=self.observatory,transitioner=transitioner)
        sequen_scheduler_SS1(blocks, seq_schedule_SS1)
        self.SS1_night_blocks = seq_schedule_SS1.to_table()  # Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
        return self.SS1_night_blocks

    def special_target_with_start_end(self, input_name):
        start = self.start_end_range[0]
        end = self.start_end_range[1]

        if (start >= Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso) or (end <= Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next')).iso):
            sys.exit('WARNING : Start time (or End time) is not on the same day.')

        dur_obs_both_target = (self.night_duration(self.day_of_night) / (2 * u.day)) * 2 * u.day
        constraints_special_target = [AltitudeConstraint(min=24 * u.deg), MoonSeparationConstraint(min=30 * u.deg), \
                                      TimeConstraint(start, end)]
        idx_to_insert_target = int(np.where((self.target_table_spc['Sp_ID'] == input_name))[0])
        blocks = []
        a = ObservingBlock(self.targets[idx_to_insert_target], dur_obs_both_target, -1,constraints=constraints_special_target,
                           configuration={'filt=' + str(self.target_table_spc['Filter'][idx_to_insert_target]),
                                          'texp=' + str(self.target_table_spc['texp_spc'][idx_to_insert_target])})
        blocks.append(a)
        print(a)
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
        blocks=[]
        a = ObservingBlock(self.targets[idx_first_target],dur_obs_both_target,-1,constraints= constraints_special_target,\
                         configuration={'filt=' + str(self.target_table_spc['Filter'][idx_first_target]),'texp=' + str(self.target_table_spc['texp_spc'][idx_first_target])})
        blocks.append(a)
        transitioner = Transitioner(slew_rate= 11*u.deg/u.second)
        seq_schedule_SS1=Schedule(self.day_of_night,self.day_of_night+1)
        sequen_scheduler_SS1=SPECULOOSScheduler(constraints=constraints_special_target, observer = self.observatory,transitioner=transitioner)
        sequen_scheduler_SS1(blocks,seq_schedule_SS1)
        self.SS1_night_blocks=seq_schedule_SS1.to_table()
        if len(self.SS1_night_blocks) == 0:
            print('WARNING = Impossible to schedule target ' + input_name + ' at this time range and/or with those constraints')
        return self.SS1_night_blocks

    def transit_follow_up(self,follow_up_list):
        df = pd.read_csv(follow_up_list,delimiter= ' ')
        constraints = [AltitudeConstraint(min=24 * u.deg), AtNightConstraint()]

        for i in range(len(df['Sp_ID'])):
            blocks = []
            # try:
            #     SS1_night_blocks
            # except NameError:
            #     SS1_night_blocks_exists = False
            # else:
            #     SS1_night_blocks_exists = True
            #     SS1_night_blocks_old = SS1_night_blocks

            epoch = Time(df['T0'][i], format='jd')
            period = df['P'][i] * u.day
            duration = df['W'][i] * u.day
            oot_time = duration.value * 1.5 * u.day
            T0_err_transit = df['T0_err'][i]
            P_err_transit = df['P_err'][i]
            W_err_transit = df['W_err'][i]
            target_transit = EclipsingSystem(primary_eclipse_time = epoch, orbital_period=period, duration=duration, name=df['Sp_ID'][i])
            print('target_transit', target_transit.next_primary_eclipse_time(self.day_of_night, n_eclipses=1))
            timing_to_obs_jd = Time(target_transit.next_primary_eclipse_time(self.day_of_night, n_eclipses=1)).jd

            target = target_list_good_coord_format(follow_up_list)
            n_transits = 1
            try:
                ing_egr = target_transit.next_primary_ingress_egress_time(self.day_of_night,n_eclipses=n_transits)
            except ValueError:
                print('No transit of ', df['Sp_ID'][i], ' on the period chosen')

            observable = is_event_observable(constraints,self.observatory, target, times_ingress_egress=ing_egr)
            if np.any(observable):
                err_T0_neg = timing_to_obs_jd[0] - (np.round((timing_to_obs_jd[0] - epoch.jd) / period.value,1) * (period.value - P_err_transit) + (epoch.jd - T0_err_transit ))
                err_T0_pos = (np.round((timing_to_obs_jd[0] - epoch.jd) / period.value,1) * (period.value + P_err_transit) + (epoch.jd + T0_err_transit )) - timing_to_obs_jd[0]
                start_transit = Time(ing_egr[0][0].value - err_T0_neg - oot_time.value - W_err_transit,format='jd') #- T0_err_transit - W_err_transit  - oot_time.value/24 - (timing_to_obs_jd[0] - epoch.jd) / period.value * P_err_transit,format='jd')
                print('start_transit', start_transit.iso)
                end_transit = Time(ing_egr[0][1].value + err_T0_pos + oot_time.value + W_err_transit,format='jd') #+ T0_err_transit + W_err_transit + oot_time.value/24 + (timing_to_obs_jd[0] - epoch.jd) / period.value * P_err_transit ,format='jd')
                print('end_transit', end_transit.iso)
                dur_obs_transit_target = (end_transit - start_transit).value * 1. * u.day
                constraints_transit_target = [AltitudeConstraint(min=24 * u.deg),
                                              MoonSeparationConstraint(min=30 * u.deg), \
                                              TimeConstraint(start_transit, end_transit)]
                idx_first_target = int(np.where((df['Sp_ID'] == df['Sp_ID'][i]))[0])
                a = ObservingBlock(target[idx_first_target], dur_obs_transit_target, -1,
                                   constraints=constraints_transit_target,
                                   configuration={'filt=' + str(df['Filter'][idx_first_target]),
                                                  'texp=' + str(df['texp_spc'][idx_first_target])})
                blocks.append(a)
                transitioner = Transitioner(slew_rate=11 * u.deg / u.second)
                seq_schedule_SS1 = Schedule(self.day_of_night, self.day_of_night+1)
                sequen_scheduler_SS1 = SPECULOOSScheduler(constraints=constraints_transit_target, observer=self.observatory,
                                                          transitioner=transitioner)
                sequen_scheduler_SS1(blocks, seq_schedule_SS1)
                print(start_transit.iso, end_transit.iso)
                if len(seq_schedule_SS1.to_table()['target']) != 0:
                    self.SS1_night_blocks = seq_schedule_SS1.to_table()  # Table.read(os.path.join(Path,tel,'special_target_test.txt'), format='ascii')#
                    return self.SS1_night_blocks
            else:
                print('pas de transit ce jour')

    def planification(self):
        for i in range(len(self.scheduled_table['target'])):
            end_before_cut = self.scheduled_table['end time (UTC)'][i]
            start_before_cut = self.scheduled_table['start time (UTC)'][i]

            if self.SS1_night_blocks is None :
                sys.exit('WARNING : No block to insert !')

            if not (self.scheduled_table_sorted is None):
                print('Several transit this night')
                if (np.any(str(self.SS1_night_blocks['target'][0]).find('Trappist-1'))==0) and (np.any(str(self.scheduled_table_sorted['target'][0]).find('Trappist-1'))==0):
                    print('deux trappist1!')
                    self.SS1_night_blocks['start time (UTC)'][0] = min(self.SS1_night_blocks['start time (UTC)'][0],self.scheduled_table_sorted['start time (UTC)'][0])
                    self.SS1_night_blocks['end time (UTC)'][0] = max(self.SS1_night_blocks['end time (UTC)'][0],self.scheduled_table_sorted['end time (UTC)'][0])

            if (self.SS1_night_blocks['start time (UTC)'][0] <= Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next')).iso): #case 1
                print('case1')
                self.SS1_night_blocks['start time (UTC)'][0] = Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next')).iso

            if (self.scheduled_table['start time (UTC)'][i] <= self.SS1_night_blocks['start time (UTC)'][0]) and (self.SS1_night_blocks['start time (UTC)'][0] <= self.scheduled_table['end time (UTC)'][i]):

                if (self.SS1_night_blocks['end time (UTC)'][0] <= self.scheduled_table['end time (UTC)'][i]): #case 2
                    print('case2')
                    self.scheduled_table['end time (UTC)'][i] = self.SS1_night_blocks['start time (UTC)'][0]
                    self.scheduled_table['duration (minutes)'][i] = (Time(self.scheduled_table['end time (UTC)'][i])-Time(self.scheduled_table['start time (UTC)'][i])).value*24*60
                    self.scheduled_table.add_row((self.SS1_night_blocks['target'][0],self.SS1_night_blocks['start time (UTC)'][0],self.SS1_night_blocks['end time (UTC)'][0],self.SS1_night_blocks['duration (minutes)'][0],self.SS1_night_blocks['ra (h)'][0],self.SS1_night_blocks['ra (m)'][0], \
                    self.SS1_night_blocks['ra (s)'][0],self.SS1_night_blocks['dec (d)'][0],self.SS1_night_blocks['dec (m)'][0],self.SS1_night_blocks['dec (s)'][0],self.SS1_night_blocks['configuration'][0]))
                    self.scheduled_table.add_row(('_2',self.SS1_night_blocks['end time (UTC)'][0],end_before_cut,(Time(end_before_cut)-Time(self.SS1_night_blocks['end time (UTC)'][0])).value*24*60,self.scheduled_table['ra (h)'][i],self.scheduled_table['ra (m)'][i],self.scheduled_table['ra (s)'][i], \
                    self.scheduled_table['dec (d)'][i],self.scheduled_table['dec (m)'][i],self.scheduled_table['dec (s)'][i],self.scheduled_table['configuration'][i]))

                elif (self.SS1_night_blocks['end time (UTC)'][0] >= end_before_cut) and \
                        (self.SS1_night_blocks['end time (UTC)'][0] <= Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso) : #case 3
                    print('case3')
                    self.scheduled_table['end time (UTC)'][i] = self.SS1_night_blocks['start time (UTC)'][0]
                    self.scheduled_table['duration (minutes)'][i]=(Time(self.scheduled_table['end time (UTC)'][i])-Time(self.scheduled_table['start time (UTC)'][i])).value*24*60
                    self.scheduled_table.add_row((self.SS1_night_blocks['target'][0],self.SS1_night_blocks['start time (UTC)'][0],self.SS1_night_blocks['end time (UTC)'][0],self.SS1_night_blocks['duration (minutes)'][0],self.SS1_night_blocks['ra (h)'],self.SS1_night_blocks['ra (m)'][0], \
                    self.SS1_night_blocks['ra (s)'][0],self.SS1_night_blocks['dec (d)'][0],self.SS1_night_blocks['dec (m)'][0],self.SS1_night_blocks['dec (s)'][0],self.SS1_night_blocks['configuration'][0]))

                elif (self.SS1_night_blocks['end time (UTC)'][0] >= Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso) : #case 4
                    print('case4')
                    self.SS1_night_blocks['end time (UTC)'][0] = Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest'))[0].iso
                    self.scheduled_table['end time (UTC)'][i] = self.SS1_night_blocks['start time (UTC)'][0]
                    self.scheduled_table['duration (minutes)'][i]=(Time(self.scheduled_table['end time (UTC)'][i])-Time(self.scheduled_table['start time (UTC)'][i])).value*24*60
                    self.scheduled_table.add_row((self.SS1_night_blocks['target'][0],self.SS1_night_blocks['start time (UTC)'][0],self.SS1_night_blocks['end time (UTC)'][0],\
                    self.SS1_night_blocks['duration (minutes)'][0],self.SS1_night_blocks['ra (h)'], self.SS1_night_blocks['ra (m)'][0], self.SS1_night_blocks['ra (s)'][0],\
                    self.SS1_night_blocks['dec (d)'][0],self.SS1_night_blocks['dec (m)'][0],self.SS1_night_blocks['dec (s)'][0],self.SS1_night_blocks['configuration'][0]))

            if (start_before_cut >= self.SS1_night_blocks['start time (UTC)'][0]):

                if (self.SS1_night_blocks['end time (UTC)'][0] <= end_before_cut):
                    if (self.SS1_night_blocks['end time (UTC)'][0] <= start_before_cut):
                        print('no change')
                    elif (self.SS1_night_blocks['start time (UTC)'][0] >= Time(self.observatory.twilight_evening_nautical(self.day_of_night,which='next'))[0].iso): #case 5
                        print('case5')
                        self.scheduled_table['start time (UTC)'][i] = self.SS1_night_blocks['end time (UTC)'][0]
                        self.scheduled_table['duration (minutes)'][i] = (Time(self.scheduled_table['end time (UTC)'][i])-Time(self.scheduled_table['start time (UTC)'][i])).value*24*60
                        self.scheduled_table.add_index('target')
                        index_to_delete = self.scheduled_table.loc[self.scheduled_table['target'][i]].index
                        self.scheduled_table.remove_row(index_to_delete)
                        self.scheduled_table.add_row((str(self.scheduled_table['target'][i]),self.SS1_night_blocks['end time (UTC)'][0],end_before_cut,self.scheduled_table['duration (minutes)'][i],self.scheduled_table['ra (h)'][i],self.scheduled_table['ra (m)'][i],self.scheduled_table['ra (s)'][i], \
                        self.scheduled_table['dec (d)'][0],self.self.scheduled_table['dec (m)'][0],self.self.scheduled_table['dec (s)'][0],self.scheduled_table['configuration'][i]))
                        scheduled_table.add_row((self.SS1_night_blocks['target'][0],self.SS1_night_blocks['start time (UTC)'][0],self.SS1_night_blocks['end time (UTC)'][0],self.SS1_night_blocks['duration (minutes)'][0],self.SS1_night_blocks['ra (h)'][0],self.SS1_night_blocks['ra (m)'][0], \
                        self.SS1_night_blocks['ra (s)'][0],self.SS1_night_blocks['dec (d)'][0],self.SS1_night_blocks['dec (m)'][0],self.SS1_night_blocks['dec (s)'][0],self.SS1_night_blocks['configuration'][0]))

                    elif (self.SS1_night_blocks['start time (UTC)'][0] <= Time(self.observatory.twilight_evening_nautical(self.day_of_night, which='next'))[0].iso):
                        self.SS1_night_blocks['start time (UTC)'][0] = Time(self.observatory.twilight_evening_nautical(self.day_of_night, which='next'))[0].iso
                        self.SS1_night_blocks['duration (minutes)'][i] = (Time(self.SS1_night_blocks['end time (UTC)'][0])-Time(self.SS1_night_blocks['start time (UTC)'][0])).value*24*60
                        self.scheduled_table['start time (UTC)'][i] = self.SS1_night_blocks['end time (UTC)'][0]
                        self.scheduled_table['duration (minutes)'][i] = (Time(self.scheduled_table['end time (UTC)'][i])-Time(self.scheduled_table['start time (UTC)'][i])).value*24*60
                        self.scheduled_table.add_index('target')
                        index_to_delete = self.scheduled_table.loc[self.scheduled_table['target'][i]].index
                        self.scheduled_table.remove_row(index_to_delete)
                        self.scheduled_table.add_row((str(self.scheduled_table['target'][i]),self.SS1_night_blocks['end time (UTC)'][0],end_before_cut,self.scheduled_table['duration (minutes)'][i],self.scheduled_table['ra (h)'][i],self.scheduled_table['ra (m)'][i],self.scheduled_table['ra (s)'][i], \
                        self.scheduled_table['dec (d)'][0],self.self.scheduled_table['dec (m)'][0],self.self.scheduled_table['dec (s)'][0],self.scheduled_table['configuration'][i]))
                        scheduled_table.add_row((self.SS1_night_blocks['target'][0],self.SS1_night_blocks['start time (UTC)'][0],self.SS1_night_blocks['end time (UTC)'][0],self.SS1_night_blocks['duration (minutes)'][0],self.SS1_night_blocks['ra (h)'][0],self.SS1_night_blocks['ra (m)'][0], \
                        self.SS1_night_blocks['ra (s)'][0],self.SS1_night_blocks['dec (d)'][0],self.SS1_night_blocks['dec (m)'][0],self.SS1_night_blocks['dec (s)'][0],self.SS1_night_blocks['configuration'][0]))

                if (self.SS1_night_blocks['end time (UTC)'][0] >= end_before_cut) and (self.SS1_night_blocks['end time (UTC)'][0] <= Time(self.observatory.twilight_morning_nautical(self.day_of_night,which='nearest')).iso):
                    self.scheduled_table.add_index('target')
                    index_to_delete=self.scheduled_table.loc[self.scheduled_table['target'][i]].index
                    self.scheduled_table.remove_row(index_to_delete)
                    self.scheduled_table.add_row((self.SS1_night_blocks['target'][0],self.SS1_night_blocks['start time (UTC)'][0],self.SS1_night_blocks['end time (UTC)'][0],self.SS1_night_blocks['duration (minutes)'][0],self.SS1_night_blocks['ra (h)'][0],self.SS1_night_blocks['ra (m)'][0], \
                    self.SS1_night_blocks['ra (s)'][0],self.SS1_night_blocks['dec (d)'][0],self.SS1_night_blocks['dec (m)'][0],self.SS1_night_blocks['dec (s)'][0],self.SS1_night_blocks['configuration'][0]))

                elif  (self.SS1_night_blocks['end time (UTC)'][0] >= Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest')).iso):
                    self.SS1_night_blocks['end time (UTC)'][0] = Time(self.observatory.twilight_morning_nautical(self.day_of_night+1,which='nearest'))[0].iso
                    self.scheduled_table.add_index('target')
                    index_to_delete = self.scheduled_table.loc[self.scheduled_table['target'][i]].index
                    self.scheduled_table.remove_row(index_to_delete)
                    self.scheduled_table.add_row((self.SS1_night_blocks['target'][0],
                                                  self.SS1_night_blocks['start time (UTC)'][0],
                                                  self.SS1_night_blocks['end time (UTC)'][0],
                                                  self.SS1_night_blocks['duration (minutes)'][0],
                                                  self.SS1_night_blocks['ra (h)'][0],
                                                  self.SS1_night_blocks['ra (m)'][0],
                                                  self.SS1_night_blocks['ra (s)'][0],
                                                  self.SS1_night_blocks['dec (d)'][0],
                                                  self.SS1_night_blocks['dec (m)'][0],
                                                  self.SS1_night_blocks['dec (s)'][0],
                                                  self.SS1_night_blocks['configuration'][0]))

        self.scheduled_table_unique=unique(self.scheduled_table, keys='target')
        index_prio=np.argsort(self.scheduled_table_unique['start time (UTC)'])
        self.scheduled_table_sorted = self.scheduled_table_unique[index_prio]

        return self.scheduled_table_sorted

    def make_scheduled_table(self):
        Path='/Users/elsaducrot/Documents/GitHub/Scheduler_global/Python'
        try:
            os.path.exists(os.path.join(Path,self.telescope,'night_blocks_' + self.telescope + '_' +  self.day_of_night.tt.datetime[0].strftime("%Y-%m-%d")+ '.txt'))
            print(os.path.join(Path,self.telescope,'night_blocks_' + self.telescope + '_' +  self.day_of_night.tt.datetime[0].strftime("%Y-%m-%d") + '.txt'))
        except NameError:
            print('no input night_block for this day')
        except FileNotFoundError:
            print('no input night_block for this day')

        if not (self.scheduled_table is None):
            return self.scheduled_table
        else:
            self.scheduled_table = Table.read(os.path.join(Path,self.telescope,'night_blocks_' + self.telescope + '_' +  self.day_of_night.tt.datetime[0].strftime("%Y-%m-%d") + '.txt'), format='ascii')
            return self.scheduled_table

    def make_night_block(self):

        Path='/Users/elsaducrot/code/spock'
        if not (self.scheduled_table_sorted is None):
            try:
                self.scheduled_table_sorted.add_index('target')
                index_to_delete = self.scheduled_table_sorted.loc['TransitionBlock'].index
                self.scheduled_table_sorted.remove_row(index_to_delete)
            except KeyError:
                print('no transition block')
            panda_table = self.scheduled_table_sorted.to_pandas()
            panda_table.to_csv(os.path.join(Path,'night_blocks_' + self.telescope + '_' +  self.day_of_night.tt.datetime[0].strftime("%Y-%m-%d") + '.txt'),sep=' ')

    def exposure_time(self, input_name):
        i = np.where((self.target_table_spc['Sp_ID'] == input_name))[0]
        if int(self.target_table_spc['SpT'][i]) <= 9:
            spt_type = 'M' + str(int(self.target_table_spc['SpT'][i]))
        elif (int(self.target_table_spc['SpT'][i]) == 12) or (int(self.target_table_spc['SpT'][i]) == 15) or (
                int(self.target_table_spc['SpT'][i]) == 18):
            spt_type = 'M' + str(int(self.target_table_spc['SpT'][i]) - 10)
        elif int(self.target_table_spc['SpT'][i]) == 10:
            spt_type = 'M9'
        elif int(self.target_table_spc['SpT'][i]) == 11:
            spt_type = 'L2'
        elif int(self.target_table_spc['SpT'][i]) == 13:
            spt_type = 'L2'
        elif int(self.target_table_spc['SpT'][i]) == 14:
            spt_type = 'L5'
        filt_ = str(self.target_table_spc['Filter'][i])
        if (filt_ == 'z\'') or (filt_ == 'r\'') or (filt_ == 'i\'') or (filt_ == 'g\''):
            filt_ = filt_.replace('\'', '')
        a = (ETC.etc(mag_val=self.target_table_spc['J'][i], mag_band='J', spt=spt_type, filt=filt_, airmass=1.2,
                     moonphase=0.6, irtf=0.8, num_tel=1, seeing=0.95))
        texp = a.exp_time_calculator(ADUpeak=30000)[0]
        return texp

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

def save_schedule(input_file,nb_observatory,save=False,over_write =True):
    with open(input_file, "r") as f:
        Inputs = yaml.load(f, Loader=yaml.FullLoader)
        df = pd.DataFrame.from_dict(Inputs['observatories'])
        observatory = charge_observatories(df[nb_observatory]['name'])[0]
        date_range = Time(Inputs['date_range'])
        date_range_in_days = int((date_range[1]- date_range[0]).value)
        telescope = df[nb_observatory]['telescopes'][0]
    for i in range(0,date_range_in_days):
        day = date_range[0] + i
        if save:
            source = './' + 'night_blocks_' + telescope + '_' +  day.tt.datetime.strftime("%Y-%m-%d") + '.txt'
            destination = '/Users/elsaducrot/Documents/GitHub/Scheduler_global/Python/' + telescope + '/'
            destination_2 = '/Users/elsaducrot/Documents/GitHub/Scheduler_global/Python/' + telescope + '/' + 'Archive_night_blocks/'
            if over_write:
                dest = shutil.copy(source, destination)
                dest2 = shutil.copy(source, destination_2)
                print('INFO : ' + '\"' + source + '\"' + ' has been over-written to ' + '\"' +  destination + '\"' )
                print('INFO : ' + '\"' + source + '\"' + ' has been copied to ' + '\"' + destination2 + '\"')
            if not over_write:
                try:
                    dest = shutil.move(source, destination)
                    dest2 = shutil.move(source, destination_2)
                    print('INFO : ' + '\"' +  source + '\"' +  ' has been copied to ' + '\"' + destination + '\"' )
                    print('INFO : ' + '\"' + source + '\"' + ' has been copied to ' + '\"' + destination2 + '\"')
                except shutil.Error:
                    print('INFO : ' + '\"' + destination + 'night_blocks_' + telescope + '_' +  day.tt.datetime.strftime("%Y-%m-%d") + '.txt' + '\"' +  ' already exists')
        if not save:
            print('INFO : Those plans have not been saved')


def make_plans(day, nb_days, telescope):
    make_np(day, nb_days, telescope)

def upload_plans(day, nb_days, telescope):
    if telescope.find('Callisto') is not -1:
        upload_np_calli(day, nb_days)
    if telescope.find('Ganymede') is not -1:
        upload_np_gany(day, nb_days)
    if telescope.find('Io') is not -1:
        upload_np_io(day, nb_days)
    if telescope.find('Europa') is not -1:
        upload_np_euro(day, nb_days)
    if telescope.find('Artemis') is not -1:
        upload_np_artemis(day, nb_days)

    # ------------------- update archive date by date plans folder  ------------------

    path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', telescope,'schedule')
    print(path_database)
    path_plans = os.path.join('/Users/elsaducrot/Documents/GitHub/Scheduler_global/Python/', telescope,'Plans_by_date/')
    print(path_plans)
    subprocess.Popen(["sshpass", "-p", 'eij7iaXi', "scp", "-r", path_plans, path_database])

    # ------------------- update archive niht blocks ------------------

    path_night_blocks = os.path.join('/Users/elsaducrot/Documents/GitHub/Scheduler_global/Python/', telescope,'Archive_night_blocks/')
    print(path_night_blocks)
    subprocess.Popen(["sshpass", "-p", 'eij7iaXi', "scp", "-r", path_night_blocks, path_database])
