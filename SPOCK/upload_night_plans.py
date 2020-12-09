#!/usr/bin/python

import subprocess
import os
from astropy.time import Time
import requests
import yaml

# ************************ Read passwords ************************

with open('passwords.csv', "r") as f:
    Inputs = yaml.load(f, Loader=yaml.FullLoader)
    pwd_appcs = Inputs['pwd_appcs'][0]
    pwd_HUB = Inputs['pwd_HUB'][0]
    user_portal = Inputs['user_portal'][0]
    pwd_portal = Inputs['pwd_portal'][0]
    pwd_appcs = Inputs['pwd_appcs'][0]
    pwd_appcs = Inputs['pwd_appcs'][0]
    pwd_SNO_Reduc1 = Inputs['pwd_SNO_Reduc1'][0]
    user_chart_studio = Inputs['user_chart_studio'][0]
    pwd_chart_studio = Inputs['pwd_chart_studio'][0]
    path_spock = Inputs['path_spock'][0]

# ************************ Create database ************************

telescopes_names = ['Io', 'Europa', 'Ganymede', 'Callisto', 'Artemis', 'Saint-Ex', 'TS_La_Silla', 'TN_Oukaimeden']
if not os.path.exists(path_spock + '/target_lists'):
    os.makedirs(path_spock + '/target_lists')
if not os.path.exists(path_spock + '/survey_hours'):
    os.makedirs(path_spock + '/survey_hours')
if not os.path.exists(path_spock + '/DATABASE'):
    os.makedirs(path_spock + '/DATABASE')
if not os.path.exists(path_spock + '/night_blocks_propositions'):
    os.makedirs(path_spock + '/night_blocks_propositions')
if not os.path.exists(path_spock + '/SPOCK_files'):
    os.makedirs(path_spock + '/SPOCK_files')
for tel in telescopes_names:
    if not os.path.exists(path_spock + '/DATABASE/' + tel):
        os.makedirs(path_spock + '/DATABASE/' + tel)
for tel in telescopes_names:
    if not os.path.exists(path_spock + '/DATABASE/' + tel + '/Archive_night_blocks'):
        os.makedirs(path_spock + '/DATABASE/' + tel + '/Archive_night_blocks')
for tel in telescopes_names:
    if not os.path.exists(path_spock + '/DATABASE/' + tel + '/Plans_by_date'):
        os.makedirs(path_spock + '/DATABASE/' + tel + '/Plans_by_date')
for tel in telescopes_names:
    if not os.path.exists(path_spock + '/DATABASE/' + tel + '/Zip_files'):
        os.makedirs(path_spock + '/DATABASE/' + tel + '/Zip_files')

# ************************ Read target lists from server ************************
target_lists = ['speculoos_target_list_v6.txt', 'target_list_special.txt', 'target_transit_follow_up.txt']
for t_list in target_lists:
    target_list_url = "http://www.mrao.cam.ac.uk/SPECULOOS/spock_files/target_lists/" + t_list
    resp = requests.get(target_list_url, auth=(user_portal, pwd_portal))
    content = resp.text.replace("\n", "")
    open(path_spock + '/target_lists/' + t_list, 'wb').write(resp.content)

survey_hours = ['ObservationHours_Saint-Ex.txt', 'ObservationHours_TRAPPIST.txt', 'ObservationHours.txt']
for file in survey_hours:
    target_list_url = "http://www.mrao.cam.ac.uk/SPECULOOS/spock_files/survey_hours/" + file
    resp = requests.get(target_list_url, auth=(user_portal, pwd_portal))
    content = resp.text.replace("\n", "")
    open(path_spock + '/survey_hours/' + file, 'wb').write(resp.content)

# **********************************************************************************************************


def upload_np_euro(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cambridge server
        ## Plans by date
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Europa',
                                     'schedule','Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Europa',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive_night_blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'Europa',
                                           'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Europa',
                                         'Archive_night_blocks/','night_blocks_Europa_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Europa',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Europa',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        #upload on data reduction computer
        ## cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Europa',
                                     'schedule','Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Europa',
                                  'Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_database_zip_file,path_local_zip_folder])
        ## Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_local_zip_file, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Europa/Plans/'])
        print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Europa')

        # Empty local folder if needed
        # files = glob.glob(path_local_zip_folder)
        # for f in files:
        #     os.remove(f)

def upload_np_calli(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cambridge server
        ## Plans by date
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
                                     'schedule','Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Callisto',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive_night_blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'Callisto',
                                           'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Callisto',
                                         'Archive_night_blocks/','night_blocks_Callisto_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Callisto',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        #upload on data reduction computer
        ## cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
                                     'schedule','Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Callisto',
                                  'Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_database_zip_file,path_local_zip_folder])
        ## Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_local_zip_file, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Callisto/Plans/'])
        print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Callisto')


def upload_np_io(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cambridge server
        ## Plans by date
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
                                     'schedule','Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Io',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive_night_blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'Io',
                                           'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Io',
                                         'Archive_night_blocks/','night_blocks_Io_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Io',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        #upload on data reduction computer
        ## cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
                                     'schedule','Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Io',
                                  'Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_database_zip_file,path_local_zip_folder])
        ## Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_local_zip_file, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Io/Plans/'])
        print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Io')

def upload_np_gany(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cambridge server
        ## Plans by date
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
                                     'schedule','Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive_night_blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'Ganymede',
                                           'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
                                         'Archive_night_blocks/','night_blocks_Ganymede_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        #upload on data reduction computer
        ## cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
                                     'schedule','Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
                                  'Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_database_zip_file,path_local_zip_folder])
        ## Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_local_zip_file, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Ganymede/Plans/'])
        print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Ganymede')

def upload_np_artemis(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cam server
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Artemis',
                                     'schedule')
        ##Plans
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Artemis',
                                     'schedule','Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Artemis',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive night blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'Artemis',
                                           'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Artemis',
                                         'Archive_night_blocks/','night_blocks_Artemis_'+str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Artemis',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Artemis',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        #upload on data reduction computer
        ## cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Artemis',
                                     'schedule','Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Artemis',
                                  'Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_SNO_Reduc1,"scp", path_database_zip_file,path_local_zip_folder])
        ## Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_SNO_Reduc1,"scp", path_local_zip_file, 'speculoos@172.16.3.11:/home/speculoos/Desktop/Plans/'])
        print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Artemis')


def upload_np_ts(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        #upload on Cam server
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'TS_La_Silla',
                                     'schedule')
        ## Plans
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'TS_La_Silla',
                                     'schedule','Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'TS_La_Silla',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ##Archive night blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'TS_La_Silla',
                                           'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'TS_La_Silla',
                                         'Archive_night_blocks/','night_blocks_TS_La_Silla_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        # ## zip_files
        # path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'TS_La_Silla',
        #                                    'Zip_files/', str(t_now) + '.zip')
        # subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database])
        # print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')


def upload_np_tn(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        #upload on Cam server
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'TN_Oukaimeden',
                                     'schedule')
        ## Plans
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'TN_Oukaimeden',
                                     'schedule','Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'TN_Oukaimeden',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ##Archive night blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'TN_Oukaimeden',
                                           'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'TN_Oukaimeden',
                                         'Archive_night_blocks/','night_blocks_TN_Oukaimeden_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        # ## zip_files
        # path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'TN_Oukaimeden',
        #                                    'Zip_files/', str(t_now) + '.zip')
        # subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database])
        # print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

def upload_np_saint_ex(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cambridge server
        ## Plans by date
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Saint-Ex',
                                     'schedule','Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Saint-Ex',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive_night_blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'Saint-Ex',
                                           'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Saint-Ex',
                                         'Archive_night_blocks/','night_blocks_Saint-Ex_'+str(t_now) +'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Saint-Ex',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Saint-Ex',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        # #upload on data reduction computer
        # ## cam server to local
        # path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Saint-Ex',
        #                              'schedule','Zip_files',str(t_now)+'.zip')
        # path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Saint-Ex',
        #                           'Zip_files/')
        # p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_database_zip_file,path_local_zip_folder])
        # ## Local to reduction computer
        # p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_local_zip_file, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Saint-Ex/Plans/'])
        # print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Saint-Ex')
