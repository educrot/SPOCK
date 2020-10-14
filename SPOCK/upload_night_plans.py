#!/usr/bin/python

import subprocess
import os
import glob
from astropy.time import Time

pwd_appcs  =  'eij7iaXi'
pwd_HUB = 'Felka_5'
pwd_SNO_Reduc1 = 'SNO_Reduc_1'

def upload_np_euro(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cambridge server
        ## Plans by date
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Europa',
                                     'schedule')
        path_plans = os.path.join('./DATABASE/', 'Europa',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive_night_blocks
        path_night_blocks = os.path.join('./DATABASE/', 'Europa',
                                         'Archive_night_blocks/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Europa',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join('./DATABASE/', 'Europa',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        #upload on data reduction computer
        ## cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Europa',
                                     'schedule','Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join('./DATABASE/', 'Europa',
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
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
                                     'schedule')
        path_plans = os.path.join('./DATABASE/', 'Callisto',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive_night_blocks
        path_night_blocks = os.path.join('./DATABASE/', 'Callisto',
                                         'Archive_night_blocks/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join('./DATABASE/', 'Callisto',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        #upload on data reduction computer
        ## cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
                                     'schedule','Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join('./DATABASE/', 'Callisto',
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
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
                                     'schedule')
        path_plans = os.path.join('./DATABASE/', 'Io',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive_night_blocks
        path_night_blocks = os.path.join('./DATABASE/', 'Io',
                                         'Archive_night_blocks/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join('./DATABASE/', 'Io',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        #upload on data reduction computer
        ## cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
                                     'schedule','Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join('./DATABASE/', 'Io',
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
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
                                     'schedule')
        path_plans = os.path.join('./DATABASE/', 'Ganymede',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive_night_blocks
        path_night_blocks = os.path.join('./DATABASE/', 'Ganymede',
                                         'Archive_night_blocks/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join('./DATABASE/', 'Ganymede',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        #upload on data reduction computer
        ## cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
                                     'schedule','Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join('./DATABASE/', 'Ganymede',
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
        path_plans = os.path.join('./DATABASE/', 'Artemis',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive night blocks
        path_night_blocks = os.path.join('./DATABASE/', 'Artemis',
                                         'Archive_night_blocks/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Artemis',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join('./DATABASE/', 'Artemis',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        #upload on data reduction computer
        ## cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Artemis',
                                     'schedule','Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join('./DATABASE/', 'Artemis',
                                  'Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_SNO_Reduc1,"scp", path_database_zip_file,path_local_zip_folder])
        ## Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_SNO_Reduc1,"scp", path_local_zip_folder, 'speculoos@172.16.3.11:/home/speculoos/Desktop/Plans/'])
        print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Artemis')

        # #upload on reduction computer
        # Path= './DATABASE/' + 'Artemis'
        # p = os.path.join(Path,str(t_now)+'.zip')
        # p = subprocess.Popen(["sshpass", "-p", pwd_SNO_Reduc1,"scp", p, 'speculoos@172.16.3.11:/home/speculoos/Desktop/Plans/'])
        # print('----->',t_now,'uploaded on the HUB for Artemis')


def upload_np_ts(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        #upload on Cam server
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'TS_La_Silla',
                                     'schedule')
        ## Plans
        path_plans = os.path.join('./DATABASE/', 'TS_La_Silla',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ##Archive night blocks
        path_night_blocks = os.path.join('./DATABASE/', 'TS_La_Silla',
                                         'Archive_night_blocks/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        # ## zip_files
        # path_local_zip_file = os.path.join('./DATABASE/', 'TS_La_Silla',
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
        path_plans = os.path.join('./DATABASE/', 'TN_Oukaimeden',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ##Archive night blocks
        path_night_blocks = os.path.join('./DATABASE/', 'TN_Oukaimeden',
                                         'Archive_night_blocks/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        # ## zip_files
        # path_local_zip_file = os.path.join('./DATABASE/', 'TN_Oukaimeden',
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
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Saint-Ex',
                                     'schedule')
        path_plans = os.path.join('./DATABASE/', 'Saint-Ex',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        ## Archive_night_blocks
        path_night_blocks = os.path.join('./DATABASE/', 'Saint-Ex',
                                         'Archive_night_blocks/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        ## zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Saint-Ex',
                                     'schedule','Zip_files')
        path_local_zip_file = os.path.join('./DATABASE/', 'Saint-Ex',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        # #upload on data reduction computer
        # ## cam server to local
        # path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Saint-Ex',
        #                              'schedule','Zip_files',str(t_now)+'.zip')
        # path_local_zip_folder = os.path.join('./DATABASE/', 'Saint-Ex',
        #                           'Zip_files/')
        # p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_database_zip_file,path_local_zip_folder])
        # ## Local to reduction computer
        # p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_local_zip_file, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Saint-Ex/Plans/'])
        # print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Saint-Ex')
