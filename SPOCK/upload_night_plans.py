#!/usr/bin/python

import subprocess
import os
from astropy.time import Time

pwd_appcs  =  'eij7iaXi'

def upload_np_euro(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso
        Path= './DATABASE' + '/Europa'
        p=os.path.join(Path,str(t_now)+'.zip')

        p = subprocess.Popen(["sshpass", "-p", 'Felka_5',"scp", p, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Europa/Plans/'])

        sts = os.waitpid(p.pid, 0)

        #update log file
        # newtext='Zip files transfered to MacEuropa'
        #
        # with open('log_update_np.txt','a') as out:
        #     out.write(newtext)
        print('----->',t_now,'uploaded on the HUB for Europa')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Europa',
                                     'schedule')
        path_plans = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'Europa',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])

        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/','Europa',
                                     'schedule')
        path_night_blocks = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'Europa',
                                         'Archive_night_blocks/')

        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])

        print('----->', t_now, 'Night plans uploaded on the Cambridge server')


def upload_np_calli(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso
        Path= './DATABASE' + '/Callisto'
        p=os.path.join(Path,str(t_now)+'.zip')

        p = subprocess.Popen(["sshpass", "-p", 'Felka_5',"scp", p, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Callisto/Plans/'])

        sts = os.waitpid(p.pid, 0)
        print('----->',t_now,'uploaded on the HUB for Callisto')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
                                     'schedule')
        path_plans = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'Callisto',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])

        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/','Callisto',
                                     'schedule')
        path_night_blocks = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'Callisto',
                                         'Archive_night_blocks/')

        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])

        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

def upload_np_io(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso
        Path= './DATABASE' + '/Io'
        p=os.path.join(Path,str(t_now)+'.zip')

        p = subprocess.Popen(["sshpass", "-p", 'Felka_5',"scp", p, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Io/Plans/'])

        sts = os.waitpid(p.pid, 0)

        #update log file
        # newtext='Zip files transfered to Macio'
        #
        # with open('log_update_np.txt','a') as out:
        #     out.write(newtext)
        print('----->',t_now,'uploaded on the HUB for Io')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
                                     'schedule')
        path_plans = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'Io',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])

        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/','Io',
                                     'schedule')
        path_night_blocks = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'Io',
                                         'Archive_night_blocks/')

        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])

        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

def upload_np_gany(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso
        Path= './DATABASE' + '/Ganymede'
        p=os.path.join(Path,str(t_now)+'.zip')

        p = subprocess.Popen(["sshpass", "-p", 'Felka_5',"scp", p, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Ganymede/Plans/'])

        sts = os.waitpid(p.pid, 0)

        #update log file
        # newtext='Zip files transfered to MacGanymede'
        #
        # with open('log_update_np.txt','a') as out:
        #     out.write(newtext)
        print('----->',t_now,'uploaded on the HUB for Ganymede')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
                                     'schedule')
        path_plans = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'Ganymede',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])

        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/','Ganymede',
                                     'schedule')
        path_night_blocks = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'Ganymede',
                                         'Archive_night_blocks/')

        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])

        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

def upload_np_artemis(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso
        Path= './DATABASE/' + 'Artemis'
        p=os.path.join(Path,str(t_now)+'.zip')

        p = subprocess.Popen(["sshpass", "-p", 'SNO_Reduc_1',"scp", p, 'speculoos@172.16.3.11:/home/speculoos/Desktop/Plans/'])

        sts = os.waitpid(p.pid, 0)

        print('----->',t_now,'uploaded on the HUB for Artemis')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Artemis',
                                     'schedule')
        path_plans = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'Artemis',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])

        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/','Artemis',
                                     'schedule')
        path_night_blocks = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'Artemis',
                                         'Archive_night_blocks/')

        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])

        print('----->', t_now, 'Night plans uploaded on the Cambridge server')


def upload_np_ts(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'TS_La_Silla',
                                     'schedule')
        path_plans = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'TS_La_Silla',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])

        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/','TS_La_Silla',
                                     'schedule')
        path_night_blocks = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'TS_La_Silla',
                                         'Archive_night_blocks/')

        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])

        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

def upload_np_tn(t_now,nb_jours):
    t0=Time(t_now)
    dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day

    for nb_day in range(0,nb_jours):
        t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'TN_Oukaimeden',
                                     'schedule')
        path_plans = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'TN_Oukaimeden',
                                  'Plans_by_date/')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database])

        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/','TN_Oukaimeden',
                                     'schedule')
        path_night_blocks = os.path.join('/Users/elsaducrot/spock_2/DATABASE/', 'TN_Oukaimeden',
                                         'Archive_night_blocks/')

        subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", "-r", path_night_blocks, path_database])

        print('----->', t_now, 'Night plans uploaded on the Cambridge server')
