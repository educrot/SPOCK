#!/usr/bin/python

import subprocess
import os
import time
import shutil
from astropy.time import Time

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
