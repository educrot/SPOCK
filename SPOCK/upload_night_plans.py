#!/usr/bin/python
import subprocess
import os
from astropy.time import Time
from SPOCK import pwd_HUB, pwd_appcs, pwd_SNO_Reduc1, path_spock


def upload_np_euro(t_now, nb_days):
    t0 = Time(t_now)
    dt = Time('2018-01-02 00:00:00', scale='tcg')-Time('2018-01-01 00:00:00', scale='tcg')  # 1 day

    for nb_day in range(0, nb_days):
        t_now = Time(t0+nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cambridge server
        # Plans by date
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Europa',
                                           'schedule', 'Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Europa',
                                  'Plans_by_date/', str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r", path_plans, path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        # Archive_night_blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                            'Europa',
                                            'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Europa',
                                         'Archive_night_blocks/', 'night_blocks_Europa_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        # zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 
                                               'Europa', 'schedule', 'Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Europa',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r",
                          path_local_zip_file, path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        # upload on HUB
        # cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 
                                              'Europa', 'schedule', 'Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Europa','Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB, "scp", path_database_zip_file,
                              path_local_zip_folder])
        # Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB, "scp", path_local_zip_file, 
                              'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Europa/Plans/'])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the HUB for Europa')


def upload_np_calli(t_now, nb_days):
    t0 = Time(t_now)
    dt = Time('2018-01-02 00:00:00', scale='tcg')-Time('2018-01-01 00:00:00', scale='tcg')  # 1 day

    for nb_day in range(0, nb_days):
        t_now = Time(t0+nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cambridge server
        # Plans by date
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
                                           'schedule', 'Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Callisto',
                                  'Plans_by_date/', str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r", path_plans, path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        # Archive_night_blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                            'Callisto',
                                            'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Callisto',
                                         'Archive_night_blocks/', 'night_blocks_Callisto_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        # zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                               'Callisto', 'schedule', 'Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Callisto',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r",
                          path_local_zip_file, path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        # upload on HUB
        # cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                              'Callisto', 'schedule', 'Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Callisto','Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB, "scp", path_database_zip_file,
                              path_local_zip_folder])
        # Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB, "scp", path_local_zip_file,
                              'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Callisto/Plans/'])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the HUB for Callisto')


def upload_np_io(t_now, nb_days):
    t0 = Time(t_now)
    dt = Time('2018-01-02 00:00:00', scale='tcg')-Time('2018-01-01 00:00:00', scale='tcg')  # 1 day

    for nb_day in range(0, nb_days):
        t_now = Time(t0+nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cambridge server
        # Plans by date
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
                                           'schedule', 'Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Io',
                                  'Plans_by_date/', str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r", path_plans, path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        # Archive_night_blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                            'Io',
                                            'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Io',
                                         'Archive_night_blocks/', 'night_blocks_Io_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        # zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                               'Io', 'schedule', 'Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Io',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r",
                          path_local_zip_file, path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        # upload on HUB
        # cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                              'Io', 'schedule', 'Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Io','Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB, "scp", path_database_zip_file,
                              path_local_zip_folder])
        # Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB, "scp", path_local_zip_file,
                              'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Io/Plans/'])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the HUB for Io')


def upload_np_gany(t_now, nb_days):
    t0 = Time(t_now)
    dt = Time('2018-01-02 00:00:00', scale='tcg')-Time('2018-01-01 00:00:00', scale='tcg')  # 1 day

    for nb_day in range(0, nb_days):
        t_now = Time(t0+nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cambridge server
        # Plans by date
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
                                           'schedule', 'Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
                                  'Plans_by_date/', str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r", path_plans, path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        # Archive_night_blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                            'Ganymede',
                                            'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
                                         'Archive_night_blocks/', 'night_blocks_Ganymede_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        # zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                               'Ganymede', 'schedule', 'Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r",
                          path_local_zip_file, path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        # upload on HUB
        # cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                              'Ganymede', 'schedule', 'Zip_files',str(t_now)+'.zip')
        path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Ganymede','Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB, "scp", path_database_zip_file,
                              path_local_zip_folder])
        # Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_HUB, "scp", path_local_zip_file,
                              'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Ganymede/Plans/'])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the HUB for Ganymede')


def upload_np_artemis(t_now, nb_days):
    t0 = Time(t_now)
    dt = Time('2018-01-02 00:00:00', scale='tcg')-Time('2018-01-01 00:00:00', scale='tcg')  # 1 day

    for nb_day in range(0, nb_days):
        t_now = Time(t0 + nb_day*dt, scale='utc', out_subfmt='date').iso

        # Plans
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'Artemis', 'schedule', 'Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Artemis',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r", path_plans, path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        # Archive night blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                            'Artemis', 'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Artemis',
                                         'Archive_night_blocks/', 'night_blocks_Artemis_'+str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        # zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                               'Artemis', 'schedule', 'Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Artemis',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r", path_local_zip_file,
                          path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')

        # Upload on data reduction computer
        # cam server to local
        path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                              'Artemis', 'schedule', 'Zip_files', str(t_now)+'.zip')
        path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Artemis', 'Zip_files/')
        p = subprocess.Popen(["sshpass", "-p", pwd_SNO_Reduc1, "scp", path_database_zip_file,
                              path_local_zip_folder])
        # Local to reduction computer
        p = subprocess.Popen(["sshpass", "-p", pwd_SNO_Reduc1, "scp", path_local_zip_file,
                              'speculoos@172.16.3.11:/home/speculoos/Desktop/Plans/'])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the HUB for Artemis')


def upload_np_ts(t_now, nb_days):
    t0 = Time(t_now)
    dt = Time('2018-01-02 00:00:00', scale='tcg')-Time('2018-01-01 00:00:00', scale='tcg')  # 1 day

    for nb_day in range(0, nb_days):
        t_now = Time(t0+nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cam server
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                     'TS_La_Silla', 'schedule')
        # Plans
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'TS_La_Silla', 'schedule', 'Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'TS_La_Silla',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        # Archive night blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                            'TS_La_Silla', 'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'TS_La_Silla',
                                         'Archive_night_blocks/','night_blocks_TS_La_Silla_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')


def upload_np_tn(t_now, nb_days):
    t0 = Time(t_now)
    dt = Time('2018-01-02 00:00:00', scale='tcg')-Time('2018-01-01 00:00:00', scale='tcg')  # 1 day

    for nb_day in range(0, nb_days):
        t_now = Time(t0+nb_day*dt, scale='utc', out_subfmt='date').iso

        # upload on Cam server
        path_database = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                     'TN_Oukaimeden', 'schedule')
        # Plans
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'TN_Oukaimeden', 'schedule', 'Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'TN_Oukaimeden',
                                  'Plans_by_date/',str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        # Archive night blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                            'TN_Oukaimeden', 'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'TN_Oukaimeden',
                                         'Archive_night_blocks/','night_blocks_TN_Oukaimeden_'+str(t_now)+'.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')


def upload_np_saint_ex(t_now, nb_days):
    """

    Parameters
    ----------
    t_now
    nb_days

    Returns
    -------

    """
    t0 = Time(t_now)
    dt = Time('2018-01-02 00:00:00', scale='tcg')-Time('2018-01-01 00:00:00', scale='tcg')  # 1 day

    for nb_day in range(0,nb_days):
        t_now = Time(t0+nb_day*dt, scale='utc', out_subfmt='date').iso
        # upload on Cambridge server
        # Plans by date
        path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                           'Saint-Ex', 'schedule', 'Plans_by_date')
        #path_database_plans_SAINT-EX_prince = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
        #                                   'Saint-Ex', 'schedule', 'Plans_by_date')
        path_plans = os.path.join(path_spock + '/DATABASE/', 'Saint-Ex','Plans_by_date/', str(t_now))
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r", path_plans, path_database_plans])
        print('----->', t_now, 'Plans uploaded on the Cambridge server')

        # Archive_night_blocks
        path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                            'Saint-Ex', 'schedule', 'Archive_night_blocks')
        path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Saint-Ex',
                                         'Archive_night_blocks/', 'night_blocks_Saint-Ex_'+str(t_now) + '.txt')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", path_night_blocks, path_database_nightb])
        print('----->', t_now, 'Night plans uploaded on the Cambridge server')

        # zip_files
        path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
                                               'Saint-Ex', 'schedule', 'Zip_files')
        path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Saint-Ex',
                                           'Zip_files/', str(t_now) + '.zip')
        subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp", "-r", path_local_zip_file,
                          path_database_zip_files])
        print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')


# def upload_np_calli_old(t_now,nb_days):
#     t0=Time(t_now)
#     dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day
#
#     for nb_day in range(0,nb_days):
#         t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso
#
#         # upload on Cambridge server
#         ## Plans by date
#         path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
#                                      'schedule','Plans_by_date')
#         path_plans = os.path.join(path_spock + '/DATABASE/', 'Callisto',
#                                   'Plans_by_date/',str(t_now))
#         subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
#         print('----->', t_now, 'Plans uploaded on the Cambridge server')
#
#         ## Archive_night_blocks
#         path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
#                                            'Callisto',
#                                            'schedule', 'Archive_night_blocks')
#         path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Callisto',
#                                          'Archive_night_blocks/','night_blocks_Callisto_'+str(t_now)+'.txt')
#         subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
#         print('----->', t_now, 'Night plans uploaded on the Cambridge server')
#
#         ## zip_files
#         path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
#                                      'schedule','Zip_files')
#         path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Callisto',
#                                            'Zip_files/', str(t_now) + '.zip')
#         subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
#         print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')
#
#         #upload on data reduction computer
#         ## cam server to local
#         path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Callisto',
#                                      'schedule','Zip_files',str(t_now)+'.zip')
#         path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Callisto',
#                                   'Zip_files/')
#         p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_database_zip_file,path_local_zip_folder])
#         ## Local to reduction computer
#         p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_local_zip_file, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Callisto/Plans/'])
#         print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Callisto')

# def upload_np_io_old(t_now,nb_days):
#     t0=Time(t_now)
#     dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day
#
#     for nb_day in range(0,nb_days):
#         t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso
#
#         # upload on Cambridge server
#         ## Plans by date
#         path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
#                                      'schedule','Plans_by_date')
#         path_plans = os.path.join(path_spock + '/DATABASE/', 'Io',
#                                   'Plans_by_date/',str(t_now))
#         subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
#         print('----->', t_now, 'Plans uploaded on the Cambridge server')
#
#         ## Archive_night_blocks
#         path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
#                                            'Io',
#                                            'schedule', 'Archive_night_blocks')
#         path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Io',
#                                          'Archive_night_blocks/','night_blocks_Io_'+str(t_now)+'.txt')
#         subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
#         print('----->', t_now, 'Night plans uploaded on the Cambridge server')
#
#         ## zip_files
#         path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
#                                      'schedule','Zip_files')
#         path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Io',
#                                            'Zip_files/', str(t_now) + '.zip')
#         subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
#         print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')
#
#         #upload on data reduction computer
#         ## cam server to local
#         path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Io',
#                                      'schedule','Zip_files',str(t_now)+'.zip')
#         path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Io',
#                                   'Zip_files/')
#         p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_database_zip_file,path_local_zip_folder])
#         ## Local to reduction computer
#         p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_local_zip_file, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Io/Plans/'])
#         print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Io')
# def upload_np_gany(t_now,nb_days):
#     t0=Time(t_now)
#     dt=Time('2018-01-02 00:00:00',scale='tcg')-Time('2018-01-01 00:00:00',scale='tcg') #1 day
#
#     for nb_day in range(0,nb_days):
#         t_now=Time(t0+ nb_day*dt, scale='utc', out_subfmt='date').iso
#
#         # upload on Cambridge server
#         ## Plans by date
#         path_database_plans = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
#                                      'schedule','Plans_by_date')
#         path_plans = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
#                                   'Plans_by_date/',str(t_now))
#         subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_plans,path_database_plans])
#         print('----->', t_now, 'Plans uploaded on the Cambridge server')
#
#         ## Archive_night_blocks
#         path_database_nightb = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/',
#                                            'Ganymede',
#                                            'schedule', 'Archive_night_blocks')
#         path_night_blocks = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
#                                          'Archive_night_blocks/','night_blocks_Ganymede_'+str(t_now)+'.txt')
#         subprocess.Popen(["sshpass", "-p", pwd_appcs , "scp", path_night_blocks, path_database_nightb])
#         print('----->', t_now, 'Night plans uploaded on the Cambridge server')
#
#         ## zip_files
#         path_database_zip_files = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
#                                      'schedule','Zip_files')
#         path_local_zip_file = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
#                                            'Zip_files/', str(t_now) + '.zip')
#         subprocess.Popen(["sshpass", "-p", pwd_appcs, "scp","-r",path_local_zip_file,path_database_zip_files])
#         print('----->', t_now, 'Zip Plans_by_dates folder uploaded on the Cambridge server')
#
#         #upload on data reduction computer
#         ## cam server to local
#         path_database_zip_file = os.path.join('speculoos@appcs.ra.phy.cam.ac.uk:/appct/data/SPECULOOSPipeline/', 'Ganymede',
#                                      'schedule','Zip_files',str(t_now)+'.zip')
#         path_local_zip_folder = os.path.join(path_spock + '/DATABASE/', 'Ganymede',
#                                   'Zip_files/')
#         p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_database_zip_file,path_local_zip_folder])
#         ## Local to reduction computer
#         p = subprocess.Popen(["sshpass", "-p", pwd_HUB,"scp", path_local_zip_file, 'speculoos@172.16.4.169:/home/speculoos/Plans_scheduler/Ganymede/Plans/'])
#         print('----->',t_now,'Zip Plans_by_dates folder uploaded on the HUB for Ganymede')
