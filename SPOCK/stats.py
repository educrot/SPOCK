from bs4 import BeautifulSoup
import requests
import pandas as pd
from astropy.time import Time
import numpy as np
from alive_progress import alive_bar
import os.path, time
import requests
import yaml

# ************************ Create database ************************

telescopes_names = ['Io', 'Europa', 'Ganymede', 'Callisto', 'Artemis', 'Saint-Ex', 'TS_La_Silla', 'TN_Oukaimeden']
if not os.path.exists('./target_lists'):
    os.makedirs('./target_lists')
if not os.path.exists('./survey_hours'):
    os.makedirs('./survey_hours')
if not os.path.exists('./DATABASE'):
    os.makedirs('./DATABASE')
if not os.path.exists('./night_blocks_propositions'):
    os.makedirs('./night_blocks_propositions')
if not os.path.exists('./SPOCK_files'):
    os.makedirs('./SPOCK_files')
for tel in telescopes_names:
    if not os.path.exists('./DATABASE/' + tel):
        os.makedirs('./DATABASE/' + tel)
for tel in telescopes_names:
    if not os.path.exists('./DATABASE/' + tel + '/Archive_night_blocks'):
        os.makedirs('./DATABASE/' + tel + '/Archive_night_blocks')
for tel in telescopes_names:
    if not os.path.exists('./DATABASE/' + tel + '/Plans_by_date'):
        os.makedirs('./DATABASE/' + tel + '/Plans_by_date')
for tel in telescopes_names:
    if not os.path.exists('./DATABASE/' + tel + '/Zip_files'):
        os.makedirs('./DATABASE/' + tel + '/Zip_files')

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

# ************************ Read target lists from server ************************
target_lists = ['speculoos_target_list_v6.txt', 'target_list_special.txt', 'target_transit_follow_up.txt']
for t_list in target_lists:
    target_list_url = "http://www.mrao.cam.ac.uk/SPECULOOS/spock_files/target_lists/" + t_list
    resp = requests.get(target_list_url, auth=(user_portal, pwd_portal))
    content = resp.text.replace("\n", "")
    open('./target_lists/' + t_list, 'wb').write(resp.content)

survey_hours = ['ObservationHours_Saint-Ex.txt', 'ObservationHours_TRAPPIST.txt', 'ObservationHours.txt']
for file in survey_hours:
    target_list_url = "http://www.mrao.cam.ac.uk/SPECULOOS/spock_files/survey_hours/" + file
    resp = requests.get(target_list_url, auth=(user_portal, pwd_portal))
    content = resp.text.replace("\n", "")
    open('./survey_hours/' + file, 'wb').write(resp.content)

# **********************************************************************************************************

target_list_v6 = pd.read_csv('./target_lists/speculoos_target_list_v6.txt',sep=' ')

def read_night_plans_server(telescope,date):
    TargetURL = "http://www.mrao.cam.ac.uk/SPECULOOS/"+telescope+\
                "/schedule/Archive_night_blocks/night_blocks_"+telescope+"_"+date+".txt"
    resp = requests.get(TargetURL, auth=(user_portal, pwd_portal))
    content = resp.text.replace("\n", "")
    open('text_file.txt', 'wb').write(resp.content)

    df =  pd.read_csv('text_file.txt', delimiter=' ', skipinitialspace=True, error_bad_lines=False)
    return df

def read_all_night_plans_server(file):
    TargetURL = file
    resp = requests.get(TargetURL, auth=(user_portal, pwd_portal))
    content = resp.text.replace("\n", "")
    open('text_file.txt', 'wb').write(resp.content)
    df = pd.read_csv('text_file.txt', delimiter=' ', skipinitialspace=True, error_bad_lines=False)
    return df

def listFD(url, ext=''):
    page = requests.get(url, auth=(user_portal, pwd_portal)).text
    # print(page)
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

def df_all_obs_scheduled(telescope):
    date_night_plan = []
    url = "http://www.mrao.cam.ac.uk/SPECULOOS/" + telescope + "/schedule/Archive_night_blocks/"
    ext = 'txt'

    with alive_bar(len(listFD(url, ext))) as bar:
        for i, file in enumerate(listFD(url, ext)):
            bar()
            time.sleep(0.001)
            if i == 0:
                df = read_all_night_plans_server(file)
            else:
                a = read_all_night_plans_server(file)

                frames = [df, a]
                df = pd.concat(frames)
                df = df.reset_index(drop=True)
                date_night_plan.append(
                    file.replace('http://www.mrao.cam.ac.uk/SPECULOOS/', '').replace(telescope, '').replace(
                        '/schedule/Archive_night_blocks//night_blocks_', '').replace('_', '').replace('.txt', ''))
        return df,date_night_plan

def date_night_start_func(df_speculoos,target):
    idx_target = np.where((df_speculoos['target'] == target))[0]
    ici = np.where((target_list_v6['Sp_ID'] == target))[0]
    if len(idx_target) != 0:
        date_oldest = min(df_speculoos['start time (UTC)'][idx_target])
        date_most_recent = max(df_speculoos['start time (UTC)'][idx_target])
        durations = np.array(df_speculoos['duration (minutes)'][idx_target])
        date_night_start  = []
        for i in range(len(idx_target)):
            date_night_start.append(Time(df_speculoos['start time (UTC)'][idx_target[i]], out_subfmt='date').iso)
    else:
        date_night_start = []
        date_oldest = 'None'
        date_most_recent = 'None'
        durations = []



    return date_night_start,date_oldest,date_most_recent,durations


def run_masterfile():
    df_all_Io,date_night_plan_Io = df_all_obs_scheduled('Io')
    df_all_Europa,date_night_plan_Europa = df_all_obs_scheduled('Europa')
    df_all_Ganymede,date_night_plan_Ganymede = df_all_obs_scheduled('Ganymede')
    df_all_Callisto,date_night_plan_Callisto = df_all_obs_scheduled('Callisto')
    df_all_Artemis,date_night_plan_Artemis = df_all_obs_scheduled('Artemis')
    df_all_TS_La_Silla,date_night_plan_TS_La_Silla = df_all_obs_scheduled('TS_La_Silla')
    df_all_TN_Oukaimeden,date_night_plan_TN_Oukaimeden = df_all_obs_scheduled('TN_Oukaimeden')
    df_all_TN_Oukaimeden,date_night_plan_TN_Oukaimeden = df_all_obs_scheduled('Saint-Ex')

    frames = [df_all_Io,df_all_Europa,df_all_Ganymede,
              df_all_Callisto,df_all_Artemis,df_all_TS_La_Silla,df_all_TN_Oukaimeden]

    df_speculoos = pd.concat(frames)
    df_speculoos = df_speculoos.sort_values('target').reset_index(drop=True)
    df_speculoos.to_csv('/Users/elsaducrot/spock_2/SPOCK_files/all_schedules.csv',sep=',',index=None)

    date_night_start_each_target = []
    date_oldest_obs = []
    date_most_recent_obs = []
    durations_all_obs = []

    idx_all = np.where((np.array([target_list_v6['telescope'][i].find('[]')
                                  for i in range(len(target_list_v6))]) == -1))[0]

    for target in target_list_v6['Sp_ID'][idx_all]:
        try:
            basic_info = date_night_start_func(df_speculoos,target)
        except UnboundLocalError:
            print('ERROR: solve')
        date_night_start_each_target.append(basic_info[0])
        date_oldest_obs.append(basic_info[1])
        date_most_recent_obs.append(basic_info[2])
        durations_all_obs.append(basic_info[3])

    df_masterfile = pd.DataFrame({'Sp_ID': target_list_v6['Sp_ID'][idx_all], 'RA': target_list_v6['RA'][idx_all],
                       'DEC': target_list_v6['DEC'][idx_all], 'telescope': target_list_v6['telescope'][idx_all],
                       'Program': target_list_v6['Program'][idx_all],
                       'nb_hours_surved': target_list_v6['nb_hours_surved'][idx_all],
                       'all_dates_scheduled':date_night_start_each_target,
                       'all_durations_scheduled':durations_all_obs,
                       'oldest_obs':date_oldest_obs,'most_recent_obs':date_most_recent_obs,
                       'Filter_spc': target_list_v6['Filter_spc'][idx_all],
                       'texp_spc': target_list_v6['texp_spc'][idx_all],
                       'Ms': target_list_v6['Ms'][idx_all], 'Rs': target_list_v6['Rs'][idx_all],
                       'SpT': target_list_v6['SpT'][idx_all]})

    df_masterfile.to_csv('/Users/elsaducrot/spock_2/SPOCK_files/spock_stats_masterfile.csv',sep=',',index=None)