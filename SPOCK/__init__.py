__all__ = ['long_term_scheduler','short_term_scheduler','make_night_plans','plots_scheduler',
           'txt_files','upload_night_plans','stats','SPECULOOSScheduler','pwd_appcs','pwd_HUB','user_portal',
           'pwd_portal','pwd_appcs','pwd_SNO_Reduc1','user_chart_studio','pwd_chart_studio','path_spock',
           'path_credential_json','login_stargate','pwd_stargate']

__version__ = "0.0.1"

import pkg_resources
import os
import requests
import yaml
from ftplib import FTP
from colorama import Fore
from datetime import date, timedelta, datetime
import pandas as pd
import sys
import numpy as np

def index_list1_list2(list1, list2):  # list 2 longer than list 1
    """ index of list1 in list2 and list2 in list1

    Parameters
    ----------
    list1 : list

    list2 : list

    Returns
    -------
    list
        list of index of list1 in list2 and list2 in list1

    """
    idx_list1_in_list2 = []
    idx_list2_in_list1 = []
    for i in range(len(list2)):
        for j in range(len(list1)):
            if list2[i] == list1[j]:
                idx_list1_in_list2.append(i)
                idx_list2_in_list1.append(j)
    return idx_list1_in_list2, idx_list2_in_list1

def _get_files():
    data_path = pkg_resources.resource_filename('SPOCK', 'credentials/')
    if not os.path.exists(data_path):
        os.makedirs(data_path)
    filename_pwd = os.path.join(data_path, 'passwords.csv')
    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' Please add password.csv file in: ' + data_path)
    if os.path.exists(filename_pwd):
        print(Fore.GREEN + 'INFO: ' + Fore.BLACK + ' OK Password file exists')
        # ************************ Read passwords ************************

        with open(filename_pwd, "r") as f:
            Inputs = yaml.load(f, Loader=yaml.FullLoader)
            pwd_appcs = Inputs['pwd_appcs'][0]
            pwd_HUB = Inputs['pwd_HUB'][0]
            user_portal = Inputs['user_portal'][0]
            pwd_portal = Inputs['pwd_portal'][0]
            pwd_appcs = Inputs['pwd_appcs'][0]
            pwd_SNO_Reduc1 = Inputs['pwd_SNO_Reduc1'][0]
            user_chart_studio = Inputs['user_chart_studio'][0]
            pwd_chart_studio = Inputs['pwd_chart_studio'][0]
            path_spock = Inputs['path_spock'][0]
            path_credential_json = Inputs['credential_json'][0]
            login_stargate = Inputs['login_stargate'][0]
            pwd_stargate = Inputs['pwd_stargate'][0]

        # ************************ Create database ************************

        telescopes_names = ['Io', 'Europa', 'Ganymede', 'Callisto', 'Artemis', 'Saint-Ex',
                            'TS_La_Silla', 'TN_Oukaimeden']
        if not os.path.exists(path_spock + '/target_lists'):
            os.makedirs(path_spock + '/target_lists')
        if not os.path.exists(path_spock + '/target_lists/stargate'):
            os.makedirs(path_spock + '/target_lists/stargate')
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

        # ********* Read target lists from server ******** UNCOMMENT IF NOT USING STARGATE  AND WG6 ***
        # target_lists = ['speculoos_target_list_v6.txt', 'target_list_special.txt', 'target_transit_follow_up.txt']
        # for t_list in target_lists:
        #     target_list_url = "http://www.mrao.cam.ac.uk/SPECULOOS/spock_files/target_lists/" + t_list
        #     resp = requests.get(target_list_url, auth=(user_portal, pwd_portal))
        #     content = resp.text.replace("\n", "")
        #     open(path_spock + '/target_lists/' + t_list, 'wb').write(resp.content)
        #
        # survey_hours = ['ObservationHours_Saint-Ex.txt', 'ObservationHours_TRAPPIST.txt',
        #                 'ObservationHours.txt','SurveyTotal.txt']
        # for file in survey_hours:
        #     target_list_url = "http://www.mrao.cam.ac.uk/SPECULOOS/spock_files/survey_hours/" + file
        #     resp = requests.get(target_list_url, auth=(user_portal, pwd_portal))
        #     content = resp.text.replace("\n", "")
        #     open(path_spock + '/survey_hours/' + file, 'wb').write(resp.content)

        return pwd_appcs, pwd_HUB, user_portal, pwd_portal, pwd_appcs, pwd_SNO_Reduc1, user_chart_studio,\
               pwd_chart_studio, path_spock, path_credential_json, login_stargate, pwd_stargate

        # **********************************************************************************************************
    else:
        print(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' No file '+ 'passwords.csv')


def get_target_list_stargate(day):
    """

    Parameters
    ----------
    day

    Returns
    -------

    """
    objdate = datetime.strptime(day, '%Y-%m-%d') #- timedelta(days=1)
    y = datetime.strftime(objdate,'%Y')
    m = datetime.strftime(objdate,'%m').lstrip('0')
    d = datetime.strftime(objdate,'%d').lstrip('0')
    ftp = FTP('z93vm.ftp.infomaniak.com')
    ftp.login(login_stargate,pwd_stargate)
    # day = Time(day, scale='utc', out_subfmt='date').iso
    file_name = 'stargate_db_'+y+'-'+m+'-'+d+'.csv'
    my_file = open(path_spock + '/target_lists/stargate/' + file_name, 'wb')
    ftp.retrbinary('RETR ' + file_name, my_file.write,8*1024)
    print(Fore.GREEN + 'INFO: ' + Fore.BLACK + 'Downloading target list from STARGATE.')
    # time.sleep(5)  # wait for the file to fully download
    return file_name


def change_fmt_stargate_TL(file_name):
    """

    Parameters
    ----------
    file_name

    Returns
    -------

    """
    # change columns names
    df = pd.read_csv(path_spock + '/target_lists/stargate/' + file_name, delimiter=';')
    df = df.rename(columns={'spc': 'Sp_ID', 'ra': 'RA', 'dc': 'DEC', 'gaia': 'Gaia_ID', 'obstime': 'nb_hours_surved',
                            'program.1': 'Program', 'spt': 'SpT', 'mag_j': 'J'})
    df['nb_hours_surved'] = df['nb_hours_surved'].fillna(0)
    df['texp_spc'] = [0] * len(df['Sp_ID'])
    df['Filter_spc'] = ['I+z'] * len(df['Sp_ID'])
    df['texp_trap'] = [0] * len(df['Sp_ID'])
    df['Filter_trap'] = ['I+z'] * len(df['Sp_ID'])
    df['telescope'] = [''] * len(df['Sp_ID'])
    df['nb_hours_threshold'] = [100] * len(df['Sp_ID'])
    df.loc[df.Program == 1, 'nb_hours_threshold'] = 200
    df['SNR_TESS_temp'] = [0]*len(df['SNR_JWST_HZ_tr'])
    df['SNR_Spec_temp'] = [0]*len(df['SNR_JWST_HZ_tr'])
    df['SNR_JWST_HZ_tr'] = df['SNR_JWST_HZ_tr'].fillna(0)
    df['SNR_TESS_temp'] = df['SNR_TESS_temp'].fillna(0)
    df['SNR_Spec_temp'] = df['SNR_Spec_temp'].fillna(0)
    df["SNR_SPIRIT"] = [0]*len(df)
    df["texp_spirit"] = [0]*len(df)
    df_spirit = pd.read_csv("/Users/ed268546/Documents/codes/SPOCK/SPIRIT/target_precision_df_1.2seeing_andorSPC_-60_I+z_pirtSPC_-60_real_zYJ_final_upgraded_2022-05-17T120501.csv", sep=',')
    idx_list1_in_list2, idx_list2_in_list1 = index_list1_list2(df_spirit["Sp_ID"], df["Sp_ID"])
    df["SNR_SPIRIT"][idx_list1_in_list2] = df_spirit["SNR_1"][idx_list2_in_list1]
    df["texp_spirit"][idx_list1_in_list2] = df_spirit["exp_time_2"][idx_list2_in_list1]

    f = open(path_spock + '/target_lists//www.mrao.cam.ac.uk/SPECULOOS/speculoos-portal/php/get_hours.php', 'r')
    f = f.read()
    line = f#.strip()
    columns = line.split('","')

    names = []
    hours = []

    for c in columns:
        if c.find('SP') != -1 and c.find('TESS') == -1 and c.find('WAS') == -1:
            info = c.split(',')
            names.append(info[0].replace(' ','').replace('SP','Sp'))
            hours.append(float(info[1]))

    df_portal = pd.DataFrame({"Sp_ID":names,"nb_hours_surved": hours})

    idx_list1_in_list2, idx_list2_in_list1 = index_list1_list2(df["Sp_ID"],df_portal["Sp_ID"])

    df["nb_hours_surved"][idx_list2_in_list1] = df_portal["nb_hours_surved"][idx_list1_in_list2]
    idx_double_stars = np.where((df["Sp_ID"] == "Sp1633-6808_2") | (df["Sp_ID"] == "Sp1633-6808_1") |
                             (df["Sp_ID"] == "Sp1953+4424_1") | (df["Sp_ID"] == "Sp1953+4424_2") |
                           (df["Sp_ID"] == "Sp0933-4353_1") | (df["Sp_ID"] == "Sp0933-4353_2"))[0]
    df["nb_hours_surved"][idx_double_stars] = 200

    df.to_csv(path_spock + '/target_lists/stargate/' + 'TL_spock_' + file_name, sep=',', index=None)

    return path_spock + '/target_lists/stargate/' + 'TL_spock_' + file_name

pwd_appcs,pwd_HUB, user_portal, pwd_portal, pwd_appcs, pwd_SNO_Reduc1, user_chart_studio, pwd_chart_studio, path_spock, path_credential_json, login_stargate, pwd_stargate = _get_files()

today = date.today() - timedelta(days=1)
today = today.strftime("%Y-%m-%d")
target_list_from_stargate_path = change_fmt_stargate_TL(get_target_list_stargate(today))


from .long_term_scheduler import *
from .short_term_scheduler import *
from .make_night_plans import *
from .plots_scheduler import *
from .txt_files import *
from .upload_night_plans import *
from .stats import *
