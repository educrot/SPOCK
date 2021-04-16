__all__ = ['long_term_scheduler','short_term_scheduler','make_night_plans','plots_scheduler',
           'txt_files','upload_night_plans','stats','SPECULOOSScheduler','pwd_appcs','pwd_HUB','user_portal',
           'pwd_portal','pwd_appcs','pwd_SNO_Reduc1','user_chart_studio','pwd_chart_studio','path_spock',
           'path_credential_json']

__version__ = "0.0.1"

import pkg_resources
import os
import requests
import yaml
from colorama import Fore


def _get_files():
    data_path = pkg_resources.resource_filename('SPOCK', 'credentials/')
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

        # ************************ Create database ************************

        telescopes_names = ['Io', 'Europa', 'Ganymede', 'Callisto', 'Artemis', 'Saint-Ex',
                            'TS_La_Silla', 'TN_Oukaimeden']
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

        survey_hours = ['ObservationHours_Saint-Ex.txt', 'ObservationHours_TRAPPIST.txt',
                        'ObservationHours.txt','SurveyTotal.txt']
        for file in survey_hours:
            target_list_url = "http://www.mrao.cam.ac.uk/SPECULOOS/spock_files/survey_hours/" + file
            resp = requests.get(target_list_url, auth=(user_portal, pwd_portal))
            content = resp.text.replace("\n", "")
            open(path_spock + '/survey_hours/' + file, 'wb').write(resp.content)

        return pwd_appcs, pwd_HUB, user_portal, pwd_portal, pwd_appcs, pwd_SNO_Reduc1, user_chart_studio,\
               pwd_chart_studio, path_spock,path_credential_json

        # **********************************************************************************************************
    else:
        print(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' No file '+ 'passwords.csv')


pwd_appcs,pwd_HUB, user_portal, pwd_portal, pwd_appcs, pwd_SNO_Reduc1, user_chart_studio, pwd_chart_studio, path_spock, path_credential_json = _get_files()


from .long_term_scheduler import *
from .short_term_scheduler import *
from .make_night_plans import *
from .plots_scheduler import *
from .txt_files import *
from .upload_night_plans import *
from .stats import *