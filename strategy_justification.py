# phase coverage simulation
import numpy as np
#from specphot.observations import *
#from specphot.lightcurve import *
#from specphot import spio
import matplotlib.pyplot as plt
#from specphot import utils
import numpy as np
import pandas as pd
import requests
import os

from os import path
import numpy as np
from tqdm import tqdm
import sys
import shutil



def phase_coverage(period, times, midpoint=0, binning=0.005):
    #if len(_times) == 1:
    _times = np.array(
        [
            ((time - (np.round((time - midpoint) / period) * period)) - midpoint)
            / period
            for time in times
        ]
    )
    # _times = [_times, _times]
    try:
        _stacked_times = np.hstack(_times)
        bin_points = len(
            utils.binning_for_cov(np.ones(len(_stacked_times)), _stacked_times, binning / period)[0]
        )
        total_bin_points = len(np.arange(-0.5, 0.5, binning / period))
        return bin_points / total_bin_points
    except ValueError:
        print('WARNING: No data for this target')

def plot_polar_phase_coverage(midpoint, period, _times, r,d):
    R = 1*d #5*(period/(2*np.pi))

    times = np.array([np.array([np.min(t), np.max(t)]) for t in _times])

    length = [(np.max(t) - np.min(t))/period for t in times]
    offseted_times =  np.array(
        [
            ((time - (np.round((time - midpoint) / period) * period)) - midpoint)
            / period
            for time in times
        ]
    )
    #fig, ax1 = plt.subplots(figsize=(20,15))
    plt.pie([1], colors=['#016CA000'], radius=r)

    for l, t in zip(length, offseted_times):
        plot_slice(t[0], l, R=R)

    my_circle=plt.Circle( (0,0), R - 0.3, color='white')
    planet = plt.Circle((-np.cos(2*np.pi*(-midpoint/period) + np.pi/2)*(R-.15),(R-.15)*np.sin(2*np.pi*(midpoint/period) + np.pi/2)), 0.1, facecolor='white', edgecolor='black')
    star = plt.Circle((0,0), 0.5, facecolor='white', edgecolor='black')
    p=plt.gcf()
    p.gca().add_artist(my_circle)
    p.gca().add_artist(star)
    p.gca().add_artist(planet)
    #plt.legend([my_circle],['For Period ' + str(round(period,2)) + ' days'],loc=4, fontsize = 'x-large')
    plt.autoscale()
    return utils.phase_coverage(np.abs(period), _times)

def plot_slice(offset, width, color='#016CA050', R=3):
    size = 1
    plt.pie([width*100, 100-width*100], colors=[color, [0,0,0,0]], counterclock=False, startangle=90-offset*360, radius=R,\
            wedgeprops=dict(width=size))

def cambridge_download_product(destination,target="",date="",telescope="",process_lc=False,user='educrot',\
                               password="58JMSGgdmzTB",target_list_path='speculoos_target_list_v3.txt',\
                               target_gaia_id=None,stack="global",):
    """
    Download photometry products from cambridge archive in the following file tree:

    .. code-block::

        Target/
         │
         └── Telescope_date_target_filter/
              ├── Telescope_date_target_filter_photometry.phots
              └── Telescope_date_target_filter_stack.fits

    Parameters
    ----------
    destination: string (path)
        local destination path
    target : string
        target name in the cambridge archive
    date : string
        date with format YYYYmmdd (e.g. 20190130)
    telescope : string
        telescope name in the cambridge archive
    process_lc: bool
        wheather to process the lighturve and store it in the ``phots``file
    user : string
        username for the cambridge archive
    password : string
        password for the cambridge archive
    target_list_path : string path
        local path of the target list file (csv)
    target_gaia_id : string
        gaia id of the target (optional, will be found from the target list)

    """
    destination = path.abspath(destination)
    if not path.exists(destination):
        os.mkdir(destination)

    url = "http://www.mrao.cam.ac.uk/SPECULOOS/portal/get_file.php?telescope={}&date={}&id={}&filter=&file=MCMC_*".format(
        telescope, date, target
    )
    resp = requests.get(url, auth=(user, password))
    assert (
        resp.status_code == 200
    ), "Wrong username or password used to access data, please check .specphot.config file"
    assert (
        resp.content != b"null" and resp.content != b"\r\nnull"
    ), "Your request is not matching any available data in the Cambridge archive. To see available data, please check http://www.mrao.cam.ac.uk/SPECULOOS/portal/"
    #output_fits_urls = np.array(
    #    [
    #        ("http://www.mrao.cam.ac.uk/SPECULOOS/" + ur[4::]).replace("\\", "")
    #        for ur in eval(resp.content)
    #    ]
    #)
    output_urls = np.array(
        [
            ("http://www.mrao.cam.ac.uk/SPECULOOS/" + ur[4::]).replace("\\", "")
            for ur in eval(resp.content)
        ]
    )
    for url in tqdm(output_urls):
        target_date_folder = path.join(destination,
                                       '{}_{}_{}_cambridge'.format(*url.split("/")[6:8][::-1], url.split("/")[4]))
        if not path.exists(target_date_folder):
            os.mkdir(target_date_folder)

        resp = requests.get(url, auth=(user, password))
        _bytes = resp.content
        open(path.join(target_date_folder, "{}.txt".format(url.split("/")[-1])), 'wb').write(_bytes)

def download_MCMC_files_from_Cambridge(destination, target="", telescope="", date="", user=None, password=None):
    """
    Download MCMC files from cam archive in the following structure:

        destination/
            ├── target_date_telescope/
            │    ├── MCMC_0aperture_1/txt
            │    ├── MCMC_aperture_2.txt
            │    └── ...
            │
            └── target2_date2_telescope2/
                 ├── MCMC_0aperture_1/txt
                 ├── MCMC_aperture_2.txt
                 └── ...

    when telescope (resp. date and target) are not specified, all telescope (resp. date and target) will be take into account

    Parameters
    ----------
    destination: string
        path of the destination folder (will be created if not existing)
    target: str
        name of the target
    telescope:
        telescope name
    date:
        date with format "yyyymmdd" (example: 20180520 for 20/05/2018)
    user:
        username for the cambridge archive
    password:
        password for the cambridge archive

    Returns
    -------

    """
    if not path.exists(destination):
        os.mkdir(destination)

    url = "http://www.mrao.cam.ac.uk/SPECULOOS/portal/get_file.php?telescope={}&date={}&id={}&filter=&file=MCMC_*".format(
        telescope, date, target
    )
    resp = requests.get(url, auth=(user, password))
    assert (
            resp.content != b"null"
    ), "Your request is not matching any available data in the Cambridge archive. To see available data, please check http://www.mrao.cam.ac.uk/SPECULOOS/portal/"

    assert (
            resp.status_code == 200
    ), "Wrong username or password used to access data, please check .specphot.config file"
    output_urls = np.array(
        [
            ("http://www.mrao.cam.ac.uk/SPECULOOS/" + ur[4::]).replace("\\", "")
            for ur in eval(resp.content)
        ]
    )

    for url in tqdm(output_urls):
        target_date_folder = path.join(destination,
                                       '{}_{}_{}_cambridge'.format(*url.split("/")[6:8][::-1], url.split("/")[4]))
        if not path.exists(target_date_folder):
            os.mkdir(target_date_folder)

        resp = requests.get(url, auth=(user, password))
        _bytes = resp.content
        open(path.join(target_date_folder, "{}.txt".format(url.split("/")[-1])), 'wb').write(_bytes)

class coverage:
    def __init__(self):
        self.name_target = None
        self.files_concat = None
        self.times = None

    def read_LCs(self):
        downloaded = []
        for i in os.listdir('/Users/elsaducrot/specphot/coverage/'):
            if i.startswith(self.name_target):
                downloaded.append(True)
                print('WARNING: folder already exists')

            else:
                downloaded.append(False)

 #       if not any(downloaded):
        cambridge_download_product('/Users/elsaducrot/specphot/coverage/',target="",date="",telescope="",process_lc=False,user="educrot",password="58JMSGgdmzTB",target_list_path='speculoos_target_list_v3.txt',target_gaia_id=None,stack="global")
        for j in os.listdir('/Users/elsaducrot/specphot/coverage/'):
            if j.startswith(self.name_target):
                shutil.move("/Users/elsaducrot/specphot/coverage/" + j, '/Users/elsaducrot/specphot/coverage' + '/' + self.name_target + '/')

        return any(downloaded)

    def upload_data(self):
        files=[]
        self.times = []
        mcmc_alone = False
        for i in os.listdir('/Users/elsaducrot/specphot/coverage/' + self.name_target):
            if (i.startswith(self.name_target))  or (i.startswith('SPEC'))  or (i.startswith('SPC')):
                try:
                    files.append(pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/'  +  i + '/' +  'MCMC_text_6.txt',sep=' ',skipinitialspace=True))
                    df = pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i  + '/' +  'MCMC_text_6.txt',sep=' ',skipinitialspace=True)
                except FileNotFoundError:
                    try:
                        files.append(pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i + '/' +  'MCMC_text_3.txt',sep=' ',skipinitialspace=True))
                        df = pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i + '/' +  'MCMC_text_3.txt',sep=' ',skipinitialspace=True)
                    except FileNotFoundError:
                        try:
                            files.append(pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i + '/' +  'MCMC_text_5.txt',sep=' ',skipinitialspace=True))
                            df = pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i + '/' +  'MCMC_text_5.txt',sep=' ',skipinitialspace=True)
                        except FileNotFoundError:
                            try:
                                files.append(pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i + '/' +  'MCMC_text_4.txt',sep=' ',skipinitialspace=True))
                                df = pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i + '/' +  'MCMC_text_4.txt',sep=' ',skipinitialspace=True)
                            except FileNotFoundError:
                                try:
                                    files.append(pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i + '/' +  'MCMC_text.txt',sep=' ',skipinitialspace=True))
                                    df = pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i + '/' +  'MCMC_text.txt',sep=' ',skipinitialspace=True)
                                except FileNotFoundError:
                                    files.append(pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i + '/' +  'MCMC_text_7.txt',sep=' ',skipinitialspace=True))
                                    df = pd.read_csv('/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i + '/' +  'MCMC_text_7.txt',sep=' ',skipinitialspace=True)

                self.times.append(df['TMID-2450000'] + 2450000)
                self.files_concat = pd.concat(files,axis =0,ignore_index=True,sort=True)
                self.files_concat = files

            if i.startswith('MCMC_'):
                if not mcmc_alone:
                    try:
                        files.append(pd.read_csv(
                            '/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i ,
                            sep=' ', skipinitialspace=True))
                        df = pd.read_csv(
                            '/Users/elsaducrot/specphot/coverage/' + self.name_target + '/' + i ,
                            sep=' ', skipinitialspace=True)
                        mcmc_alone = True
                        self.times.append(df['TMID-2450000'] + 2450000)
                        self.files_concat = pd.concat(files, axis=0, ignore_index=True, sort=True)
                        self.files_concat = files
                        print(mcmc_alone)
                    except FileNotFoundError:
                        print('ok')
                else:
                    print('nop')



    def observed_time(self):
        observed_time = np.array([t + 2450000  for t in self.files_concat['TMID-2450000']])
        return observed_time
