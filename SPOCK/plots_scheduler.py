from astroplan import Observer,FixedTarget
from astroplan.plots import plot_airmass
from astroplan.utils import time_grid_from_range
from astropy.coordinates import SkyCoord, EarthLocation, get_moon
from astropy.time import Time
from astropy import units as u
import chart_studio
from colorama import Fore
from datetime import date, timedelta
import io
import json
from json.decoder import JSONDecodeError
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
import os
import pandas as pd
from plotly import graph_objs as go
import plotly.figure_factory as ff
from plotly import offline
import requests
from requests.auth import HTTPBasicAuth
import sys
from SPOCK import user_portal, pwd_portal, user_chart_studio, pwd_chart_studio, \
    path_spock, target_list_from_stargate_path
import SPOCK.short_term_scheduler as SPOCKST

chart_studio.tools.set_credentials_file(username=user_chart_studio, api_key=pwd_chart_studio)
chart_studio.tools.set_config_file(world_readable=True, sharing='public')


def charge_observatories(Name):
    """ charge the observatory

    Parameters
    ----------
    Name : str
        Name of the observatory

    Returns
    -------
    list
        list of the observatory
    """

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
        location_TNOuka = EarthLocation.from_geodetic( -7.862263*u.deg,31.20516*u.deg, 2751*u.m)
        observatories.append(Observer(location=location_TNOuka, name="TNOuka", timezone="UTC"))

    if 'Munich' in str(Name):
        location_munich= EarthLocation.from_geodetic(48.2*u.deg, -11.6*u.deg, 600*u.m)
        observatories.append(Observer(location=location_munich, name="Munich", timezone="UTC"))

    return observatories

def airmass_plot_saved(name_observatory,telescope,day):
    """

    Parameters
    ----------
    name_observatory : str
        name of the observatory (ex : 'SSO')
    telescope : str
        name of the telescope (exx : 'Io')
    day : date
        date of day in fmt "yyyy-mm-dd"


    Returns
    -------
    plot
         visibility plot on telescope at a given day

    """
    night_block = pd.read_csv(os.path.join(path_spock + '/DATABASE/', telescope,
                                              'Archive_night_blocks','night_blocks_' + telescope + '_' + day.tt.datetime.strftime("%Y-%m-%d") + '.txt'),\
                              sep=' ', skipinitialspace=True)
    observatory = charge_observatories(name_observatory)[0]
    day = Time(night_block['start time (UTC)'][0]) - 5 *u.hour
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day,which='next').jd - 0.07, \
                                      observatory.twilight_morning_nautical(day+1,which='nearest').jd + 0.07, 100),format='jd')
    colors_start_new_target = ['black', 'darkgray', 'lightgray']
    sun_set =observatory.twilight_evening_nautical(day,which='next').iso
    sun_rise = observatory.twilight_morning_nautical(day+1,which='nearest').iso
    for i in range(len(night_block)):

        dec = str(int(float(night_block['dec (d)'][i]))) + ' ' + str(
            int(abs(float(night_block['dec (m)'][i])))) + ' ' + str(int(abs(float(night_block['dec (s)'][i]))))
        ra = str(int(float(night_block['ra (h)'][i]))) + ' ' + str(
            int(abs(float(night_block['ra (m)'][i])))) + ' ' + str(int(abs(float(night_block['ra (s)'][i]))))
        plot_airmass(FixedTarget(coord=SkyCoord(ra=ra,dec=dec,unit=(u.hourangle, u.deg)),\
                                 name=night_block['target'][i]), observatory, delta_midnight)
        t = Time(night_block['start time (UTC)'][i])
        plt.vlines(t.iso, 3, 1,linestyle='--',color=colors_start_new_target[i],alpha = 0.8,label='start ' + str(night_block['target'][i]))
        plt.vlines(sun_set, 3, 1, linestyle=':', color='yellow', alpha=0.9)
        plt.vlines(sun_rise, 3, 1, linestyle=':', color='yellow', alpha=0.9)
        plt.grid(color='gainsboro', linestyle='-', linewidth=1,alpha=0.3)
        plt.title('Visibility plot for the night of the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' on ' + str(telescope))
        plt.legend(loc=2)

def airmass_plot_proposition(name_observatory,telescope,day):
    """

    Parameters
    ----------
    name_observatory : str
        name of the observatory (ex : 'SSO')
    telescope : str
        name of the telescope (exx : 'Io')
    day : date
        date of day in fmt "yyyy-mm-dd"


    Returns
    -------
    plot
         visibility plot on telescope at a given day

    """
    night_block = pd.read_csv(os.path.join(path_spock + 'night_blocks_propositions/','night_blocks_' + telescope + '_' + day.tt.datetime.strftime("%Y-%m-%d") + '.txt'),\
                              sep=' ', skipinitialspace=True)
    observatory = charge_observatories(name_observatory)[0]
    day = Time(night_block['start time (UTC)'][0]) - 5 *u.hour
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day,which='next').jd - 0.07,\
                                      observatory.twilight_morning_nautical(day+1,which='nearest').jd + 0.07, 100),format='jd')
    colors_start_new_target = ['black', 'darkgray', 'lightgray']
    sun_set =observatory.twilight_evening_nautical(day,which='next').iso
    sun_rise = observatory.twilight_morning_nautical(day+1,which='nearest').iso
    for i in range(len(night_block)):

        dec = str(int(float(night_block['dec (d)'][i]))) + ' ' + str(
            int(abs(float(night_block['dec (m)'][i])))) + ' ' + str(int(abs(float(night_block['dec (s)'][i]))))
        ra = str(int(float(night_block['ra (h)'][i]))) + ' ' + str(
            int(abs(float(night_block['ra (m)'][i])))) + ' ' + str(int(abs(float(night_block['ra (s)'][i]))))
        plot_airmass(FixedTarget(coord=SkyCoord(ra=ra,dec=dec,unit=(u.hourangle, u.deg)),\
                                 name=night_block['target'][i]), observatory, delta_midnight)
        t = Time(night_block['start time (UTC)'][i])
        plt.vlines(t.iso, 3, 1,linestyle='--',color=colors_start_new_target[i],alpha = 0.8,label='start ' + str(night_block['target'][i]))
        plt.vlines(sun_set, 3, 1, linestyle=':', color='yellow', alpha=0.9)
        plt.vlines(sun_rise, 3, 1, linestyle=':', color='yellow', alpha=0.9)
        plt.grid(color='gainsboro', linestyle='-', linewidth=1,alpha=0.3)
        plt.title('Visibility plot for the night of the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' on ' + str(telescope))
        plt.legend(shadow=True, loc=2)

def airmass_altitude_plot_saved(name_observatory,telescope,day):
    """

    Parameters
    ----------
    name_observatory : str
        name of the observatory (ex : 'SSO')
    telescope : str
        name of the telescope (exx : 'Io')
    day : date
        date of day in fmt "yyyy-mm-dd"


    Returns
    -------
    plot
         visibility plot on telescope at a given day

    """
    night_block = pd.read_csv(os.path.join(path_spock + '/DATABASE/', telescope,
                                              'Archive_night_blocks','night_blocks_' + telescope + '_' + day.tt.datetime.strftime("%Y-%m-%d") + '.txt'),\
                              sep=' ', skipinitialspace=True)
    observatory = charge_observatories(name_observatory)[0]
    day = Time(night_block['start time (UTC)'][0]) - 5 *u.hour
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day,which='next').jd - 0.07,\
                                      observatory.twilight_morning_nautical(day+1,which='nearest').jd + 0.07, 100),format='jd')
    sun_set =observatory.twilight_evening_nautical(day,which='next').iso
    sun_rise = observatory.twilight_morning_nautical(day+1,which='nearest').iso
    fig, axs = plt.subplots(1)
    colors_start_new_target = ['black', 'darkgray', 'lightgray']
    for i in range(len(night_block)):
        dec = str(int(float(night_block['dec (d)'][i]))) + ' ' + str(
            int(abs(float(night_block['dec (m)'][i])))) + ' ' + str(int(abs(float(night_block['dec (s)'][i]))))
        ra = str(int(float(night_block['ra (h)'][i]))) + ' ' + str(
            int(abs(float(night_block['ra (m)'][i])))) + ' ' + str(int(abs(float(night_block['ra (s)'][i]))))
        plot_airmass(FixedTarget(coord=SkyCoord(ra=ra,dec=dec,unit=(u.hourangle, u.deg)),\
                                 name=night_block['target'][i]), observatory, delta_midnight, altitude_yaxis=True)
        plt.ylabel('Altitude (degrees)')
        t = Time(night_block['start time (UTC)'][i])
        axs.vlines(t.iso, 3, 1,linestyle='--',color=colors_start_new_target[i],alpha = 0.7,label='start ' + str(night_block['target'][i]))
        axs.vlines(sun_set, 3, 1, linestyle=':', color='yellow', alpha=0.9)
        axs.vlines(sun_rise, 3, 1, linestyle=':', color='yellow', alpha=0.9)
        #plt.legend(loc=2)
        plt.grid(color='gainsboro', linestyle='-', linewidth=1, alpha=0.3)
        plt.title('Visibility plot for the night of the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' on ' + str(telescope))

def airmass_altitude_plot_proposition(name_observatory,telescope,day):
    """

    Parameters
    ----------
    name_observatory : str
        name of the observatory (ex : 'SSO')
    telescope : str
        name of the telescope (exx : 'Io')
    day : date
        date of day in fmt "yyyy-mm-dd"


    Returns
    -------
    plot
         visibility plot on telescope at a given day

    """
    night_block = pd.read_csv(os.path.join(path_spock + '/night_blocks_propositions/','night_blocks_' + telescope + '_' + day.tt.datetime.strftime("%Y-%m-%d") + '.txt'),\
                              sep=' ', skipinitialspace=True)
    observatory = charge_observatories(name_observatory)[0]
    day = Time(night_block['start time (UTC)'][0]) - 5 *u.hour
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day,which='next').jd - 0.07,\
                                      observatory.twilight_morning_nautical(day+1,which='nearest').jd + 0.07, 100),format='jd')
    sun_set =observatory.twilight_evening_nautical(day,which='next').iso
    sun_rise = observatory.twilight_morning_nautical(day+1,which='nearest').iso
    fig, axs = plt.subplots(1)
    colors_start_new_target = ['black', 'darkgray', 'lightgray']
    for i in range(len(night_block)):
        dec = str(int(float(night_block['dec (d)'][i]))) + ' ' + str(
            int(abs(float(night_block['dec (m)'][i])))) + ' ' + str(int(abs(float(night_block['dec (s)'][i]))))
        ra = str(int(float(night_block['ra (h)'][i]))) + ' ' + str(
            int(abs(float(night_block['ra (m)'][i])))) + ' ' + str(int(abs(float(night_block['ra (s)'][i]))))
        plot_airmass(FixedTarget(coord=SkyCoord(ra=ra,dec=dec,unit=(u.hourangle, u.deg)),\
                                 name=night_block['target'][i]), observatory, delta_midnight, altitude_yaxis=True)
        plt.ylabel('Altitude (degrees)')
        t = Time(night_block['start time (UTC)'][i])
        axs.vlines(t.iso, 3, 1,linestyle='--',color=colors_start_new_target[i],alpha = 0.7,label='start ' + str(night_block['target'][i]))
        axs.vlines(sun_set, 3, 1, linestyle=':', color='yellow', alpha=0.9)
        axs.vlines(sun_rise, 3, 1, linestyle=':', color='yellow', alpha=0.9)
        #plt.legend(loc=2)
        plt.grid(color='gainsboro', linestyle='-', linewidth=1, alpha=0.3)
        plt.title('Visibility plot for the night of the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' on ' + str(telescope))


def gantt_chart_all(target_list=None):
    """

    Parameters
    ----------
    target_list : path
        path of the target list

    Returns
    -------
    plot
        gant chart of all scheduled targets so far

    """
    if target_list is None:
        target_list = target_list_from_stargate_path
    target_table_spc = pd.read_csv(target_list, delimiter=',')
    all_targets = target_table_spc['Sp_ID']
    files = []
    start = []
    finish = []
    targets = []
    configuration = []
    telescope = []
    location = EarthLocation.from_geodetic(-70.40300000000002 * u.deg, -24.625199999999996 * u.deg,
                                           2635.0000000009704 * u.m)
    telescopes = ['Io', 'Europa', 'Ganymede', 'Callisto', 'Artemis', 'Saint-Ex','TS_La_Silla','TN_Oukaimeden']
    for tel in telescopes:
        for i in os.listdir(os.path.join(path_spock + '/DATABASE/', tel,
                                         'Archive_night_blocks')):
            if i.endswith('.txt'):
                df = pd.read_csv(os.path.join(path_spock + '/DATABASE/', tel,
                                              'Archive_night_blocks', i), sep=' ', skipinitialspace=True)
                start.append(list(df["start time (UTC)"]))
                finish.append(list(df["end time (UTC)"]))
                targets.append(list(df['target']))
                configuration.append(list(df['configuration']))
                telescope.append(tel)
    df2 = pd.DataFrame(
        {'start': start, 'finish': finish, 'targets': targets, 'telescope': telescope, 'configuration': configuration})

    Resource=df2['telescope']
    Task=df2['targets']
    Start=df2['start']
    Finish=df2['finish']
    Description=df2['configuration']
    df=[]
    colors = {'Io': 'rgba(220, 0, 0,0.3)',
              'Europa': 'rgba(0, 0, 255,0.75)',
              'Callisto': 'rgba(0, 255, 255,0.75)',
              'Ganymede': 'rgba(255, 128, 0,0.75)',
              'Artemis':'rgba(107,142,35,0.75)',
              'Saint-Ex':'rgba(255,215,0,0.9)',
              'Io_s': 'rgba(255, 182, 193, .9)',
              'Europa_s': 'rgba(28,134,238,0.9)',
              'Ganymede_s': 'rgba(255,160,122,0.9)',
              'Callisto_s': 'rgba(152,245,255,0.9)',
              'TS_La_Silla': 'rgba(255,0,255,0.9)',
              'TN_Oukaimeden':'rgba(0,128,128,0.9)'}

    for i in range(0,len(Task)):
        for j in range(0,len(Task[i])):
            df.append(dict(Task=Task[i][j], Start=Start[i][j], Finish=Finish[i][j],Resource=Resource[i],Description=Description[i][j]))

    fig = ff.create_gantt(df, colors=colors, index_col='Resource', show_colorbar=True, showgrid_x=True, showgrid_y=True, group_tasks=True)

    fig['layout'].update(height=800, width=1300, title='Preview: Plans sent to telescopes')
    fig['layout'].update(autosize=False, margin= go.layout.Margin(l=100))
    config = {
        'scrollZoom': True
    }
    offline.plot(fig,auto_open=True,filename=path_spock + '/SPOCK_Figures/Preview_schedule.html',config=config)


def gantt_chart(date_start, date_end, telescope, local=False):
    """

    Parameters
    ----------
    date_start : date
        date start 'yyyy-mm-dd'
    date_end : date
        date start 'yyyy-mm-dd'
    telescope : str
        name of the  telescope

    Returns
    -------
    plot
        gant chart on a given range of  days

    """
    start = []
    finish = []
    targets = []
    configuration = []
    list_night_blocks = []
    telescopes = []
    if isinstance(telescope, str):
        telescope = [telescope]
    date_range_in_days = int((Time(date_end) - Time(date_start)).value)

    # Read night block from archive
    if local is False:
        for tel in telescope:
            for i in range(0, date_range_in_days):
                day = date_start + i
                df = SPOCKST.read_night_block(telescope=tel, day=day.tt.datetime.strftime("%Y-%m-%d"))
                start.append(list(df["start time (UTC)"]))
                finish.append(list(df["end time (UTC)"]))
                targets.append(list(df['target']))
                configuration.append(list(df['configuration']))
                telescopes.append(tel)

        df2 = pd.DataFrame(
        {'start': start, 'finish': finish, 'targets': targets, 'telescope': telescopes, 'configuration': configuration})

        Resource = df2['telescope']
        Task = df2['targets']
        Start = df2['start']
        Finish = df2['finish']
        Description = df2['configuration']
        df = []
        colors = {'Io': 'rgba(220, 0, 0,0.3)',
                  'Europa': 'rgba(0, 0, 255,0.75)',
                  'Callisto': 'rgba(0, 255, 255,0.75)',
                  'Ganymede': 'rgba(255, 128, 0,0.75)',
                  'Artemis':'rgba(107,142,35,0.75)',
                  'Saint-Ex':'rgba(255,255,0,0.75)',
                  'Io_s': 'rgba(255, 182, 193, .9)',
                  'Europa_s': 'rgba(28,134,238,0.9)',
                  'Ganymede_s': 'rgba(255,160,122,0.9)',
                  'Callisto_s': 'rgba(152,245,255,.9)',
                  'TS_La_Silla': 'rgba(255,0,255,0.9)',
                  'TN_Oukaimeden': 'rgba(0,128,128,0.9)'}

        for i in range(0,len(Task)):
            for j in range(0,len(Task[i])):
                df.append(dict(Task=Task[i][j], Start=Start[i][j], Finish=Finish[i][j],Resource=Resource[i],Description=Description[i][j]))

        fig = ff.create_gantt(df, colors=colors, index_col='Resource', show_colorbar=True, showgrid_x=True, showgrid_y=True, group_tasks=True)

        fig['layout'].update(height=800, width=1300, title='Preview: Plans sent to telescopes')
        fig['layout'].update(autosize=False, margin=go.layout.Margin(l=100))
        config = {
            'scrollZoom': True}
        offline.plot(fig,auto_open=True,filename=path_spock + '/SPOCK_Figures/Preview_schedule.html', config=config)

    else:
        for tel in telescope:
            for i in range(0, date_range_in_days):
                day = date_start + i
                list_night_blocks.append('night_blocks_' + tel + '_' + day.tt.datetime.strftime("%Y-%m-%d") + '.txt')
            for i in os.listdir(os.path.join(path_spock + '/DATABASE/', tel,
                                             'Archive_night_blocks')):
                if i.endswith('.txt'):
                    if any(i in s for s in list_night_blocks):
                        df = pd.read_csv(os.path.join(path_spock + '/DATABASE/', tel,
                                                      'Archive_night_blocks', i), sep=' ', skipinitialspace=True)
                        start.append(list(df["start time (UTC)"]))
                        finish.append(list(df["end time (UTC)"]))
                        targets.append(list(df['target']))
                        configuration.append(list(df['configuration']))
                        telescopes.append(tel)

        df2 = pd.DataFrame(   {'start': start, 'finish': finish, 'targets': targets, 'telescope': telescopes,
                               'configuration': configuration})
        Resource = df2['telescope']
        Task = df2['targets']
        Start = df2['start']
        Finish = df2['finish']
        Description = df2['configuration']
        df = []
        colors = {'Io': 'rgba(220, 0, 0,0.3)',
                  'Europa': 'rgba(0, 0, 255,0.75)',
                  'Callisto': 'rgba(0, 255, 255,0.75)',
                  'Ganymede': 'rgba(255, 128, 0,0.75)',
                  'Artemis':'rgba(107,142,35,0.75)',
                  'Saint-Ex':'rgba(255,255,0,0.75)',
                  'Io_s': 'rgba(255, 182, 193, .9)',
                  'Europa_s': 'rgba(28,134,238,0.9)',
                  'Ganymede_s': 'rgba(255,160,122,0.9)',
                  'Callisto_s': 'rgba(152,245,255,.9)',
                  'TS_La_Silla': 'rgba(255,0,255,0.9)',
                  'TN_Oukaimeden': 'rgba(0,128,128,0.9)'}

        for i in range(0,len(Task)):
            for j in range(0,len(Task[i])):
                df.append(dict(Task=Task[i][j], Start=Start[i][j], Finish=Finish[i][j],Resource=Resource[i],Description=Description[i][j]))

        fig = ff.create_gantt(df, colors=colors, index_col='Resource', show_colorbar=True, showgrid_x=True, showgrid_y=True, group_tasks=True)

        fig['layout'].update(height=800, width=1300, title='Preview: Plans sent to telescopes')
        fig['layout'].update(autosize=False, margin=go.layout.Margin(l=100))
        config = {
            'scrollZoom': True }
        offline.plot(fig,auto_open=True,filename=path_spock + '/SPOCK_Figures/Preview_schedule.html',config=config)


def airmass_altitude_plot_given_target(name_observatory, day, target, path_target_list=None):
    """
    Parameters
    ----------
    name_observatory : str
        name of the observatory
    day  : date
        date in format  'yyyy-mm-dd'
    target :  str
         name of  target
    path_target_list  : path
        path of the target  list

    Returns
    -------

    """
    if path_target_list is None:
        path_target_list = target_list_from_stargate_path
    observatory = charge_observatories(name_observatory)[0]
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day, which='next').jd - 0.07, \
                                      observatory.twilight_morning_nautical(day + 1, which='nearest').jd + 0.07, 100),
                          format='jd')
    sun_set = observatory.twilight_evening_nautical(day, which='next').iso
    sun_rise = observatory.twilight_morning_nautical(day + 1, which='nearest').iso
    target_list = pd.read_csv(path_target_list, delimiter=',')
    idx_target_list = list(target_list['Sp_ID']).index(target)

    fig, axs = plt.subplots(1, figsize=(8, 6))
    colors_start_new_target = ['black', 'darkgray', 'lightgray']

    dec = target_list['DEC'][idx_target_list]
    ra = target_list['RA'][idx_target_list]
    # Moon distance
    t_midnight = Time(sun_rise) - ((Time(sun_rise) - Time(sun_set)) / 2)
    moon = get_moon(time=t_midnight, location=observatory.location)
    distance_moon = round(moon.separation(SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))).value, 2)

    plot_styles = {'linestyle': '-', 'color': 'k'}
    if name_observatory == 'SSO':
        plot_styles = {'linestyle': '-', 'color': 'skyblue'}
    if name_observatory == 'SNO':
        plot_styles = {'linestyle': '-', 'color': 'teal'}
    if name_observatory == 'Saint-Ex':
        plot_styles = {'linestyle': '-', 'color': 'gold'}
    plot_airmass(FixedTarget(coord=SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg)), \
                             name=target), observatory, delta_midnight, brightness_shading=True, altitude_yaxis=True,
                 style_kwargs=plot_styles)
    # axs.vlines(sun_set, 3, 1, linestyle='--', color='orange', alpha=0.9, linewidth=2)
    # axs.vlines(sun_rise, 3, 1, linestyle='--', color='orange', alpha=0.9, linewidth=2)

    plt.ylabel('Altitude (degrees)')
    plt.grid(color='gainsboro', linestyle='-', linewidth=1, alpha=0.3)
    if distance_moon < 30.:
        plt.title('Visibility  plot for ' + target + ' on the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' at '
                  + name_observatory + '\n' + r'Moon is at: ' + str(distance_moon) + " degress"
                  , y=-0.01, color='purple')
    else:
        plt.title('Visibility  plot for ' + target + ' on the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' at '
                  + name_observatory + '\n' + r'Moon is at: ' + str(distance_moon) + " degress"
                  , y=-0.01, color='black')
    plt.show()


def airmass_altitude_plot_nolist(name_observatory, day, target, ra, dec):
    """
    Parameters
    ----------
    name_observatory : str
        name of the observatory
    day  : date
        date in format  'yyyy-mm-dd'
    target :  str
         name of  target
    path_target_list  : path
        path of the target  list

    Returns
    -------

    """

    observatory = charge_observatories(name_observatory)[0]
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day, which='next').jd - 0.07,
                                      observatory.twilight_morning_nautical(day + 1, which='nearest').jd + 0.07, 100),
                          format='jd')
    sun_set = observatory.twilight_evening_nautical(day, which='next').iso
    sun_rise = observatory.twilight_morning_nautical(day + 1, which='nearest').iso

    fig, axs = plt.subplots(1, figsize=(8, 6))
    colors_start_new_target = ['black', 'darkgray', 'lightgray']

    # Moon distance
    t_midnight = Time(sun_rise) - ((Time(sun_rise) - Time(sun_set)) / 2)
    moon = get_moon(time=t_midnight, location=observatory.location)
    distance_moon = round(float(moon.separation(SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))).value), 2)

    plot_styles = {'linestyle': '-', 'color': 'k'}
    if name_observatory == 'SSO':
        plot_styles = {'linestyle': '-', 'color': 'skyblue'}
    if name_observatory == 'SNO':
        plot_styles = {'linestyle': '-', 'color': 'teal'}
    if name_observatory == 'Saint-Ex':
        plot_styles = {'linestyle': '-', 'color': 'gold'}
    plot_airmass(FixedTarget(coord=SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg)), name=target),
                 observatory, delta_midnight, brightness_shading=True, altitude_yaxis=True,
                 style_kwargs=plot_styles)

    plt.ylabel('Altitude (degrees)')
    plt.grid(color='gainsboro', linestyle='-', linewidth=1, alpha=0.3)
    if distance_moon < 30.:
        plt.title('Visibility  plot for ' + target + ' on the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' at '
                  + name_observatory + '\n' + r'Moon is at: ' + str(distance_moon) + " degress"
                  , y=-0.01, color='red')
    else:
        plt.title('Visibility  plot for ' + target + ' on the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' at '
                  + name_observatory + '\n' + r'Moon is at: ' + str(distance_moon) + " degress"
                  , y=-0.01, color='black')
    plt.show()


def constraints_scores(constraints,target,observatory,start,end):
    """

    Parameters
    ----------
    constraints : astroplan.constraints
         constraints
    target : astroplan.FixedTarget
        target
    observatory : observatory
        astroplan.observed
    start : date
        start  date
    end : date
        end date

    Returns
    -------
    plot
        plot of the constraints scores

    """
    time_resolution = 0.5 * u.hour
    time_grid = time_grid_from_range([start, end],time_resolution=time_resolution)

    observability_grid = np.zeros((len(constraints), len(time_grid)))

    for i, constraint in enumerate(constraints):
        # Evaluate each constraint
        observability_grid[i, :] = constraint(observatory, target, times=time_grid)

    print(observability_grid)

    # Create plot showing observability of the target:

    extent = [-0.5, -0.5 + len(time_grid), -0.5, 2.5]

    fig, ax = plt.subplots()
    ax.imshow(observability_grid, extent=extent)

    ax.set_yticks(range(0, 3))
    ax.set_yticklabels([c.__class__.__name__ for c in constraints])

    ax.set_xticks(range(len(time_grid)))
    ax.set_xticklabels([t.datetime.strftime("%H:%M") for t in time_grid])

    ax.set_xticks(np.arange(extent[0], extent[1]), minor=True)
    ax.set_yticks(np.arange(extent[2], extent[3]), minor=True)

    ax.grid(which='minor', color='w', linestyle='-', linewidth=2)
    ax.tick_params(axis='x', which='minor', bottom='off')
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
    plt.setp(ax.get_yticklabels(), rotation=45, ha='right')

    ax.tick_params(axis='y', which='minor', left='off')
    ax.set_xlabel('Time on {0} UTC'.format(time_grid[0].datetime.date()))
    fig.subplots_adjust(left=0.35, right=0.9, top=0.9, bottom=0.1)


def coverage(t, p):
    if p == 0:
        return 1
    else:
        ph = ((t + 0.5*p)%p - (0.5*p))
        sph_in = np.sort(ph)
        sph_out = sph_in[::-1]

        sph_in_diff = np.abs(np.diff(sph_in))
        sph_out_diff = np.abs(np.diff(sph_out))

        df = np.min(np.diff(t))

        spaces_in = np.sort(sph_in[np.hstack([*np.argwhere(sph_in_diff > 4*df).T, len(sph_in)-1])])
        spaces_out = np.sort(sph_out[np.hstack([*np.argwhere(sph_out_diff > 4*df).T, len(sph_in)-1])])

        return np.sum(spaces_in - spaces_out)/p


def getSPClcV2(target, ap = '', pwvCorr = 0, user= user_portal, password= pwd_portal):
    urlGet = "http://www.mrao.cam.ac.uk/SPECULOOS/portal/get_tls_prep_v2.php?date=*&id=" + target + \
            "&filter=&telescope=&ap=" + str(ap) + "&pwvCorr=" + str(pwvCorr)
    rGET = requests.get(urlGet, auth=(user, password))
    names = ['TMID-2450000', 'BJDMID-2450000', 'DIFF_FLUX', 'ERROR', 'DIFF_FLUX_PWV',
             'RA_MOVE', 'DEC_MOVE', 'FWHM', 'PSF_a_5', 'PSF_b_5', 'SKYLEVEL', 'AIRMASS',
             'EXPOSURE']
    lc = pd.read_csv(io.StringIO(rGET.text), header=None, names=names, delimiter='\s+', index_col=None)
    lc.reset_index(drop=True)
    return lc


def getSPCdata(target, date, telescope='any', ap=6, user= user_portal, password=pwd_portal):
    if telescope == 'any':
        urlGet = "http://www.mrao.cam.ac.uk/SPECULOOS/portal/get_file.php?telescope=&date=" + date + "&id=" + target + "&filter=&file=MCMC_text_" + str(ap)
    else:
        urlGet = "http://www.mrao.cam.ac.uk/SPECULOOS/portal/get_file.php?telescope=" + telescope + '&date=' + date + "&id=" + target + "&filter=&file=MCMC_text_" + str(ap)
    rGET = requests.get(urlGet, auth=(user, password))
    rFile = json.loads(rGET.text)
    # retrieve available file
    try:
        urlMCMC = 'http://www.mrao.cam.ac.uk/SPECULOOS/portal/' + rFile[0]
        rMCMC = requests.get(urlMCMC, auth=(user, password))
        # read files
        MCMC = pd.read_csv(io.StringIO(rMCMC.text), sep=' ')
        MCMC['JD'] = MCMC['TMID-2450000'] + 2450000
        targetdf = pd.DataFrame({'TMID-2450000':MCMC['JD'].values,'DIFF_FLUX': MCMC['DIFF_FLUX'].values,
                                 'ERROR': MCMC['ERROR'].values,'dx_MOVE':MCMC['RA_MOVE'].values,
                                 'dy_MOVE':MCMC['DEC_MOVE'].values,'FWHM':MCMC['FWHM'].values,
                                 'AIRMASS': MCMC['AIRMASS'].values,'EXPOSURE':MCMC['EXPOSURE'].values})
        return targetdf
    except TypeError:
        targetdf = pd.DataFrame({'JD':[],'DIFF_FLUX': [], 'ERROR': [], 'AIRMASS': []})
        return targetdf.reset_index(drop=True)


def get_all_LCS(gaia_id_target, fix_expt=None):
    times = []
    diff_fluxes = []
    exposures = []
    dates = []

    user = 'educrot'
    password = '9UCExnjwes'
    url = "http://www.mrao.cam.ac.uk/SPECULOOS/speculoos-portal/php/get_observations.php"
    res = requests.get(url, auth=HTTPBasicAuth(user_portal, pwd_portal))
    cache = json.loads(res.content)

    for i in range(len(cache[0])):
        if cache[0][i]["gaia_id"] == str(gaia_id_target):
            print(Fore.GREEN + 'INFO:  ' + Fore.BLACK + "Downloading LC of " + cache[0][i]["sp_id"] +
                  " on " + cache[0][i]["telescope"] + " the " + str(cache[0][i]["date"]))
            url = "http://www.mrao.cam.ac.uk/SPECULOOS/portal-night-data/" + str(cache[0][i]["telescope"]) + \
                  '_' + str(cache[0][i]["date"]) + '_' + str(cache[0][i]["filter"]) + '_' + str(gaia_id_target) + '_' + \
                  str(cache[0][i]["sp_id"].split(" ")[0]) + '.json'
            res = requests.get(url, auth=HTTPBasicAuth(user_portal, pwd_portal))
            try:
                data = json.loads(res.content)
                # counts = collections.Counter(exposures)
                if fix_expt is None:
                    times.append(np.array(data['environment']['BJD-OBS']))
                    diff_fluxes.append(data["stars"][str(data['best_ap'])][0]['DIFF_FLUX'])
                    exposures.append(data['environment']['EXPOSURE'])
                    dates.append(data['date'])
                else:
                    if data['environment']['EXPOSURE'] == fix_expt:
                        times.append(np.array(data['environment']['BJD-OBS']))
                        diff_fluxes.append(data["stars"][str(data['best_ap'])][0]['DIFF_FLUX'])
                        exposures.append(data['environment']['EXPOSURE'])
                        dates.append(data['date'])
                        print(exposures)
            except json.decoder.JSONDecodeError:
                print(Fore.RED + 'ERROR:  ' + Fore.BLACK + "Can not download LC of " + cache[0][i]["sp_id"] +
                      " on " + cache[0][i]["telescope"] + " the " + str(cache[0][i]["date"]))
                pass

    time = np.sort(np.concatenate(times))
    diff_flux = np.concatenate(diff_fluxes)
    return time, diff_flux, exposures, dates, str(cache[0][i]["sp_id"].split(" ")[0]), times


def phase_coverage_given_target(target, pmin, pmax, fix_expt=None, path_target_list=None, times=None):
    if path_target_list is None:
        path_target_list = target_list_from_stargate_path
    schedules_st = SPOCKST.Schedules()
    schedules_st.load_parameters()
    target_list_follow_up = schedules_st.target_table_spc_follow_up
    target_list_follow_up = target_list_follow_up.sort_values(by="Sp_ID")
    target_list_speculoos = pd.read_csv(path_target_list, sep=',')
    target_list_speculoos = target_list_speculoos.sort_values(by="Sp_ID")

    all_targets = pd.DataFrame({'Sp_ID': pd.concat([target_list_follow_up['Sp_ID'],
                                                   target_list_speculoos["Sp_ID"]]),
                                "Gaia_ID": pd.concat([target_list_follow_up['Gaia_ID'],
                                                     target_list_speculoos["Gaia_ID"]])})
    all_targets.reset_index(inplace=True)

    idx_target_list = list(all_targets['Sp_ID']).index(target)
    if times is None:
        # data = getSPClcV2(target=target, ap='', pwvCorr=0)
        #
        # t = data['BJDMID-2450000']
        # t = t.fillna(0)
        # t = np.sort(t)
        t, diff_flux, exposures, dates, target_name, times = get_all_LCS(
            gaia_id_target=all_targets['Gaia_ID'][idx_target_list],
            fix_expt=fix_expt)
    else:
        t = times

    P_min = pmin
    P_max = pmax
    periods = np.arange(P_min, P_max, 0.01)

    try:
        covs = [coverage(t, period) for period in periods]
        mean_cov = np.mean(covs)*100
        fig, ax = plt.subplots(1, figsize=(9, 7))
        if target[0:2] == "Sp":
            anchored_text = AnchoredText(r'$SNR_{JWST}$ = ' +
                                         str(round(target_list_speculoos['SNR_JWST_HZ_tr'][idx_target_list], 3))
                                         + "\n" +
                                         "Hours observed = " +
                                         str(round(target_list_speculoos['nb_hours_surved'][idx_target_list], 2)) +
                                         ' hours',  loc=3)
            ax.add_artist(anchored_text)

        plt.plot(periods, np.array(covs)*100, c="silver", label='Effective cov = ' + str(round(mean_cov, 1)) + ' %')
        plt.plot(periods, np.array(covs)*100, ".", c="k",)
        plt.ylabel('Phase coverage in %')
        plt.xlabel('Period in days')
        plt.legend(fontsize=16)
        plt.grid(color='gainsboro', linestyle='-', linewidth=1, alpha=0.3)
        plt.title(r'Phase coverage for ' + target + r'with periods $\in$' + str(pmin) + ' - ' + str(pmax))
        plt.show()
        print(Fore.YELLOW + 'WARNING:  ' + Fore.BLACK + ' If  you feel the coverage is not consistent with ' +
              'the number of hours observed check if the exposure' +
              'time has been changed along the observations. ')

        # exposure time vs dates
        if len(dates) > 10:
            plt.figure(figsize=(int(len(dates) / 4), 4))
        else:
            plt.figure(figsize=(8, 7))
        plt.plot(dates, exposures, 'H', color='goldenrod', alpha=0.8)
        plt.ylim(min(exposures)-1, max(exposures)+1)
        plt.ylabel('Exposure time (seconds)')
        plt.xticks(rotation=60)
        plt.xlabel('Dates', )
        plt.show()

    except ValueError:
        print(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' No data for this target ! ')


def phase_coverage(period, times, midpoint=0, binning=0.005):
    _times = np.array(
        [
            ((time - (np.round((time - midpoint) / period) * period)) - midpoint)
            / period
            for time in times
        ]
    )
    _stacked_times = np.hstack(_times)
    bin_points = len(
        binning_for_cov(np.ones(len(_stacked_times)), _stacked_times, binning / period)[0]
    )
    total_bin_points = len(np.arange(-0.5, 0.5, binning / period))
    return bin_points / total_bin_points


def plot_polar_phase_coverage(midpoint, period, _times, r, d):
    R = 1 * d  # 5*(period/(2*np.pi))

    times = np.array([np.array([np.min(t), np.max(t)]) for t in _times])

    length = [(np.max(t) - np.min(t)) / period for t in times]
    offseted_times = np.array(
        [
            ((time - (np.round((time - midpoint) / period) * period)) - midpoint)
            / period
            for time in times
        ]
    )
    fig, ax1 = plt.subplots(figsize=(4, 4))
    plt.pie([1], colors=['#016CA000'], radius=r)

    for l, t in zip(length, offseted_times):
        plot_slice(t[0], l, R=R)

    my_circle = plt.Circle((0, 0), R - 0.3, color='white')
    planet = plt.Circle((-np.cos(2 * np.pi * (-midpoint / period) + np.pi / 2) * (R - .15),
                         (R - .15) * np.sin(2 * np.pi * (midpoint / period) + np.pi / 2)), 0.1,
                        facecolor='black', edgecolor='black')
    star = plt.Circle((0, 0), 0.5, facecolor='brown', edgecolor='brown', alpha=0.5)
    p = plt.gcf()
    p.gca().add_artist(my_circle)
    p.gca().add_artist(star)
    p.gca().add_artist(planet)
    # plt.legend([my_circle],['For Period ' + str(round(period,2)) + ' days'],loc=4, fontsize = 'x-large')
    plt.autoscale()
    return phase_coverage(np.abs(period), _times)


def plot_slice(offset, width, color='#016CA050', R=3):
    size = 1
    plt.pie([width * 100, 100 - width * 100], colors=[color, [0, 0, 0, 0]], counterclock=False,
            startangle=90 - offset * 360, radius=R, \
            wedgeprops=dict(width=size))


def binning_for_cov(
        flux,
        jd,
        bin_size,
        error=None,
        std=False,
        mean_method=np.mean,
        mean_error_method=lambda l: np.sqrt(np.sum(np.power(l, 2))) / len(l), ):
    # TODO: take the middle of the binned time
    bins = np.arange(np.min(jd), np.max(jd), bin_size)
    d = np.digitize(jd, bins)

    final_bins = []
    binned_flux = []
    if error is not None:
        binned_error = []
    _std = []

    for i in range(1, np.max(d) + 1):
        s = np.where(d == i)
        if len(s[0]) > 0:
            try:
                flux[s[0]]
            except:
                t = 8
            binned_flux.append(mean_method(flux[s[0]]))
            final_bins.append(np.mean(jd[s[0]]))
            _std.append(np.std(flux[s[0]]) / np.sqrt(len(s[0])))
            if error is not None:
                binned_error.append(mean_error_method(error[s[0]]))
        # else:
        #     np.delete(bins, np.where(bins == bins[i])[0])
        #     print("deleted {}".format(np.where(bins == bins[i])))

    if std:
        return final_bins, binned_flux, np.array(_std)
    elif error is not None and type(error) in [np.array, np.ndarray, list]:
        return final_bins, binned_flux, binned_error
    else:
        return final_bins, binned_flux


def plot_annulus_phase_covered(times, period, target_name, t0):
    planet_pos = (t0 - min(np.concatenate(times))) % period
    cov = plot_polar_phase_coverage(midpoint=planet_pos, period=period, _times=times, r=1.5, d=3)
    plt.title('For period ' + str(period) + ' days ' + '\n coverage is ' + str(cov) + "%", fontsize=18, y=-0.2)
    # plt.savefig('phase_coverage_' + str(period) + '_' + str(target_name) + '.png', bbox_inches='tight')


def phase_folded_LC(t, diff_flux, period, t0, x_lim_phase=0.02):
    x_fold = (t - t0 + 0.5 * period) % period - 0.5 * period
    plt.subplots(1, figsize=(8, 5))
    plt.errorbar(x_fold, diff_flux, fmt='.', alpha=0.5,)
    plt.vlines(t0,0.98,1.02,'r',alpha=0.3,zorder=3)
    plt.title("Phase folded on period " + str(period) + " days")
    plt.ylim(0.98, 1.02)
    plt.xlim(-x_lim_phase, x_lim_phase)
    plt.show()
