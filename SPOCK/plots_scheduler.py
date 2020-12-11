import chart_studio
from plotly import graph_objs as go
import pandas as pd
from plotly import offline
from astroplan import Observer
from astropy.time import Time
from astropy import units as u
from astroplan.plots import dark_style_sheet, plot_airmass
from astroplan import FixedTarget
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, get_sun, AltAz, EarthLocation
import numpy as np
import plotly.figure_factory as ff
from astroplan.utils import time_grid_from_range
from colorama import Fore
import os
import pkg_resources
import requests
import SPOCK.long_term_scheduler as SPOCKLT
import yaml
import io
from matplotlib.offsetbox import AnchoredText
from SPOCK  import pwd_appcs,pwd_HUB,user_portal,pwd_portal,pwd_appcs,pwd_SNO_Reduc1,user_chart_studio,pwd_chart_studio,path_spock



# pwd_appcs,pwd_HUB,user_portal,pwd_portal,pwd_appcs,pwd_SNO_Reduc1,user_chart_studio,pwd_chart_studio,path_spock = SPOCKLT._get_files()
chart_studio.tools.set_credentials_file(username=user_chart_studio, api_key=pwd_chart_studio)
chart_studio.tools.set_config_file(world_readable=True,sharing='public')

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
        target_list = path_spock + '/target_lists/speculoos_target_list_v6.txt'
    target_table_spc = pd.read_csv(target_list,delimiter=' ')
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

def gantt_chart(date_start,date_end,telescope):
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

    df2 = pd.DataFrame(
        {'start': start, 'finish': finish, 'targets': targets, 'telescope': telescopes, 'configuration': configuration})

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
        'scrollZoom': True
    }
    offline.plot(fig,auto_open=True,filename=path_spock + '/SPOCK_Figures/Preview_schedule.html',config=config)

def airmass_altitude_plot_given_target(name_observatory,day,target,path_target_list=None):
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
        path_target_list = path_spock + '/target_lists/speculoos_target_list_v6.txt'
    observatory = charge_observatories(name_observatory)[0]
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day, which='next').jd - 0.07, \
                                      observatory.twilight_morning_nautical(day + 1, which='nearest').jd + 0.07, 100),
                          format='jd')
    sun_set = observatory.twilight_evening_nautical(day, which='next').iso
    sun_rise = observatory.twilight_morning_nautical(day + 1, which='nearest').iso
    target_list = pd.read_csv(path_target_list, delimiter=' ')
    idx_target_list = list(target_list['Sp_ID']).index(target)

    fig, axs = plt.subplots(1,figsize=(9,7))
    colors_start_new_target = ['black', 'darkgray', 'lightgray']

    dec = target_list['DEC'][idx_target_list]
    ra = target_list['RA'][idx_target_list]
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
    #axs.vlines(sun_set, 3, 1, linestyle='--', color='orange', alpha=0.9, linewidth=2)
    #axs.vlines(sun_rise, 3, 1, linestyle='--', color='orange', alpha=0.9, linewidth=2)

    plt.ylabel('Altitude (degrees)')
    plt.grid(color='gainsboro', linestyle='-', linewidth=1, alpha=0.3)
    plt.title('Visibility  plot for ' + target + ' on the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' at '
              + name_observatory,y=-0.01)
    #plt.title('Visibility plot for target ' + target + ' on the ' + str(day.tt.datetime.strftime("%Y-%m-%d")))
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
    if p==0:
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
    urlGet ="http://www.mrao.cam.ac.uk/SPECULOOS/portal/get_tls_prep_v2.php?date=*&id=" + target + "&filter=&telescope=&ap=" + str(ap) + "&pwvCorr=" + str(pwvCorr)
    rGET = requests.get(urlGet, auth=(user, password))
    names = ['TMID-2450000', 'BJDMID-2450000', 'DIFF_FLUX', 'ERROR', 'DIFF_FLUX_PWV',
    'RA_MOVE', 'DEC_MOVE', 'FWHM', 'PSF_a_5', 'PSF_b_5', 'SKYLEVEL', 'AIRMASS',
    'EXPOSURE']
    lc = pd.read_csv(io.StringIO(rGET.text), header=None, names=names, delimiter='\s+',index_col=None)
    lc.reset_index(drop=True)
    return lc

def getSPCdata(target, date, telescope='any', ap=6, user= user_portal, password= pwd_portal):
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
        targetdf = pd.DataFrame({'TMID-2450000':MCMC['JD'].values,'DIFF_FLUX': MCMC['DIFF_FLUX'].values, 'ERROR': MCMC['ERROR'].values,
                 'dx_MOVE':MCMC['RA_MOVE'].values,'dy_MOVE':MCMC['DEC_MOVE'].values,'FWHM':MCMC['FWHM'].values,
                   'FWHMx':MCMC['PSF_a_5'].values,'FWHMy':MCMC['PSF_b_5'].values,'SKYLEVEL':MCMC['SKYLEVEL'].values,
                               'AIRMASS': MCMC['AIRMASS'].values,'EXPOSURE':MCMC['EXPOSURE'].values})
        return targetdf
    except TypeError:
        targetdf = pd.DataFrame({'JD':[],'DIFF_FLUX': [], 'ERROR': [], 'AIRMASS': []})
        return targetdf.reset_index(drop=True)

def phase_coverage_given_target(target,pmin,pmax,path_target_list=None):

    if path_target_list is None:
        path_target_list = path_spock + '/target_lists/speculoos_target_list_v6.txt'

    target_list = pd.read_csv(path_target_list,sep=' ')
    idx_target_list = list(target_list['Sp_ID']).index(target)
    data = getSPClcV2(target=target, ap='', pwvCorr=0)

    t = data['BJDMID-2450000']
    t = t.fillna(0)

    colors_start_new_target = ['black', 'darkgray', 'lightgray']


    P_min = pmin
    P_max = pmax
    periods = np.arange(P_min, P_max, 0.01)

    try:
        covs = [coverage(t, period) for period in periods]
        mean_cov = np.mean(covs)*100
        fig, ax = plt.subplots(1, figsize=(9, 7))
        anchored_text = AnchoredText(r'$SNR_{JWST}$ = ' +  str(round(target_list['SNR_JWST_HZ_tr'][idx_target_list],3))
                                     + "\n" +
                                     "Hours observed = " + str(round(target_list['nb_hours_surved'][idx_target_list],2)) + ' hours',
                                     loc=3)

        plt.plot(periods, covs, c="silver",label='Effective cov = ' +  str(round(mean_cov,1)) + ' %')
        plt.plot(periods, covs, ".", c="k",)
        # ax.annotate(r'$SNR_{JWST}$ = ' +  str(round(target_list['SNR_JWST_HZ_tr'][idx_target_list],3)), (3,1),
        #             xytext=(0.8, 0.85), textcoords='axes fraction',
        #             fontsize=16,
        #             horizontalalignment='right', verticalalignment='top')
        # ax.annotate(r'Effective cov = ' +  str(round(mean_cov,1)) + ' %', (3,1),
        #             xytext=(0.8, 0.8), textcoords='axes fraction',
        #             fontsize=16,
        #             horizontalalignment='right', verticalalignment='top')
        ax.add_artist(anchored_text)
        plt.ylabel('Phase coverage in %')
        plt.xlabel('Period in days')
        plt.legend(fontsize=16)
        plt.grid(color='gainsboro', linestyle='-', linewidth=1, alpha=0.3)
        plt.title(r'Phase coverage for ' + target + ' with periods $\in$ ' + str(pmin) + ' - ' + str(pmax) )
        plt.show()
    except ValueError:
        print(Fore.RED + 'ERROR:  ' + Fore.BLACK + ' No data for this target ! ')
