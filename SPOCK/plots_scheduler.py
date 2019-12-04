import chart_studio
chart_studio.tools.set_credentials_file(username='ed510', api_key='lN1JDlEfs0FPrLHqPScL')
chart_studio.tools.set_config_file(world_readable=True,sharing='public')
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
import  os

def charge_observatories(Name):
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
        location_TNOuka = EarthLocation.from_geodetic(31.20516*u.deg, -7.862263*u.deg, 2751*u.m)
        observatories.append(Observer(location=location_TNOuka, name="TNOuka", timezone="UTC"))

    if 'Munich' in str(Name):
        location_munich= EarthLocation.from_geodetic(48.2*u.deg, -11.6*u.deg, 600*u.m)
        observatories.append(Observer(location=location_munich, name="Munich", timezone="UTC"))

    return observatories

def airmass_plot_saved(name_observatory,telescope,day):
    night_block = pd.read_csv(os.path.join('./DATABASE/', telescope,
                                              'Archive_night_blocks','night_blocks_' + telescope + '_' + day.tt.datetime.strftime("%Y-%m-%d") + '.txt'),\
                              sep=' ', skipinitialspace=True)
    observatory = charge_observatories(name_observatory)[0]
    day = Time(night_block['start time (UTC)'][0]) - 5 *u.hour
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day,which='next').jd - 0.05,\
                                      observatory.twilight_morning_nautical(day+1,which='nearest').jd + 0.05, 100),format='jd')
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
    night_block = pd.read_csv(os.path.join('./night_blocks_propositions/','night_blocks_' + telescope + '_' + day.tt.datetime.strftime("%Y-%m-%d") + '.txt'),\
                              sep=' ', skipinitialspace=True)
    observatory = charge_observatories(name_observatory)[0]
    day = Time(night_block['start time (UTC)'][0]) - 5 *u.hour
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day,which='next').jd - 0.05,\
                                      observatory.twilight_morning_nautical(day+1,which='nearest').jd + 0.05, 100),format='jd')
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

def airmass_altitude_plot_saved(name_observatory,telescope,day):
    night_block = pd.read_csv(os.path.join('./DATABASE/', telescope,
                                              'Archive_night_blocks','night_blocks_' + telescope + '_' + day.tt.datetime.strftime("%Y-%m-%d") + '.txt'),\
                              sep=' ', skipinitialspace=True)
    observatory = charge_observatories(name_observatory)[0]
    day = Time(night_block['start time (UTC)'][0]) - 5 *u.hour
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day,which='next').jd - 0.05,\
                                      observatory.twilight_morning_nautical(day+1,which='nearest').jd + 0.05, 100),format='jd')
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
        plt.legend(loc=2)
        plt.grid(color='gainsboro', linestyle='-', linewidth=1, alpha=0.3)
        plt.title('Visibility plot for the night of the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' on ' + str(telescope))

def airmass_altitude_plot_proposition(name_observatory,telescope,day):
    night_block = pd.read_csv(os.path.join('./night_blocks_propositions/','night_blocks_' + telescope + '_' + day.tt.datetime.strftime("%Y-%m-%d") + '.txt'),\
                              sep=' ', skipinitialspace=True)
    observatory = charge_observatories(name_observatory)[0]
    day = Time(night_block['start time (UTC)'][0]) - 5 *u.hour
    delta_midnight = Time(np.linspace(observatory.twilight_evening_nautical(day,which='next').jd - 0.05,\
                                      observatory.twilight_morning_nautical(day+1,which='nearest').jd + 0.05, 100),format='jd')
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
        plt.legend(loc=2)
        plt.grid(color='gainsboro', linestyle='-', linewidth=1, alpha=0.3)
        plt.title('Visibility plot for the night of the ' + str(day.tt.datetime.strftime("%Y-%m-%d")) + ' on ' + str(telescope))


def gantt_chart_all(target_list):
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
    telescopes = ['Io', 'Europa', 'Ganymede', 'Callisto', 'Artemis', 'Saint-Ex']
    for tel in telescopes:
        for i in os.listdir(os.path.join('./DATABASE/', tel,
                                         'Archive_night_blocks')):
            if i.endswith('.txt'):
                df = pd.read_csv(os.path.join('./DATABASE/', tel,
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
              'Saint-Ex':'rgba(255,255,0,0.75)',
              'Io_s': 'rgba(255, 182, 193, .9)',
              'Europa_s': 'rgba(28,134,238,0.9)',
              'Ganymede_s': 'rgba(255,160,122,0.9)',
              'Callisto_s': 'rgba(152,245,255,.9)'}

    for i in range(0,len(Task)):
        for j in range(0,len(Task[i])):
            df.append(dict(Task=Task[i][j], Start=Start[i][j], Finish=Finish[i][j],Resource=Resource[i],Description=Description[i][j]))

    fig = ff.create_gantt(df, colors=colors, index_col='Resource', show_colorbar=True, showgrid_x=True, showgrid_y=True, group_tasks=True)

    fig['layout'].update(height=800, width=1300, title='Preview: Plans sent to telescopes')
    fig['layout'].update(autosize=False, margin= go.layout.Margin(l=100))
    config = {
        'scrollZoom': True
    }
    offline.plot(fig,auto_open=True,filename='./SPOCK_Figures/Preview_schedule.html',config=config)

def gantt_chart(date_start,date_end,telescope):
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
        for i in os.listdir(os.path.join('./DATABASE/', tel,
                                         'Archive_night_blocks')):
            if i.endswith('.txt'):
                if any(i in s for s in list_night_blocks):
                    df = pd.read_csv(os.path.join('./DATABASE/', tel,
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
              'Callisto_s': 'rgba(152,245,255,.9)'}

    for i in range(0,len(Task)):
        for j in range(0,len(Task[i])):
            df.append(dict(Task=Task[i][j], Start=Start[i][j], Finish=Finish[i][j],Resource=Resource[i],Description=Description[i][j]))

    fig = ff.create_gantt(df, colors=colors, index_col='Resource', show_colorbar=True, showgrid_x=True, showgrid_y=True, group_tasks=True)

    fig['layout'].update(height=800, width=1300, title='Preview: Plans sent to telescopes')
    fig['layout'].update(autosize=False, margin=go.layout.Margin(l=100))
    config = {
        'scrollZoom': True
    }
    offline.plot(fig,auto_open=True,filename='./SPOCK_Figures/Preview_schedule.html',config=config)
