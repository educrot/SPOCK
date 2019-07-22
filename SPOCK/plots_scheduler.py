import plotly
plotly.tools.set_credentials_file(username='ed510', api_key='lN1JDlEfs0FPrLHqPScL')
plotly.tools.set_config_file(world_readable=True,sharing='public')
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import pandas as pd
from plotly import offline
from astroplan import Observer
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, get_sun, AltAz, EarthLocation
import numpy as np
import plotly.plotly as py
import plotly.figure_factory as ff

init_notebook_mode(connected=True)
pd.set_option('display.max_rows', 500)

def gantt_chart(plan_file,observatory):
    df=[]

    df2 = pd.read_csv(plan_file+'.txt', delimiter=' ',skipinitialspace=True)
    Task=df2['Name']
    Start=df2['Start']
    Finish=df2['Finish']
    #Complete=df2['Priority']
    Resource=df2['Telescope']

    for i in range(0,len(Task)):
        df.append(dict(Task=Task[i], Start=Start[i], Finish=Finish[i],Resource=Resource[i],Description=[df2['nb_hours_start'][i],df2['nb_hours_end'][i]])) #Complete=Complete[i])) #Description=Description


    colors = {'Io': 'rgb(220, 0, 0)',
              'Europa': 'rgb(0, 0, 255)',
              'Callisto': 'rgb(0, 255, 255)',
              'Ganymede': 'rgb(255, 128, 0)',
              'Io_s': 'rgba(255, 182, 193, .9)',
              'Europa_s': 'rgba(28,134,238,0.9)',
              'Ganymede_s': 'rgba(255,160,122,0.9)',
              'Callisto_s': 'rgba(152,245,255,.9)'}

    fig = ff.create_gantt(df, colors=colors, index_col='Resource', show_colorbar=True, showgrid_x=True, showgrid_y=True, group_tasks=True)

    location = EarthLocation.from_geodetic(-70.40300000000002*u.deg, -24.625199999999996*u.deg,2635.0000000009704*u.m)
    observatory = Observer(location=location, name="observatory", timezone="UTC")
    for i in range(0,len(df2)):
        dt_nautical_civil_evening_2=(observatory.twilight_evening_nautical(Time(Start[i],format='iso',scale='utc'),which='nearest')-observatory.twilight_evening_civil(Time(Start[i],format='iso',scale='utc'),which='nearest'))/2
        dt_civil_nautical_morning_2=(observatory.twilight_morning_civil(Time(Start[i],format='iso',scale='utc'),which='next')-observatory.twilight_morning_nautical(Time(Start[i],format='iso',scale='utc'),which='next'))/2
        twilight_evening_between_civil_nautic=Time(observatory.twilight_evening_nautical(Time(Start[i],format='iso',scale='utc'),which='nearest')-dt_nautical_civil_evening_2)
        twilight_morning_between_civil_nautic=Time(observatory.twilight_morning_civil(Time(Start[i],format='iso',scale='utc'),which='next') - dt_civil_nautical_morning_2)
        fig['layout']['shapes'].append(
            {
                'type':'rect',
                'x0': twilight_evening_between_civil_nautic.iso,
                'y0': -1,
                'x1': twilight_morning_between_civil_nautic.iso,
                'y1': 16,
                'line': {
                    'color':'rgba(128,0,128,1)',
                    'width':2,
                },
                'fillcolor': 'rgba(128,0,128,0.05)',
            })

    return fig
