from astropy.time import Time
import io
import pandas as pd
import plotly.figure_factory as ff
from plotly import graph_objs as go
from plotly import offline
import requests
from SPOCK import user_portal, pwd_portal
import sys


def gantt_chart(date_start, date_end, telescope_to_plot):
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
    date_start = Time(date_start)
    date_end = Time(date_end)
    start = []
    finish = []
    targets = []
    configuration = []
    telescopes = []

    if isinstance(telescope_to_plot, str):
        telescope = [telescope_to_plot]
    if telescope_to_plot == 'all':
        telescope_to_plot = ['Io', 'Europa', 'Ganymede', 'Callisto', 'Artemis', 'Saint-Ex']
    date_range_in_days = int((Time(date_end) - Time(date_start)).value)

    for tel in telescope_to_plot:
        list_night_blocks = []
        for i in range(0, date_range_in_days):
            day = date_start + i
            list_night_blocks.append('night_blocks_' + tel + '_' + day.tt.datetime.strftime("%Y-%m-%d") + '.txt')
        for l in list_night_blocks:
            TargetURL = "http://www.mrao.cam.ac.uk/SPECULOOS/"+tel + \
                "/schedule/Archive_night_blocks/" + l
            resp = requests.get(TargetURL, auth=(user_portal, pwd_portal))
            content = resp.text.replace("\n", "")
            df = pd.read_csv(io.StringIO(resp.text), sep=' ')

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
            df.append(dict(Task=Task[i][j], Start=Start[i][j], Finish=Finish[i][j], Resource=Resource[i],
                           Description=Description[i][j]))

    fig = ff.create_gantt(df, colors=colors, index_col='Resource', show_colorbar=True, showgrid_x=True,
                          showgrid_y=True, group_tasks=True)

    fig['layout'].update(height=800, width=1300, title='Preview: Plans sent to telescopes')
    fig['layout'].update(autosize=False, margin=go.layout.Margin(l=100))
    config = {
        'scrollZoom': True
    }
    offline.plot(fig, auto_open=True, filename='./SPOCK_Figures/Preview_schedule.html', config=config)


if __name__ == "__main__":
    gantt_chart(sys.argv[1], sys.argv[2], sys.argv[3])
