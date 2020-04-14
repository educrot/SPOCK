import SPOCK.long_term_scheduler as SPOCKLT
import SPOCK.short_term_scheduler as SPOCKST
import SPOCK.plots_scheduler as SPOCKplot
from astropy.time import Time
import SPOCK.ETC as ETC
from astroplan import Observer, FixedTarget, is_observable, AtNightConstraint,AirmassConstraint,AltitudeConstraint,MoonSeparationConstraint
from eScheduler.constraints_spc import CelestialPoleSeparationConstraint
from astropy.time import Time
from astropy import units as u
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
from astroplan.utils import time_grid_from_range


# ---------------------- LONG TERM SCHEDULER ---------------------

obs = 5
schedule = SPOCKLT.Schedules()
schedule.load_parameters('./input.csv',obs)
#schedule.make_schedule(Altitude_constraint = 25, Moon_constraint = 30)
SPOCKLT.make_docx_schedule(schedule.observatory, schedule.telescope, schedule.date_range,'Manu',
                           schedule.target_list)

#SPOCKLT.make_np(schedule.date_range[0],schedule.date_range_in_days,schedule.telescope)
#SPOCKLT.upload_plans(schedule.date_range[0], nb_days=schedule.date_range_in_days,telescope = schedule.telescope)

print()


# ---------------------- SHORT TERM SCHEDULER ---------------------
obs = 1
schedule = SPOCKST.Schedules()
schedule.load_parameters('./input_short_term.csv',obs)

if schedule.use == 'follow_up':
    schedule.transit_follow_up('target_transit_follow_up.txt')
if schedule.use == 'special_start_end':
    input_name = 'HW_Vir'
    schedule.special_target_with_start_end(input_name)
if schedule.use == 'special':
    input_name = 'GJ-299'
    schedule.special_target(input_name)
if schedule.use == 'monitoring':
    input_name = 'Sp0755-2404'
    schedule.monitoring(input_name,5,61)
if schedule.use == 'dome_rotation':
    schedule.dome_rotation()

#schedule.make_scheduled_table()
#schedule.planification()
#schedule.make_night_block()


SPOCKST.make_plans(day=schedule.day_of_night,nb_days=1,telescope=schedule.telescope)



print()


subaru = Observer.at_site("Paranal")
time_range = Time(["2020-03-18 06:00", "2020-03-18 12:00"])

df = pd.read_csv('speculoos_target_list_v3.txt', delimiter=' ')
target_table_spc = Table.from_pandas(df)
targets = [FixedTarget(coord=SkyCoord(ra=target_table_spc['RA'][i] * u.degree, dec=target_table_spc['DEC'][i] * u.degree),
                name=target_table_spc['Sp_ID'][i]) for i in range(len(target_table_spc['RA']))]

# from astropy.io import ascii
# target_table = ascii.read('targets.txt')
#
# from astropy.coordinates import SkyCoord
# import astropy.units as u
# targets = [FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name)
#            for name, ra, dec in target_table]

constraints = [CelestialPoleSeparationConstraint(min=30*u.deg),
               AltitudeConstraint(min=23*u.deg, max=40*u.deg),MoonSeparationConstraint(min=30*u.deg)]
#observability = is_observable(constraints, subaru, targets,time_range=time_range)


paranal = Observer.at_site('Paranal')
# Define range of times to observe between
start_time = Time('2020-03-17 23:00:01')
end_time = Time('2020-03-18 11:00:01')
time_resolution = 1 * u.hour

# Create grid of times from ``start_time`` to ``end_time``
# with resolution ``time_resolution``
time_grid = time_grid_from_range([start_time, end_time],
                                 time_resolution=time_resolution)

observability_grid = np.zeros((len(constraints), len(time_grid)))

for i, constraint in enumerate(constraints):
    # Evaluate each constraint
    observability_grid[i, :] = constraint(paranal, targets[861], times=time_grid)

# Create plot showing observability of the target:

extent = [-0.5, -0.5+len(time_grid), -0.5, 2.5]

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


print()





schedule = SPOCKLT.Schedules()
# 1 for SSO , 2 for SNO and 3 for Saint-Ex
schedule.load_parameters('./input.csv',obs)

schedule.make_schedule(Altitude_constraint = 25, Moon_constraint = 30)

#SPOCKLT.save_schedule('/Users/elsaducrot/code/spock/input.csv',3,save=False,over_write=False)

#SPOCKLT.make_np(schedule.date_range[0],1,schedule.telescope)

#schedule.reference_table()

print()


schedule = SPOCKST.Schedules()
schedule.load_parameters('./input_short_term.csv',obs)
input_name = 'TI206544316'
#schedule.special_target(input_name)
schedule.special_target_with_start_end(input_name)
schedule.make_scheduled_table()
schedule.planification()
schedule.make_night_block()

print()
SPOCKST.make_np(schedule.day_of_night[0],1,schedule.telescope)






schedule = SPOCKLT.Schedules()
# 1 for SSO , 2 for SNO and 3 for Saint-Ex
schedule.load_parameters('/Users/elsaducrot/code/spock/input.csv',3)

SPOCKplot.airmass_plot_saved('SSO','Io',schedule.date_range[0])
SPOCKplot.airmass_altitude_plot_saved('SSO','Io',schedule.date_range[0])
#SPOCKplot.gantt_chart(schedule.date_range[0],schedule.date_range[1],['Io','Europa','Artemis'])

print()