import SPOCK.long_term_scheduler as SPOCKLT
import SPOCK.short_term_scheduler as SPOCKST
import SPOCK.plots_scheduler as SPOCKplot
from astropy.time import Time


schedule = SPOCKLT.schedules()
# 1 for SSO , 2 for SNO and 3 for Saint-Ex
schedule.load_parameters('/Users/elsaducrot/code/spock/input.csv',3)

SPOCKplot.airmass_plot_saved('SSO','Io',schedule.date_range[0])
SPOCKplot.airmass_altitude_plot_saved('SSO','Io',schedule.date_range[0])
#SPOCKplot.gantt_chart(schedule.date_range[0],schedule.date_range[1],['Io','Europa','Artemis'])

print()

# ---------------------- LONG TERM SCHEDULER ---------------------

schedule = SPOCKLT.schedules()
# 1 for SSO , 2 for SNO and 3 for Saint-Ex
schedule.load_parameters('/Users/elsaducrot/code/spock/input.csv',3)

schedule.make_schedule(Altitude_constraint = 25, Moon_constraint = 30)

#SPOCKLT.save_schedule('/Users/elsaducrot/code/spock/input.csv',3,save=False,over_write=False)

#SPOCKLT.make_np(schedule.date_range[0],1,schedule.telescope)

#schedule.reference_table()
day = Time('2019-11-13 15:00:00.000')
schedule.visibility_plot(day)

print()

# ---------------------- SHORT TERM SCHEDULER ---------------------

schedule = SPOCKST.schedules()
schedule.load_parameters('input_short_term.csv',1)

if schedule.use == 'follow_up':
    schedule.transit_follow_up('target_dynamical_follow_up.txt')
if schedule.use == 'special_start_end':
    input_name = 'Sp0755-2404'
    schedule.special_target_with_start_end(input_name)
if schedule.use == 'special':
    input_name = 'Sp0000-1245'
    schedule.special_target(input_name)
if schedule.use == 'monitoring':
    input_name = 'Sp0755-2404'
    schedule.monitoring(input_name,5,61)

schedule.make_scheduled_table()
schedule.planification()
schedule.make_night_block()

SPOCKST.make_np(schedule.day_of_night,1,schedule.telescope)

# schedule.save_schedule('input_short_term.csv',1,save=True,over_write=False)


print()
