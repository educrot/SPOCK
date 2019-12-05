import SPOCK.long_term_scheduler as SPOCKLT
import SPOCK.short_term_scheduler as SPOCKST
import SPOCK.plots_scheduler as SPOCKplot
from astropy.time import Time

obs = 1

schedule = SPOCKST.schedules()
schedule.load_parameters('./input_short_term.csv',obs)

if schedule.use == 'follow_up':
    schedule.transit_follow_up('traget_list_special.txt')
if schedule.use == 'special_start_end':
    input_name = 'Sp0755-2404'
    schedule.special_target_with_start_end(input_name)
if schedule.use == 'special':
    input_name = 'TI206544316'
    schedule.special_target(input_name)
if schedule.use == 'monitoring':
    input_name = 'Sp0755-2404'
    schedule.monitoring(input_name,5,61)

schedule.make_scheduled_table()
schedule.planification()
schedule.make_night_block()
schedule.scheduled_table_sorted

print()
day = Time('2019-12-05 15:00:00.000')
SPOCKST.save_schedule('./input_short_term.csv',obs,save=True,over_write =True)

SPOCKST.make_np(day,1,schedule.telescope)
SPOCKST.upload_plans(day, nb_days=1,telescope = schedule.telescope)


# ---------------------- LONG TERM SCHEDULER ---------------------

schedule = SPOCKLT.schedules()
# 1 for SSO , 2 for SNO and 3 for Saint-Ex
schedule.load_parameters('/Users/elsaducrot/code/spock/input.csv',2)

schedule.make_schedule(Altitude_constraint = 25, Moon_constraint = 30)

#SPOCKLT.save_schedule('/Users/elsaducrot/code/spock/input.csv',3,save=False,over_write=False)

#SPOCKLT.make_np(schedule.date_range[0],1,schedule.telescope)

#schedule.reference_table()

print()


# ---------------------- SHORT TERM SCHEDULER ---------------------
obs = 1

schedule = SPOCKST.schedules()
schedule.load_parameters('./input_short_term.csv',obs)

if schedule.use == 'follow_up':
    schedule.transit_follow_up('traget_list_special.txt')
if schedule.use == 'special_start_end':
    input_name = 'Sp0755-2404'
    schedule.special_target_with_start_end(input_name)
if schedule.use == 'special':
    input_name = 'TIC_206544316'
    schedule.special_target(input_name)
if schedule.use == 'monitoring':
    input_name = 'Sp0755-2404'
    schedule.monitoring(input_name,5,61)

schedule.make_scheduled_table()
schedule.planification()
schedule.make_night_block()
schedule.scheduled_table_sorted

#SPOCKST.make_plans(day=schedule.day_of_night[0],nb_days=1,telescope=schedule.telescope)


# schedule.save_schedule('input_short_term.csv',1,save=True,over_write=False)


print()



schedule = SPOCKLT.schedules()
# 1 for SSO , 2 for SNO and 3 for Saint-Ex
schedule.load_parameters('/Users/elsaducrot/code/spock/input.csv',3)

SPOCKplot.airmass_plot_saved('SSO','Io',schedule.date_range[0])
SPOCKplot.airmass_altitude_plot_saved('SSO','Io',schedule.date_range[0])
#SPOCKplot.gantt_chart(schedule.date_range[0],schedule.date_range[1],['Io','Europa','Artemis'])

print()