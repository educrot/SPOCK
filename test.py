import SPOCK.long_term_scheduler as SPOCKLT
import SPOCK.short_term_scheduler as SPOCKST
import SPOCK.plots_scheduler as SPOCKplot
import SPOCK.ETC as ETC


obs = 1
schedule = SPOCKLT.Schedules()
schedule.load_parameters('./input.csv',obs)
schedule.make_schedule(Altitude_constraint = 25, Moon_constraint = 30)



print()

# ---------------------- SHORT TERM SCHEDULER ---------------------

if schedule.use == 'follow_up':
    schedule.transit_follow_up('target_transit_follow_up.txt')
if schedule.use == 'special_start_end':
    input_name = 'Sp0755-2404'
    schedule.special_target_with_start_end(input_name)
if schedule.use == 'special':
    input_name = 'TI425933644'
    schedule.special_target(input_name)
if schedule.use == 'monitoring':
    input_name = 'Sp0755-2404'
    schedule.monitoring(input_name,5,61)

schedule.make_scheduled_table()
schedule.planification()
schedule.make_night_block()
SPOCKST.make_plans(day=schedule.day_of_night[0],nb_days=1,telescope=schedule.telescope)

# ---------------------- LONG TERM SCHEDULER ---------------------

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