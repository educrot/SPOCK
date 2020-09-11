import SPOCK.long_term_scheduler as SPOCKLT
import SPOCK.short_term_scheduler as SPOCKST
import SPOCK.plots_scheduler as SPOCKplot


# ---------------------- SHORT TERM SCHEDULER ---------------------
obs = 2
schedule = SPOCKST.Schedules()
schedule.load_parameters('./SPOCK/input_short_term.csv',obs)

if schedule.use == 'follow_up':
    schedule.transit_follow_up('./SPOCK/target_transit_follow_up.txt')
if schedule.use == 'special_start_end':
    input_name = 'HW_Vir'
    schedule.special_target_with_start_end(input_name)
if schedule.use == 'special':
    input_name = 'Sp0025+5422'
    schedule.special_target(input_name)
if schedule.use == 'monitoring':
    input_name = 'Sp0755-2404'
    schedule.monitoring(input_name,5,61)
if schedule.use == 'dome_rotation':
    schedule.dome_rotation()

schedule.make_scheduled_table()
schedule.planification()
schedule.make_night_block()


SPOCKST.make_plans(day=schedule.day_of_night,nb_days=1,telescope=schedule.telescope)

print()

# ---------------------- LONG TERM SCHEDULER ---------------------

obs = 1
schedule = SPOCKLT.Schedules()
schedule.load_parameters('./SPOCK/input.csv',obs)

schedule.make_schedule(Altitude_constraint = 24, Moon_constraint = 30)

#SPOCKLT.make_docx_schedule(schedule.observatory, schedule.telescope, schedule.date_range,'Manu',
#                           schedule.target_list)

#SPOCKLT.make_np(schedule.date_range[0],schedule.date_range_in_days,schedule.telescope)

print()