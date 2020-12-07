import SPOCK.long_term_scheduler as SPOCKLT
import SPOCK.short_term_scheduler as SPOCKST
import SPOCK.plots_scheduler as SPOCKplot
import SPOCK.make_night_plans as SPOCKmkpl
import SPOCK.stats as SPOCKstats
from astropy import units as u
from astropy.time import Time
import SPOCK.ETC as ETC

obs = 4
schedule = SPOCKST.Schedules()
schedule.load_parameters('./input.csv',obs)
SPOCKST.make_docx_schedule(schedule.observatory, schedule.telescope, schedule.date_range,'Manu')


print()

obs = 1
schedule = SPOCKST.Schedules()
schedule.load_parameters('./input_short_term.csv',obs)
schedule.use = 'special'
# schedule.start_end_range = Time(['2020-12-13 15:00:00','2020-12-15 15:00:00'])
schedule.day_of_night = Time('2020-12-12 15:00:00')
schedule.telescope = 'Europa'
schedule.special_target(input_name="Trappist-1")
schedule.make_scheduled_table()
schedule.planification()
schedule.make_night_block()

print()

# ---------------------- SHORT TERM SCHEDULER ---------------------
obs = 1
schedule = SPOCKST.Schedules()
schedule.load_parameters('./input_short_term.csv',obs)

if schedule.use == 'follow_up':
    schedule.transit_follow_up('./target_lists/target_transit_follow_up.txt',name="Trappist-1b")
if schedule.use == 'special_start_end':
    input_name = 'TOI-1847'
    schedule.special_target_with_start_end(input_name)
if schedule.use == 'special':
    input_name = 'Trappist-1'
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

obs = 5
schedule = SPOCKLT.Schedules()
schedule.load_parameters('./input.csv',obs)


#SPOCKLT.make_np(schedule.date_range[0], schedule.date_range_in_days, schedule.telescope)
#SPOCKLT.upload_plans(schedule.date_range[0], nb_days=schedule.date_range_in_days, telescope=schedule.telescope)

#schedule.exposure_time_table(day=None)

schedule.make_schedule(Altitude_constraint = 28, Moon_constraint = 30)

#SPOCKLT.make_docx_schedule(schedule.observatory, schedule.telescope, schedule.date_range,'Manu',
    #                       schedule.target_list)


#SPOCKLT.make_np(schedule.date_range[0],schedule.date_range_in_days,schedule.telescope)
#SPOCKLT.upload_plans(schedule.date_range[0], nb_days=schedule.date_range_in_days,telescope = schedule.telescope)

print()

# ---------------------- Short scripts ---------------------

#SPOCKplot.phase_coverage_given_target(name_observatory='SSO', target='Sp0004-2058',
                                    #  path_target_list='speculoos_target_list_v6.txt', pmin=0, pmax=4)

# a = SPOCKstats.read_night_plans_server(telescope="Io",date="20201006")

#SPOCKLT.upload_plans('2020-11-01 00:00:00', nb_days=30,telescope = 'Europa')

#SPOCKmkpl.make_np('2020-11-01 00:00:00',30,'Europa')


schedule = SPOCKLT.Schedules()
schedule.load_parameters('./input.csv',1)
schedule.exposure_time_table(day=None)

print()






# ---------------------- LONG TERM SCHEDULER ---------------------

obs = 2
schedule = SPOCKLT.Schedules()
schedule.load_parameters('./input.csv',obs)


#SPOCKLT.make_np(schedule.date_range[0], schedule.date_range_in_days, schedule.telescope)
#SPOCKLT.upload_plans(schedule.date_range[0], nb_days=schedule.date_range_in_days, telescope=schedule.telescope)

#schedule.exposure_time_table(day=None)

schedule.make_schedule(Altitude_constraint = 28, Moon_constraint = 30)

#SPOCKLT.make_docx_schedule(schedule.observatory, schedule.telescope, schedule.date_range,'Manu',
    #                       schedule.target_list)


#SPOCKLT.make_np(schedule.date_range[0],schedule.date_range_in_days,schedule.telescope)
#SPOCKLT.upload_plans(schedule.date_range[0], nb_days=schedule.date_range_in_days,telescope = schedule.telescope)

print()

# ---------------------- Plots ---------------------

SPOCKplot.phase_coverage_given_target(name_observatory='SSO',target='Sp0004-2058',path_target_list=None,pmin=0.1,pmax=6)

print()

