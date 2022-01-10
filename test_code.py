import SPOCK.short_term_scheduler as SPOCKST
from astropy.time import Time
import datetime

Target = "Sp0227-6447.01"
nb_transits = 10


schedules_st = SPOCKST.Schedules()
schedules_st.load_parameters()

day = Time("2021-11-25 15:00:00").iso
ra, dec, timing, period, duration = SPOCKST.get_info_follow_up_target(Target,
                                                                      target_list_follow_up=schedules_st.target_table_spc_follow_up)
table_prediction = SPOCKST.prediction(name=Target,ra=ra,dec=dec,
                   timing=timing,
                   period=period,duration=duration,
                   start_date=day,ntr=nb_transits)
print(table_prediction)




