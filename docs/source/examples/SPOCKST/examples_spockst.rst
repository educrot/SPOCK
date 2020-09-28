.. _examples_spockst:

Short term scheduling
----------------------------

For short term planification of any targets (SPECULOOS or other)  observations

Plans can be created with:

.. code:: ipython3

    import SPOCK.short_term_scheduler as SPOCKST

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

    schedule.make_scheduled_table()
    schedule.planification()
    schedule.make_night_block()

    SPOCKST.make_plans(day=schedule.day_of_night,nb_days=1,telescope=schedule.telescope)

.. note::

    To make the plans if you agree with what SPOCK output

.. code:: ipython3

    SPOCKST.make_plans(day=schedule.day_of_night,nb_days=1,telescope=schedule.telescope)
