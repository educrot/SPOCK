.. _examples_spocklt:

Long term  scheduling
----------------------------

Long term lanification of SPECULOOS targets observations

Plans can be created with:

.. code:: ipython3

    import SPOCK.long_term_scheduler as SPOCKLT

    schedule = SPOCKLT.schedules()
    obs = 1 # 1 for SSO , 2 for SNO and 3 for Saint-Ex
    schedule.load_parameters('./SPOCK/input.csv',obs)
    schedule.make_schedule(Altitude_constraint = 25, Moon_constraint = 30)

