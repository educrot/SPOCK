**SPOCK** (Speculoos Observatory SChedule maKer) is a Python library for dealing with the planification of SPECULOOS targets observations

Schedule targets on several criteria:

*  Visibility of the target

*  Priority (from different metrics)

*  number of hours already performed

Documentation for RTD
---------------------

You will find complete documentation for setting up your project at `SPOCK Read
the Docs site`_.

.. _SPOCK Read the Docs site: http://speculoos7.astro.ulg.ac.be/gitlab/-/ide/project/eDucrot/spock/edit/master/-/docs/source/index.rst

Installation
---------------------

Use the package manager [git clone]() to install SPOCK::

    git clone http://speculoos7.astro.ulg.ac.be/gitlab/eDucrot/spock.git

    cd spock



Input files
---------------------

For ``long_term_scheduler`` reate your *'input_file.csv'* file in the following format::


    date_range: 
      - "2020-05-11 15:00:00"
      - "2020-05-31 15:00:00"
    observatories:
      1:
        name: SSO
        telescopes: [Callisto,Europa,Io]
      2:
        name: SNO
        telescopes: [Artemis]
      3: 
        name: Saint-Ex
        telescopes: [Saint-Ex]
      4: 
        name: TS_La_Silla
        telescopes: [TS_La_Silla]
      5: 
        name: TN_Oukaimeden
        telescopes: [TN_Oukaimeden]
    strategy: "continuous"
    duration_segments: 20
    nb_segments: 3
    target_list: speculoos_target_list_v6.txt


For ``short_term_scheduler`` create your *'input_file.csv'* file in the following format::

    day_of_night: `
      - "2019-11-20 15:00:00"
    start_end_range: 
      - "2019-11-21 04:00:00"
      - "2019-11-21 10:30:00"
    use: "follow_up"
    observatories:
      1:
        name: SSO
        telescopes: [Io,Ganymede,Callisto,Io,Europa]
      2:
        name: SNO
        telescopes: [Artemis]
      3: 
        name: Saint-Ex
        telescopes: [Saint-Ex]
      4: 
        name: TS_La_Silla
        telescopes: [TS_La_Silla]
      5: 
        name: TN_Oukaimeden
        telescopes: [TN_Oukaimeden]
    target_list: target_list_special.txt



Contributing
---------------------
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

License
---------------------

<span style=“color:red;”> text </span>