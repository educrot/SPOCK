.. image:: https://img.shields.io/badge/docs-dev-blue.svg
    :target: https://educrot.github.io/SPOCK/index.html

.. image:: figures/logo_SPOCK_2.png
   :width: 600

**SPOCK** (Speculoos Observatory SChedule maKer) is a python package developed to handle
the planification of observation of the SPECULOOS telescopes. The project SPECULOOS -Search for habitable Planets EClipsing ULtra-cOOl Stars –
searches for potentially habitable exoplanets around the smallest and coolest stars
of the solar neighborhood `Link to site <https://www.speculoos.uliege.be/cms/c_4259452/fr/speculoos>`_.


Targets are scheduled following several criteria:

*  Observability

*  Priority

*  Number of hours already performed

Documentation SPOCK
---------------------

You will find complete documentation for setting up your project at `SPOCK Read
the Docs site`_.

.. _SPOCK Read the Docs site : https://educrot.github.io/SPOCK/index.html


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
Suggestions are welcome. For major changes, please open an issue first to discuss what you would like to change.


License
---------------------

<span style=“color:red;”> text </span>