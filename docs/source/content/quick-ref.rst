.. _quick-ref:

Quick reference
===============

Long term scheduling
--------------------

*SPOCK* has a module :mod:`SPOCK.long_term_scheduler` that contain all the  required
functions and classes to create appropriate plans (in the form of night blocks)for the SPECULOOS telescopes.

.. currentmodule:: SPOCK.long_term_scheduler

.. rubric:: Schedules class

.. autosummary::
   :nosignatures:

    Schedules


.. rubric:: Others functions

.. autosummary::
   :nosignatures:

    save_schedule
    make_plans
    upload_plans
    make_docx_schedule

Short term scheduling
----------------------

*SPOCK* also has a module :mod:`SPOCK.short_term_scheduler` that contain all the  required
functions and classes to add any external observations to existing plans. This include all
observations unrelated to SPECULOOS main programs.

.. currentmodule:: SPOCK.short_term_scheduler

.. rubric:: Schedules class

.. autosummary::
   :nosignatures:

    Schedules


.. rubric:: Others functions

.. autosummary::
   :nosignatures:

    save_schedule
    make_plans
    upload_plans
    make_docx_schedule

Plots
----------

*SPOCK* also has a module :mod:`SPOCK.plots_scheduler` that contain all the  required
functions to show visibility plots or gant chart, useful to keep track of scheduled observations.

.. currentmodule:: SPOCK.plots_scheduler

.. rubric:: Visibility plots

.. autosummary::
   :nosignatures:

    airmass_altitude_plot_given_target

.. rubric:: Gantt chart

.. autosummary::
   :nosignatures:

    gantt_chart_all
    gantt_chart

ACP files
----------

*SPOCK* also has modules :mod:`SPOCK.test_txtfiles` and :mod:`SPOCK.make_night_plans`  that contain all the  required
functions to convert the night blocks to ACP readable files.

.. currentmodule:: SPOCK.txt_files

.. rubric:: ACP files

.. autosummary::
   :nosignatures:

    startup
    target
    flatdawn
    flatexo_io
    biasdark
    shutdown


.. currentmodule:: SPOCK.make_night_plans

.. rubric:: Convert night blocks

.. autosummary::
   :nosignatures:

    make_np
