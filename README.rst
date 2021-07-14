.. image:: https://travis-ci.com/educrot/SPOCK.svg?branch=master&status=passed
    :target: https://travis-ci.com/educrot/SPOCK

.. image:: https://img.shields.io/badge/docs-dev-green.svg
    :target: https://educrot.github.io/SPOCK/index.html

.. image:: ./SPOCK_Figures/logo_SPOCK_2.png
   :width: 600

**SPOCK** (Speculoos Observatory SChedule maKer) is a python package developed to handle the
planification of observation of the SPECULOOS telescopes. The project SPECULOOS -Search for habitable Planets EClipsing ULtra-cOOl Stars –
searches for potentially habitable exoplanets around the smallest and coolest stars
of the solar neighborhood `Link to site <https://www.speculoos.uliege.be/cms/c_4259452/fr/speculoos>`_.

*SPOCK* allows you to schedule SPECULOOS core program targets on several criteria:

*  Visibility of the target

*  Priority (calculated from stellar parameters)

*  Number of hours already performed

*  Coordination between different site

as well as external program targets (planetary candidates, eclipsing binaries, complex rotators, etc)

Documentation SPOCK
---------------------

You will find complete documentation (in dev) for setting up your project at `SPOCK Read
the Docs site <https://educrot.github.io/SPOCK/index.html>`_.


Installation
---------------------

.. _installation:


.. warning::
    You must be part of the SPECULOOS consortium to use *SPOCK*.


Install *SPOCK* locally::

    git clone https://github.com/educrot/SPOCK.git

    cd spock
    python setup.py install




More about *SPOCK*
---------------------

*SPOCK* is presented in more details in `Sebastian et al. 2020 <http://arxiv.org/abs/2011.02069>`_.
