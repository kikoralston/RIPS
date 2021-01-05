.. test_documentation documentation master file, created by
   sphinx-quickstart on Thu Apr 16 11:34:36 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:github_url: https://github.com/kikoralston/RIPS


Capacity expansion / unit commitment and economic dispatch under climate change
********************************************************************************

.. toctree::
   :hidden:
   :maxdepth: 2

   example
   Input data/index
   GAMS/index
   Parameters/index
   Renewables/index
   Demand/index
   Setup Fleet/index
   Thermal deratings/index
   Main Scripts/index


Introduction
=============

The capacity expansion / unit commitment and economic dispatch under climate change model uses an integrated,
multi-model framework to examine how climate-related risks in the power sector could affect power system's planning and
operation decisions. The figure below illustrates our modeling framework.

.. figure::  ./_static/diagram_models.png
   :align:   center

   *Modeling framework*


We used data sets with projections of high-resolution (1/8th degree, ~12 km) daily surface meteorological variables
from twenty different Global Circulation Models (GCM) from the Coupled Model Intercomparison Project 5 :cite:`Taylor:2012jg`,
spatially downscaled using the Multivariate Adaptive Constructed Analogs (MACA) method :cite:`Abatzoglou:2012kc`.
These projections of climate forcings were used by an econometric model to simulate future regional hourly electricity
demand under climate change scenarios :cite:`RalstonFonseca2019`. Also, the data from the GCMs was used to simulate
daily river flows and water temperatures in our study region :cite:`Cheng2020`. To simulate these daily river flows
and water temperatures, we combined a macro-scale hydrological model, the Variable Infiltration Capacity model (VIC)
:cite:`Liang2004`; a coupled routing and water management model, the Model for Scale Adaptive River Transport
Large-Scale Water Management Model (MOSART-WM) :cite:`Li2013,Voisin2013, Voisin2017`; and a stream temperature model,
the River Basin Model (RBM) :cite:`Yearsley2009`. Next, the generated river flows
and water temperatures were used to estimate hydropower potential from existing hydro
generators. To simulate how thermoelectric generators are affected by climate change, we used the
Integrated Environmental Control Model (IECM) :cite:`IECM:2018, Zhai2011`.
The IECM outputs plant performance characteristics and costs for different combinations of power plant technologies and
cooling systems. We used this model to estimate typical operating response curves of thermal generators to changes in
climate and hydrological conditions :cite:`Loew2020`. Then, we used these operating response curves with the
projections of climate forcings and hydrological variables to simulate operating conditions of thermoelectric power
plants under climate change.  Additionally, we used the data set of simulated wind and solar generation
profiles from the U.S. National Renewable Energy Laboratory (NREL) :cite:`USNationalRenewableEnergyLaboratory2010, Draxl2015`
to define upper bounds on wind and solar power plants hourly generation in our model. We combined all these simulated data
sets into a capacity expansion model to simulate future generator additions :cite:`RalstonFonseca2020`.
Also, we integrated the resulting generator fleets in a unit commitment and economic dispatch model to simulate these fleet
configurations under consistent climate change conditions.


To run a CE, UCED, or combined CE/UCED simulation the function ``masterFunction()`` (see :ref:`main-scripts`) must be called within a python script.
This python script must also load the parameters needed to execute the simulation and pass them as arguments to ``masterFunction()``.
A `GAMS <https://www.gams.com/>`_ license is necessary to run this model.

See :ref:`example` for a minimal example of how to run the model.

..
   This is q inline equation If :math:`\sigma_{1}` equals :math:`\sigma_{2}` then etc, etc.

   This is a normal equation

   .. math::
      \alpha+\beta=\gamma


Acknowledgments
==================

Funding for this work came from the National Science Foundation (NSF) as part of the Resilient Interdependent
Infrastructure Processes and Systems (RIPS) program via Grant Number EFRI-1441131.


References
==================

.. bibliography::
   :style: plain
   :cited:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
