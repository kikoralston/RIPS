.. _example:

Example
********

To run a CE, UCED, or combined CE/UCED simulation the function ``masterFunction()`` (see :ref:`main-scripts`) must be
called within a python script.

This python script must also load the parameters needed to execute the simulation and pass them as arguments to ``masterFunction()``.

A minimal example of how to execute a capacity expansion only analysis from 2015 to 2050 is illustrated below.

The txt files are formatted csv files that contain the parameters data.

The parameters are python objects in which the fields store the values of specific parameters. After loading the values of
the different parameters, they can also be changed by setting the values of the fields directly.::

   from RIPSMasterScript import *
   from Parameters import *

   # Load parameters
   # assumes that all parameter files are in the current working directory

   genparam = Generalparameters()
   genparam.load(fname='generalparameters.txt')

   reserveparam = Reserveparameters()
   reserveparam.load(fname='reserveparameters.txt')

   curtailparam = Curtailmentparameters()
   curtailparam.load(fname='curtailmentparameters.txt')

   # change some parameters fields

   # run only capacity expansion (CE) simulation
   genparam.runCE = True
   genparam.runFirstUCYear = False
   genparam.runUC = False

   # define start and end years of CE simulation
   genparam.startYear = 2015
   genparam.endYear = 2050

   masterFunction(genparam, reserveparam, curtailparam)


The results from a capacity expansion simulation are stored in the results folder (defined in the field ``resultsDir`` inside the
``Generalparameters`` object). Most results are saved as csv files. There are also `gdx` files which have all results
from the optimization model.

For example the file *genFleetAfterCE2050.csv* has the final composition of the generator fleet after the CE model simulation.