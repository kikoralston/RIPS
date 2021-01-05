Renewables
===================

This module has the functions used to compute available generation for wind and solar generators and potential energy for hydro generators.

Impacts of climate change on wind and solar generators was not simulated in this application. Wind and solar available generation is imported
from NREL databases (:cite:`USNationalRenewableEnergyLaboratory2010, Draxl2015`). Impacts of climate change on hydro generators was simulated
only for existing generators, since we assumed that there would be no new hydro generators built in the southeast U.S.


Wind and Solar
----------------

.. automodule:: renewables.GetRenewableCFs
    :members:

Hydro
------

.. automodule:: renewables.GetHydroMaxGenPotential
    :members:
