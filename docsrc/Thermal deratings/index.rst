.. _thermal-deratings:

Thermal Deratings
===================

The thermal deratings modules includes all functions used to simulate capacity deratings of thermal generators (existing or hypothetical).

This modules uses data from climate models (in netcdf format); hydrological simulations from the VIC, RBM, and MOSART-WM models (netcdf format);
and the pre-defined parameters of the operating response functions of thermal generators to ambient conditions (in json format).

While these functions were implemented within the CE/UCED simulation models, they can be called separately from the main model in order
to just compute spatially- and temporally diferentialted thermal deratings of existing or hypothetical thermal power plants.

.. automodule:: thermalderatings.LoadEligibleCellWaterTs
    :members:

