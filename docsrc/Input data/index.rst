Input data
**********

The simulation uses data from different sources. Input data is located in two different directories: 

* a directory for general input data;
* a directory for climate and hydrology data (which because of its large size may need to be stored in a different location)

General Input data
=====================

General input data is located in a directory with a specific subfolder structure. The path for this directory must be specified in the field `dataRoot` of the
:mod:`Generalparameters` object.

The diagram below shows the structure of the subfolders in the directory.

::

   dataRoot/
      |--Curtailment/
      |--DemandData/
      |--eGrid/
      |--EIA860/
      |--FuelPricesCapacityExpansion/
      |--GAMS/
      |--HydroMonthlyDataPNNL/
      |--NEEDS/
      |--NewPlantData/
      |--NRELSolarPVData/
      |--PHORUM/
      |--StateShapeFiles/
      |--Transmission/
      |--WINDSERCData/
      |--cells2zones.pk
      |--co2values.json

        
Climate and hydrology data
============================

Because of their potential large size, climate and hydrological data can be stored in a separate folder. This directory must be specified in the :mod:`Curtailmentparameters` object. 
