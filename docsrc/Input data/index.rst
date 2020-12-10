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


* Curtailment: this folder contains the json files with the parameters used to simulate climate induced capacity deratings in thermal generators
* DemandData: this folder contains the csv files with parameters of the regression function used to simulate electricity demand
* eGrid: this folder contains csv files extracted from the eGrid model that are used to create the initial fleet
* EIA860: this folder contains csv files extracted from the EIA860 with cooling data that are used to create the initial fleet
* FuelPricesCapacityExpansion: this folder contains the csv file with fuel prices projections
* GAMS: this folder contains the GAMS files with the optimization models
* HydroMonthlyDataPNNL: this folder contains pre-processed pickle files with hydro potential data for existing hydro plants
* NEEDS: this folder contains csv files extracted from the NEEDS model that are used to create the initial fleet
* NewPlantData: this folder contains csv files with the parameters of new technologies used in the CE model
* NRELSolarPVData: this folder contains csv files with solar potential data (source: NREL)
* PHORUM: this folder contains csv files with PHORUM data used in the UCED model
* StateShapeFiles: this folder contains shapefiles with the
* Transmission: this folder contains a csv file with transmission line capacity
* WINDSERCData: this folder contains csv files with wind potential data (source: NREL)
* cells2zones.pk: this pickle file contains a pre-processed dictionary that maps each cell to its IPM zone
* co2values.json: this json file defines the configuration of the CO2 cap (see function readInitialCondCO2)


Climate and hydrology data
============================

Because of their potential large size, climate and hydrological data can be stored in a separate folder. This directory must be specified in the :mod:`Curtailmentparameters` object. 
