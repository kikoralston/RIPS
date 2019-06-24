"""
This script adds FIPS codes and IPM zone names to the metadata csv files for wind and solar

This way this assignment will not have to be recomputed everytime by the function getPlantInfoInZone
"""

import pandas as pd
import numpy as np
import os
from Parameters import *
from AuxFuncs import str_elapsedtime
from AssignCellsToIPMZones import getFIPSOfPt
from GetRenewableCFs import getCoordsFromFilename


# auxiliary functions
def f(row):
    return getFIPSOfPt(fipsToPolys, row['latitude'], row['longitude'])


def g(row):
    if row['fips'] in fipsToZones:
        zone = fipsToZones[row['fips']]
    else:
        zone = None  # several lat/lons on coast are not included in IPM zones

    return zone


cwd = os.getcwd()

genparam = Generalparameters.Generalparameters()
genparam.load(fname=os.path.join(cwd, 'generalparameters.txt'))

fipsToZones, fipsToPolys = genparam.fipsToZones, genparam.fipsToPolys

for typerenew in ['wind', 'solar']:

    t_step = time.time()
    print(typerenew)

    if typerenew.lower() == 'wind':
        renewableDir = os.path.join(genparam.dataRoot, 'WINDSERCData')
        path_metadata_old = os.path.join(renewableDir, 'toolkit_sites_v7_SERC.csv')
        path_metadata_new = path_metadata_old.replace('.csv', '_zones.csv')
        metadata = pd.read_csv(path_metadata_old, dtype={'site_id': np.int, 'longitude': np.float,
                                                         'latitude': np.float, 'State': object,
                                                         'County': object, 'fraction_of_usable_area': np.float,
                                                         'power_curve': np.float, 'capacity': np.float,
                                                         'wind_speed': np.float,
                                                         'capacity_factor': np.float})
    else:
        renewableDir = os.path.join(genparam.dataRoot, 'NRELSolarPVData', 'SERC')

        path_metadata_old = os.path.join(renewableDir, 'SolarCapacityFactorsNRELSERC.csv')
        path_metadata_new = path_metadata_old.replace('.csv', '_zones.csv')
        metadata = pd.read_csv(path_metadata_old, dtype={'State': object, 'File': object, 'CF': np.float,
                                                         'PlantSize': np.float})

        metadata['latitude'] = metadata['File'].apply(lambda x: float(getCoordsFromFilename(x)[0]))
        metadata['longitude'] = metadata['File'].apply(lambda x: float(getCoordsFromFilename(x)[1]))

    metadata['fips'] = metadata.apply(f, axis=1)
    metadata['zone'] = metadata.apply(g, axis=1)

    metadata.to_csv(path_metadata_new)

    print('Done! Elapsed time: {}'.format(str_elapsedtime(t_step)))
    print()
