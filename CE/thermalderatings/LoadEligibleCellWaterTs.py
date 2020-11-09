# Michael Craig, 17 May 2017
# Functions load cell water T data from cells eligible for new generators by first setting
# the folder names with future data for cells (each folder = 1 cell), then loading water temperatures
# from those folders.
#
# August 2018
# UW changed the format of the files with water data to NETCDF. I added a flag 'netcdf' to some functions to read them.


import os, copy, sys
import numpy as np
from AssignCellsToIPMZones import mapCellToIPMZone
from thermalderatings.ModifyGeneratorCapacityWithWaterTData import getGenToCellAndCellToGenDictionaries
from thermalderatings.AuxCurtailmentFuncs import (get_all_cells_from_netcdf, order_cells_by_flow, get_all_cells_in_zone,
                                                  loadCellWaterTs, createBaseFilenameToReadOrWrite)
import pandas as pd


def loadEligibleCellWaterTs(genFleet, currYear, genparam, curtailparam, netcdf=True, n_cells=100):
    """Load daily average water temperatures for grid cells for current simulation year

    This function selects the eligible grid cells that will be used in the analysis and loads the water data of current
    simulation year.

    Current options for which cells can host new plants:

    * 'All' (get all cell data)
    * 'WithGen' (get data for cells that have generator inside at beginning of CE runs).
    * 'max_flow' (get data for the first ``n_cells`` cells with max average flow in current year.).
    * 'macrocell' (get data for cells defined in a file 'macrocells.csv' that defines the centers of spatial clusters of cells. This was done to decrease computational burden.).

    This option must be defined in the object ``curtailparam`` in the field ``cellsEligibleForNewPlants``

    :param genFleet: (list) 2d list with generator fleet
    :param currYear: (integer) current year of simulation
    :param genparam: object of type :mod:`Generalparameters`
    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :param netcdf: (boolean) True if files are in netcdf format
    :param n_cells: (integer) number of cells to consider in 'max_flow' option
    :return: dict of {cell folder name : [[Datetime],[AverageWaterT(degC)], [AirT], [flow]]}
    """
    allCellFoldersInZone, eligibleCellFolders = setCellFolders(genFleet, currYear, genparam, curtailparam,
                                                               netcdf=netcdf, n_cells=n_cells)

    dict_out = dict()
    for gcm in curtailparam.listgcms:
        dict_out[gcm] = loadCellWaterTs(eligibleCellFolders, allCellFoldersInZone, curtailparam, gcm, currYear)

    return dict_out


def setCellFolders(genFleet, currYear, genparam, curtailparam, netcdf=True, n_cells=100):
    """Get list of eligible cells with future data in each ipm zone

    Current options for which cells can host new plants:

    * 'All' (get all cell data)
    * 'WithGen' (get data for cells that have generator inside at beginning of CE runs).
    * 'max_flow' (get data for the first ``n_cells`` cells with max average flow in current year.).
    * 'macrocell' (get data for cells defined in a file 'macrocells.csv' that defines the centers of spatial clusters of cells. This was done to decrease computational burden.).

    This option must be defined in the object ``curtailparam`` in the field ``cellsEligibleForNewPlants``

    :param genFleet: 2d list with generator fleet
    :param currYear: (integer) current year of simulation
    :param genparam: object of type :mod:`Generalparameters`
    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :param netcdf: (boolean) if true reads water data from NETCDF files (new format from August 2018)
    :param n_cells: (integer)  number of cells to consider in 'max_flow' option
    :return:
    """

    if netcdf:
        # get cells from first gcm in list
        gcm = curtailparam.listgcms[0]

        allCellFolders = get_all_cells_from_netcdf(curtailparam, gcm)
        allCellFoldersInZone = get_all_cells_in_zone(allCellFolders, genparam, curtailparam)

    else:
        allCellFolders = [name for name in os.listdir(curtailparam.rbmOutputDir)
                          if os.path.isdir(os.path.join(curtailparam.rbmOutputDir, name))]

        allCellFoldersInZone = isolateCellsInZones(allCellFolders, genparam)

    if genparam.cellsEligibleForNewPlants == 'all':
        eligibleCellFolders = copy.deepcopy(allCellFoldersInZone)

    elif genparam.cellsEligibleForNewPlants == 'maxflow':
        # get all cells in each zone ordered by flow (descending order)
        best_cells_zone_dict = order_cells_by_flow(genparam, curtailparam, currYear, sys.maxsize)

        # get list of cells with existing generators in it (including build decisions from previous CE runs)
        if genFleet is not None:
            cellLatLongToGenDict = getGenToCellAndCellToGenDictionaries(genFleet)[1]
            cells_list = list(cellLatLongToGenDict.keys())
            cells_list_str = list(map(lambda x: createBaseFilenameToReadOrWrite(curtailparam.locPrecision, x[0], x[1]),
                                      cells_list))
        else:
            cells_list_str = []

        # for each zone remove cells that already have a generator in it and keep only 'n_cells' best cells
        for z in best_cells_zone_dict:
            auxList = [x for x in best_cells_zone_dict[z] if x not in cells_list_str]
            auxList = auxList[:n_cells]
            best_cells_zone_dict[z] = auxList

        # combine into one single list with all eligible cells
        eligibleCellFolders = []
        for z in best_cells_zone_dict:
            eligibleCellFolders = eligibleCellFolders + best_cells_zone_dict[z]

    elif genparam.cellsEligibleForNewPlants == 'withGens':
        # (genToCellLatLongsDictValues,genToCellLatLongsDict,cellLatLongToGenDict,
        #             genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
        (genToCellLatLongsDict, cellLatLongToGenDict, genToCellLatLongsList) = \
            getGenToCellAndCellToGenDictionaries(genFleet)

        eligibleCellFolders = []
        for (cellLat, cellLong) in cellLatLongToGenDict:
            eligibleCellFolders.append(createBaseFilenameToReadOrWrite(curtailparam.locPrecision, cellLat, cellLong))

    elif genparam.cellsEligibleForNewPlants == 'macrocell':
        # uses definition of macro cells

        # read csv file previously compiled in R assigning cells to "macro cells" (larger group of neighboring cells)
        df_macro = pd.read_csv(os.path.join(curtailparam.rbmRootDir, 'macrocells.csv'))

        # create column with lat_lon name of grid cell
        df_macro['cell'] = df_macro.apply(lambda row: createBaseFilenameToReadOrWrite(curtailparam.locPrecision,
                                                                                      row['lat'], row['lon']), axis=1)

        best_cells_zone_dict = order_cells_by_flow(genparam, curtailparam, currYear, sys.maxsize, output_list=False)

        eligibleCellFolders = []

        df_total = None

        # for each zone, choose cell in each macro cell with greatest annual flow in current year
        for z in genparam.ipmZones:
            df_zone_1 = df_macro[df_macro['zone'] == z]
            df_zone_2 = best_cells_zone_dict[z]

            # merge both dfs
            df_zone = pd.merge(df_zone_1, df_zone_2, how='left', on='cell')

            df_zone = df_zone.groupby('allocated.cell').agg(np.max).reset_index()

            print('Number of candidate sites for {0}: {1:4d}'.format(z, df_zone.shape[0]))
            print(list(df_zone['cell']))

            if df_total is None:
                df_total = df_zone
            else:
                df_total = pd.concat([df_total, df_zone])

            # combine into one single list with all eligible cells among all zones
            eligibleCellFolders = eligibleCellFolders + list(df_zone['cell'])

    else:
        print('Parameter \'cellsEligibleForNewPlants\' must be either \'all\' or \'maxflow\' or \'withGens\' or '
              '\'macrocell\'....')
        print('Ending simulation!')
        sys.exit()

    return allCellFoldersInZone, eligibleCellFolders


def isolateCellsInZones(allCellFolders, genparam):
    # """Isolates ALL cells in zones of analysis and returns list of cell folders in zones
    #
    # THIS FUNCTION IS NOT USED ANYMORE WITH THE NETCDF DATA FORMAT. IT WAS USED WITH THE PREVIOUS VERSION OF THE UW DATA.
    #
    # :param allCellFolders:
    # :param genparam:
    # :return:
    # """

    allCellFoldersInZone = list()

    for cell in allCellFolders:
        cellZone = mapCellToIPMZone(cell, genparam.fipsToZones, genparam.fipsToPolys)
        if cellZone in genparam.impZones: allCellFoldersInZone.append(cell)

    return allCellFoldersInZone