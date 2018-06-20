"""
Scripts to read all csv files with wind data and compress them into Data frames that can be read by the python scripts
"""

import os
import sys
import numpy as np
import pandas as pd
import time
import pickle
import gc


def read_csv_wind(path_csv, list_files):

    #nfiles = 100
    #list_files = list_files[:nfiles]

    nfiles = len(list_files)
    int_print = 20

    start_time = time.time()
    t_last = start_time
    dffinal = None

    for i, f in enumerate(list_files):

        if i % int_print == 0:
            t_curr = time.time()

            delta = t_curr - start_time
            s1 = int(delta % 60)
            m1 = int(delta / 60)

            remaining = ((t_curr - t_last)/int_print)*(nfiles - i)
            s2 = int(remaining % 60)
            m2 = int(remaining / 60)

            print('{0:4d}. Processing {1}. Elapsed Time: {2:d}:{3:02d}. '
                  'Estimated remaining time: {4:d}:{5:02d}'.format(i, f, m1, s1, m2, s2))
            t_last = time.time()

        dftemp = pd.read_csv(filepath_or_buffer=os.path.join(path_csv, f))

        if i == 0:
            dffinal = dftemp
        else:
            dffinal = dffinal.merge(dftemp, on='Datetime', how='outer')

    end_time = time.time()

    delta = end_time - start_time
    s = int(delta % 60)
    m = int(delta / 60)

    print()
    print('Elapsed time: {0:d}:{1:02d}'.format(m, s))

    return dffinal


def read_csv_wind_save_pickle(path_csv, path_save=os.path.expanduser('~/'), size_batch=None):

    list_files = os.listdir(path_csv)
    nfiles = len(list_files)

    if size_batch is None:
        list_batches = [list_files]
    else:
        # split indexes into lists of approximately size size_batch
        list_idx = np.arange(nfiles)
        list_idx = np.array_split(list_idx, nfiles // size_batch)
        # create list with sublists of files of apporx. size size_batch
        list_batches = [[list_files[i] for i in idx] for idx in list_idx]

    for i, batch in enumerate(list_batches):
        a = read_csv_wind(path_csv=path_csv, list_files=batch)

        filename = os.path.join(path_save, 'wind_batch{0:02d}.pkl'.format(i))
        a.to_pickle(path=filename)

        print()
        print('Data frame saved to {}'.format(filename))
        print()

        # delete big data frame and run garbage collection (just in case)
        del a
        gc.collect()