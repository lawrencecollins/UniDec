from unidec_modules.ChromEng import *
from unidec_modules import unidectools as ud
import numpy as np
import pandas as pd 
import os

def get_data(eng, data_type, platemap):
    """Returns 2 DataFrames or list (if no platemap)containing x and y data of a defined datatype (e.g. massdat) for each well. Can also concatenate with a platemap"""

    datax, datay = [], []
    index = []

    for s in eng.data.spectra:
        temp = getattr(s, data_type)
        datax.append(temp[:, 0])
        datay.append(temp[:, 1])
        index.append(s.name) 

    x, y = pd.DataFrame(datax, index = index), pd.DataFrame(datay, index = index)

    return pd.concat([platemap, x], axis = 1), pd.concat([platemap, y], axis = 1)

# def to_hdf5(hdf5_name, dir, filenames, var1, var2, overwrite = True):
#     """Exports text files in directory to a single HDF5 file"""

#     if overwrite == True: # prevents errors
#         if os.path.exists(hdf5_name) == True:
#             os.remove(hdf5_name)

def get_window(data, thresh):
    """Returns list of (x1, x2) corresponding to window above intensity threshold (% of max y value where data is 2x2 array)"""

    x, y = data[:, 0], data[:, 1]
    thresh = y.max()*thresh

    filt = y>thresh

    windows = []
    temp = []
    for i, j in enumerate(filt):
        if j == True and filt[i-1] == False:
            temp.append(x[i])
        elif j == True and filt[i+1] == False:
            temp.append(x[i])
            windows.append(temp)
            temp = []

    return windows, thresh

