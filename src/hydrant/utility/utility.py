"""
Common functions for general tasks such as conversion of a csv file to nc for future
manupulaiton
"""


import xarray as xr
import numpy as np
import glob
import os
import sys

from typing import (
    Optional,
    Dict,
    Tuple,
    Union,
    Sequence,
)

import xarray as xr
import glob
import geopandas as gpd
import sys
from   easymore import Easymore
import pandas as pd
import pint_xarray



def mesh_output_txt_to_nc(csv_name,
                          data_frame_DateTime_column = 'time',
                          variable_name = 'data',
                          variable_dim_name = 'n',
                          unit_of_variable = ' ',
                          variable_long_name = ' ',
                          Fill_value = '-9999',
                          ddb_file = None,
                          rank_var_ddb = None,
                          segID_var_ddb = None,
                          nc_file_to_save = None):
    
    
    RFF_csv = pd.read_csv(csv_name, header=None)
    RFF_csv = RFF_csv.rename(columns={0: 'time'})
    RFF_csv['time'] = pd.to_datetime(RFF_csv['time'])
    RFF_csv = RFF_csv.set_index('time')
    RFF_csv = RFF_csv.iloc[:, :-1] # drop the last column which is empty
    

    # initializing EASYMORE object
    esmr = Easymore()

    # convert csv files of discharge, its flag and station info to netCDF
    Data = esmr.dataframe_to_netcdf_xr(RFF_csv,
                                       data_frame_DateTime_column = data_frame_DateTime_column,
                                       variable_name = variable_name,
                                       variable_dim_name = variable_dim_name,
                                       unit_of_variable = unit_of_variable,
                                       variable_long_name = variable_long_name,
                                       Fill_value = Fill_value)
    
    if not ddb_file is None:
        ddb = xr.open_dataset(ddb_file)
    
        if not rank_var_ddb is None:
            Data['Rank'] = xr.DataArray(ddb[rank_var_ddb].values, dims=(variable_dim_name,))
            
        if not segID_var_ddb is None:
            Data['segID'] = xr.DataArray(ddb[segID_var_ddb].values, dims=(variable_dim_name,))
            
    if not nc_file_to_save is None:
        if os.path.isfile(nc_file_to_save):
            os.remove(nc_file_to_save)
        Data.to_netcdf(nc_file_to_save)
    
    return Data