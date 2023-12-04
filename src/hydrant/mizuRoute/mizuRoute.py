"""
Common functions for creating network topology for mizuRoute and read,
re-order mizuRoute output for analysis and visualzaition
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



def create_nc_ntopo(riv,
                    cat,
                    network): # can be hdma, or tdx
    
    if network.lower() == 'merit':
        
        # sort by COMID
        riv.sort_values(by='COMID').reset_index(drop=True)
        cat.sort_values(by='COMID').reset_index(drop=True)
        
        # check if the riv and cat IDs are similar
        if sum(riv['COMID'].values - cat['COMID'].values) != 0:
            sys.exit('The COMID of riv and cat should be the same')
        
        # rename lengthkm to length
        riv = riv.rename(columns={'lengthkm': 'length'})
        
        # drop geometry to make gdf to df
        ntopo = riv.drop(columns = 'geometry')
        
        # convert index to n
        ntopo = ntopo.rename_axis('n')
        
        # to xarray object,
        ntopo = ntopo.to_xarray()
        
        # add the units to the variables using pint
        ntopo = ntopo.pint.quantify({'length': 'km',
                                     'area': 'km**2',
                                     'uparea': 'km**2'})
        
        # convert dictionary
        convert = {'length': 'm',
                   'area': 'm**2',
                   'uparea': 'm**2'}
        
        # covert the units in ntopo
        ntopo = ntopo.pint.to(convert)
        ntopo = ntopo.pint.dequantify()
        
    elif network.lower() == 'hdma':
        
        # sort by COMID
        riv.sort_values(by='seg_id').reset_index(drop=True)
        cat.sort_values(by='hruid').reset_index(drop=True)
        
        # check if the riv and cat IDs are similar
        if sum(riv['seg_id'].values - cat['hruid'].values) != 0:
            sys.exit('The seg_id and hruid of riv and cat should be the same')
        
        # drop geometry to make gdf to df
        ntopo = riv.drop(columns = 'geometry')

        # convert index to n
        ntopo = ntopo.rename_axis('n')

        # to xarray object,
        ntopo = ntopo.to_xarray()
        
    # return the ntopo xarray object
    return ntopo
                    


def reorder_output(file_name,
                   order_ids,
                   var_id,
                   dim_id,
                   var_time,
                   dim_time,
                   var_to_keep = None,
                   save_reordered = False,
                   reorder_pre = 'reorder_',
                   output_folder = 'reorder'):
    
    """reorder a netcdf file based on an ID dim and variable

    Parameters
    ----------
    file_name : str
        The path to input forcing file(s). Can be specify with *.
    order_ids : Sequence[float]
        A sequence of ID that the nc will be reordered against.
    var_id : str
        Name of varibale that includes the original nc file IDs.
    dim_id : str
        Dimension of var_id.
    var_time : str
        Name of variable time.
    dim_time: Dict[str, str], optional
        Dimention of variable time.
    var_to_keep: list[str], optional
        collection of varibales to be included in reordered nc file(s).
    save_reordered: logical[], optional
        Flag to save the xarray object as nc file(s). If multiple file is 
        provided it will be set to True internally.
    reorder_pre: str, optional
        String that will added to the beginning of the saved reordered files
    output_folder: str, optional
        The path that the reorder file will be saved at.

    Returns
    -------
    xarray.Dataset
        Returns a reordered xarray.Dataset.
        If more files are identified, None.

    [FIXME]: The merge functionality could be added in the future if needed. 
    Currently the suggestion is to merge file using cdo mergetime outside of
    this function before or after reordering.
    """
    
    #
    files = sorted(glob.glob(file_name))
    
    if len(files)>1:
        print('The number of files passed to function is larger than 1. '+\
              'The function output will be set to None. '+\
              'The reorder files are going to be set written in '+\
              output_folder+' folder with prefix of '+reorder_pre)
        print('It is suggested to use packages such as xarray or cdo to '+\
              'merge the generated nc files if needed. Examples are: \n'+\
              'cdo mergetime input*.nc merged.nc # example of directly using cdo librarty \n'+\
              'or \n'+\
              'import cdo # binary required \n'+\
              'cdo_obj = cdo.Cdo()  # CDO object \n'+\
              'ds = cdo_obj.mergetime(input=input*.nc, returnXArray=variables)  # Mergeing')
        save_reordered = True

    for file in files:
              
        # open the nc file
        ds = xr.open_dataset(file)
        
        # drop unecessarily variables that are identified
        if not var_to_keep is None:
            variables_to_drop = [var for var in ds.variables if var not in var_to_keep]
            variables_to_drop = [elem for elem in variables_to_drop if elem not in [var_id,var_time]]
            print(variables_to_drop)
            ds = ds.drop_vars(variables_to_drop)
        
        # find the index in var_id and rearrange the nc file on that dimension
        idx = np.where(np.in1d(np.array(ds[var_id][:]), order_ids))[0]
        ds = ds.isel({dim_id:idx})
        
        # Save the rearranged dataset
        if save_reordered:
            if not os.path.isdir(output_folder):
                os.makedirs(output_folder)
            file_name_to_save = os.path.join(output_folder, reorder_pre+os.path.basename(file))
            if os.path.isfile(file_name_to_save):
                os.remove(file_name_to_save)
            ds.to_netcdf(file_name_to_save)
        
        # close the oppened file
        ds.close()
    
    if len(files)>1:
        ds = None
    
    # return
    return ds