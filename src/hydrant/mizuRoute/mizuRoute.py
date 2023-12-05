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


def remapping_config (config):
    
    
    esmr = Easymore()

    # config from dictionary
    esmr = esmr.from_dict(config)
    
    # make sure the remapping file creation is correct
    esmr.only_create_remap_csv = True
    
    # execute EASYMORE
    esmr.nc_remapper()
    
    # read the remapping file and process
    ds = xr.open_dataset(config.get('temp_dir')+'_remapping.nc')
    
    # convert ds to df so the manupulation is easier
    df = ds.to_dataframe()
    #print(df)
    
    # sort based on ID_t
    df = df.sort_values(by=['ID_t'])
    
    # sort and get the case
    case = df['easymore_case'].values[0].item()
    #print(case)
    
    # river network ID and its frequency in intersection with runoff field grid
    RN_id, RN_frequency_in_intersection = np.unique(df['ID_t'].values, return_counts=True)
    print("The dimension of overlap between river network and hydrological units:")
    for dim in ds.dims:
        print(len(ds[dim]))
    if len(RN_id) == len(RN_frequency_in_intersection):
        print("The dimension of river network and its frequancy are:")
        print(len(RN_id))
        print(len(RN_frequency_in_intersection))
    else:
        sys.exit("the dimension of river network and its frequency are not the same")
        
    # Create a empty ntopo
    ntopo = xr.Dataset()
    
    # Populate the xarray dataarrays for polyid and their frequancy
    ntopo['polyid'] = xr.DataArray(RN_id, dims=('polyid',))
    ntopo['frequency'] = xr.DataArray(RN_frequency_in_intersection, dims=('polyid',))
    
    # populate the other varibales
    ntopo['IDmask'] = xr.DataArray(df['ID_t'].values, dims=('intersect',))
    ntopo['weight'] = xr.DataArray(df['weight'].values, dims=('intersect',))
    
    if case == 1 or case == 2:
        ntopo['cols'] = xr.DataArray(df['cols'].values+1, dims=('intersect',)) # python-EASYMORE to index for fortran
        ntopo['rows'] = xr.DataArray(df['rows'].values+1, dims=('intersect',)) # python-EASYMORE to index for fortran
    elif case == 3:
        ntopo['ID_s'] = xr.DataArray(df['ID_s'].values, dims=('intersect',))
        
    #print(ntopo)
    
    return ntopo

#     # get remapping information from easymore remap file
#     IDmask  = np.array(df['ID_t']) # the ID of river network
#     weight  = np.array(df['weight']) # the weight of each hydrological unit in river network
#     i_index = np.array(df['cols']) # cols for case 1 and 2
#     j_index = np.array(df['rows']) # rows for case 1 and 2
#     ID_s    = np.array(df['ID_s']) # ID from unstructure mesh case 3


#     # write the netcdf file
#     if os.path.isfile(file_to_save):
#         os.remove(file_to_save)

#     with nc4.Dataset(file_to_save, "w", format="NETCDF4") as ncid: # creating the NetCDF file

#         # define the dimensions
#         dimid_ID  = ncid.createDimension('polyid', len(RN_id))  # limited dimensiton equal the number of hruID
#         dimid_RMP = ncid.createDimension('intersect', len(df))   # intersection

#         # Variables RN_id, RN_frequency_in_intersection
#         #
#         RN_IDvar   = ncid.createVariable('polyid', 'int',  ('polyid', ), fill_value = -9999, zlib=True, complevel=9)
#         RN_IDvar.long_name = 'ID of River Network subbasins'
#         RN_IDvar.standard_name = 'ID of River Network subbasins'
#         RN_IDvar.units = '1'
#         RN_IDvar[:] = RN_id
#         #
#         RN_FRvar   = ncid.createVariable('frequency', 'int', ('polyid', ), fill_value = -9999, zlib=True, complevel=9)
#         RN_FRvar.long_name = 'Frequancy of intersection River Network subbasins with hydrological subbasins'
#         RN_FRvar.standard_name = 'Frequancy of intersection River Network subbasins with hydrological subbasins'
#         RN_FRvar.units = '1'
#         RN_FRvar[:] = RN_frequency_in_intersection

#         # ID_mask, weight, i_index, j_index, ID_s
#         #
#         IDmaskvar  = ncid.createVariable('IDmask', 'int',  ('intersect', ), fill_value = -9999, zlib=True, complevel=9)
#         IDmaskvar.long_name = 'ID of rive network subbasins'
#         IDmaskvar.standard_name = 'ID of rive network subbasins'
#         IDmaskvar.units = '1'
#         IDmaskvar[:] = IDmask
#         #
#         weightvar   = ncid.createVariable('weight', 'f8',  ('intersect', ), fill_value = -9999, zlib=True, complevel=9)
#         weightvar.long_name = 'weight of each hydrological unit in river network subbasins'
#         weightvar.standard_name = 'weight of each hydrological unit in river network subbasins'
#         weightvar.units = '1'
#         weightvar[:] = weight

#         if case == 1 or case == 2:
#             #
#             i_indexvar  = ncid.createVariable('i_index', 'int',  ('intersect', ), fill_value = -9999, zlib=True, complevel=9)
#             i_indexvar.long_name = 'cols from the source nc file'
#             i_indexvar.standard_name = 'cols from the source nc file'
#             i_indexvar.units = '1'
#             i_indexvar[:] = i_index + 1
#             #
#             j_indexvar   = ncid.createVariable('j_index', 'int',  ('intersect', ), fill_value = -9999, zlib=True, complevel=9)
#             j_indexvar.long_name = 'rows from the source nc file'
#             j_indexvar.standard_name = 'rows from the source nc file'
#             j_indexvar.units = '1'
#             j_indexvar[:] = j_index + 1

#         if case == 3:
#             #
#             ID_HRvar   = ncid.createVariable('ID_HR', 'int',  ('intersect', ), fill_value = -9999, zlib=True, complevel=9)
#             ID_HRvar.long_name = 'river network ID'
#             ID_HRvar.standard_name = 'river network ID'
#             ID_HRvar.units = '1'
#             ID_HRvar[:] = ID_s

#         # general attributes for NetCDF file
#         ncid.Conventions = 'CF-1.6'
#         ncid.Author = 'The data were written by easymore_codes'
#         ncid.License = 'MIT'
#         ncid.History = 'Created '
#         ncid.Source = 'Case: ; remapped by script from library of Shervan Gharari (https://github.com/ShervanGharari/EASYMORE).'


#     # show the remap file
#     ds = xr.open_dataset(file_to_save)
#     ds