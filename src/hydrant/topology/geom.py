"""
Common functions for working with river and sub-basin geometries
in the context of hydrological modelling
"""

import geopandas as gpd
import pandas as pd
import os

from typing import (
    Optional,
    Dict,
    Tuple,
    Union,
    Sequence,
)

from .river_graph import find_upstream

def merit_read_file (
    pfafs: list,
    path_riv: str,
    riv_file_template: str,
    path_cat: str,
    cat_file_template: str,
    path_cst: Optional [str] = None,
    cst_file_template: Optional [str] = None):
    
    cat_all = gpd.GeoDataFrame()
    riv_all = gpd.GeoDataFrame()
    
    for pafaf in pfafs:
        
        # read files cat, riv, cst
        riv = gpd.read_file(os.path.join(path_riv, riv_file_template.replace('*', pafaf)))
        cat = gpd.read_file(os.path.join(path_cat, cat_file_template.replace('*', pafaf)))
        if not path_cst is None and not cst_file_template is None:
            cst = gpd.read_file(os.path.join(path_cst, cst_file_template.replace('*', pafaf)))
            # add cat and cst
            cst = merit_cst_prepare(cst,
                                    {'id':'COMID','area':'unitarea'},
                                    cat = cat,
                                    cat_col_id = 'COMID')
        else:
            cst = None
        
        # merge the cat and cst
        cat = merit_cat_cst_merge (cat,
                                   cst = cst)
        
        # append the files
        riv_all = pd.concat([riv_all, riv])
        cat_all = pd.concat([cat_all, cat])
    
    # sort COMID
    riv_all.sort_values(by='COMID', axis='index', inplace=True)
    riv_all.reset_index(drop=True, inplace=True)
    
    # sort COMID
    cat_all.sort_values(by='COMID', axis='index', inplace=True)
    cat_all.reset_index(drop=True, inplace=True)
    
    # set the projection
    riv_all.set_crs(epsg=4326, inplace=True, allow_override=True)
    cat_all.set_crs(epsg=4326, inplace=True, allow_override=True)
    
    return riv_all, cat_all
    

def merit_cst_prepare(
    cst: gpd.GeoDataFrame,
    cst_col: Optional[Dict[str, str]] = None,
    cat: Optional[gpd.GeoDataFrame] = None,
    cat_col_id: Optional[str] = None,
    cst_col_id_reset: bool = True,
    crs: int = 4326,
    *args,
    ) -> gpd.GeoDataFrame:

    # get the possible existing id, area if exists
    cst_col_id = 'COMID'
    cst_col_area = 'unitarea'
    if cst_col is not None:
        cst_col_id = cst_col.get('id')
        cst_col_area = cst_col.get('area')

    if not cst.crs:
        cst.set_crs(epsg=4326, inplace=True)
        warnings.warn('CRS of the coastal hillslope Shapefile has been assumed to be EPSG:4326')

    if cst_col_id_reset:
        max_cat_id = 0
        if cat is not None:
            max_cat_id = cat[cat_col_id].max()
        cst[cst_col_id] = range(max_cat_id+1,
                                max_cat_id+1+len(cst))
    else:
        if not cst_col_id in cst.columns:
            sys.exit('the corresponding id is not given for cosatl hillslope')
        else:
            max_cat_id = 0
            if cat is not None:
                max_cat_id = cat[cat_col_id].max()
            min_cst_id = cst[cst_col_id].min()
            if min_cst_id < max_cat_id:
                sys.exit('there is some mixed up COMID between the cat and costal hillslope')

    if not cst_col_area in cst.columns: # then we need to populate the id
        cst[cst_col_area] = cst.to_crs(epsg=6933).area / 1e6

    # assign `coastal_hillslope` flag
    cst['hillslope'] = 1
    
    # drop FID column
    cst = cst.drop(columns = ['FID'])

    # return
    return cst

def merit_cat_cst_merge (cat,
                         cst: Optional[gpd.GeoDataFrame] = None,
                         crs: int = 4326) -> gpd.GeoDataFrame:
    
    cat['hillslope'] = 0
    if not cst is None:
        cst['hillslope'] = 1
        catchment = gpd.GeoDataFrame(pd.concat([cat, cst]))
    else:
        catchment = cat

    # assign crs
    catchment.set_crs(epsg=crs, inplace=True)
    catchment.reset_index(drop=True, inplace=True)
    
    return catchment
    

def prepare_ntopo(
    riv: gpd.GeoDataFrame,
    riv_cols: Dict[str, str],
    cat: gpd.GeoDataFrame,
    cat_cols: Dict[str, str],
    network: str,
    outlet_val: int = -9999,
    crs: int = 4326,
    *args,
    **kwargs,
) -> gpd.GeoDataFrame:
    '''prepares and match the river with costal and other hillslope which are presented in cat
    but not presented in riv. It it so make sure the routing model has the needed variables
    Parameters
    ----------
    riv: geopandas.GeoDataFrame
        The river network object containing river segments
    riv_cols: dict
        A dictionary containing 'id', 'next_id', keys and values
        corresponding to the river segement IDs and their downstream
        river segments, respectively. Also, other keys includes 'slope',
        'length', 'length_direct'
    cat: geopandas.GeoDataFrame
        The catchment topology object derived from `prepare_cat`
        function
    cat_cols: dict
        A dictionary containing 'id', 'hillslope' keys and values
        corresponding to the catchment IDs and also the flag indicating
        whether the catchment is a coastal hillslope one (1) or not (0)
    outlet_val: int, optional [defaults to -9999]
        Value being assigned to the 'next_id' cell of the each
        river segments indicating the segment is an outlet for the network;
        before any subsetting, if any coastal hillslope segments exist, they
        will only be set to the outlet value.
    Returns
    -------
    river: geopandas.GeoDataFrame
        The prepared river network topology GeoDataFrame object
    '''

    # necessary assertions
    assert isinstance(riv, gpd.GeoDataFrame), "`riv` object must be of type geopandas.GeoDataFrame"
    assert isinstance(riv_cols, dict), "`riv_cols` object must be of type dict"
    assert isinstance(cat, gpd.GeoDataFrame), "`cat` object must be of type geopandas.GeoDataFrame"
    assert isinstance(cat_cols, dict), "`cat_cols` object must be of type dict"

    # `id` and `next_id` keys must be included in the `riv_cols` dictionary
    if not 'id' and 'next_id' in riv_cols:
        raise ValueError("`id` and `next_id` must be included in the `riv_cols` dictionary")

    # `id`, `geom`, and `hillslope` keys must be included in the `cat_cols` dictionary
    if not 'id' and 'geom' and 'hillslope' in cat_cols:
        raise ValueError("`id`, `geom` and `hillslope` must be included in the `cat_cols` dictionary")

    # necessary definitions
    riv_col_id = riv_cols.get('id')
    riv_col_next_id = riv_cols.get('next_id')
    hil_col_id = cat_cols.get('hillslope')
    cat_col_id = cat_cols.get('id')
    cat_col_area = cat_cols.get('area')

    # extract the hillslope column from `cat` and assign `id` as index
    hills = cat.set_index(cat_col_id)
    hills = hills[hil_col_id]

    # add the hillslope flag from `cat`
    # `cat` has necessary greater than or equal number of items than `riv`
    river = gpd.GeoDataFrame(riv.join(other=hills, on=riv_col_id, how='outer'),
                             geometry='geometry')

    # assign `outlet_val` to the `hillslope_id` river segments
    river.loc[river[hil_col_id] == 1, riv_col_next_id] = outlet_val

    # sort based on `riv_col_id` and reset_index
    river.sort_values(by=riv_col_id, axis='index', inplace=True)
    river.reset_index(drop=True, inplace=True)
    
    # sort based on `cat_col_id` and reset_index
    cat.sort_values(by=cat_col_id, axis='index', inplace=True)
    cat.reset_index(drop=True, inplace=True)

    # assign crs
    river.set_crs(epsg=crs, inplace=True)
    
    
    if network == 'merit':
        
        # fix cicular segemnts, in other words, those whose 'id' and
        # 'next_id's are the same
        river.loc[river[riv_col_id] == river[riv_col_next_id], riv_col_next_id] = outlet_val
        
        # pass the unit area from cat
        river['unitarea'] = 0
        river['unitarea'] = cat['unitarea']
        
        # fill NAs with `riv_na_val` values in `rivers`
        river.loc[river[hil_col_id] == 1, ['maxup','up1','up2','up3','up4']] = 0
        river.loc[river[hil_col_id] == 1, ['lengthkm','lengthdir','sinuosity','slope','slope_taud']] = 0.001
        river.loc[river[hil_col_id] == 1, ['strmDrop_t','order']] = 1
        river.loc[river[hil_col_id] == 1, 'uparea'] = river.loc[river[hil_col_id] == 1, 'unitarea']
        
        # approximate width based on uparea
        
        # dictionary
        data_types = {
            'COMID': int,       
            'order': int,     
            'NextDownID': int,       
            'maxup': int,
            'up1': int,
            'up2': int,
            'up3': int,
            'up4': int,
            'hillslope': int,
        }

        # Convert columns to specified data types
        river = river.astype(data_types)

        
#         # necessary columns being chosen for downcasting of dtypes
#         # and other necessary assignments
#         col_list = list(river.columns)
#         for col in [cat_col_id, cat_col_geom]:
#             col_list.remove(col)

#         # fill NAs with `riv_na_val` values in `rivers`
#         river.loc[river[hillslope_id] == 1, col_list] = riv_na_val

        
#         # fix zero and negative values for fields that are not acceptable
#         river_sec = [riv_col_len, riv_col_len_dir, riv_col_slope]

#         # remove possible `None` values from `river_sec`
#         while None in river_sec:
#             river_sec.remove(None)

#         # creating `river_sec_df`
#         river_sec_df = river.loc[:, river_sec]
#         river_sec_df.where(~(river_sec_df <= 0), riv_na_val, inplace=True)
#         river_sec_df.where(~(river_sec_df.isnull()), riv_na_val, inplace=True)
#         river.loc[:, river_sec] = river_sec_df

#         # downcast dtype of necessary columns to integer, if possible
#         fcols = river[col_list].select_dtypes('float').columns
#         river[fcols] = river[fcols].apply(pd.to_numeric, downcast='integer')
    
    # return
    return river, cat


def intersect_topology(
    cat: gpd.GeoDataFrame,
    cat_cols: Dict[str, str],
    riv: gpd.GeoDataFrame,
    riv_cols: Dict[str, str],
    shapefile: gpd.GeoDataFrame = None,
    outlet_id: Optional[Union[int, Sequence[int]]] = None,
    outlet_val: Optional[int] = -9999,
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    '''Intesecting `cat` and `riv` with a given `shapefile`
    or a specified `outlet_id` of interest.

    Parameters
    ----------
    cat: geopandas.GeoDataFrame
        The catchment GeoDataFrame output of `prepare_cat` function
    cat_cols: dict
        A dictionary with 'id' key corresponding to the `cat` ID
        column
    riv: geopandas.GeoDataFrame
        The river GeoDataFrame output of `prepare_riv` function
    riv_cols: int
        A dictionary with 'id' and 'next_id' keys corresponding to
        the `riv` ID and next downstream ID column
    shapefile: geopandas.GeoDataFrame, optional [defaults to None]
        A GeoDataFrame specifying the limits of the river network
        topology of interest
    outlet_id: int, or sequence of ints, optional [defaults to None]
        the ID or list of IDs specifying the most downstream river segment
        to be specified as the outlet of the river network topology
        in addition to any non-contributing subbasins
    outlet_val: int, optional [defailts to -9999]
        An integer being set to the most downstream river segments

    Returns
    -------
    cat_clipped: geopandas.GeoDataFrame
        The clipped catchments object
    riv_clipped: geopandas.GeoDataFrame
        The clipped rivers object
    '''
    # necessary definitions
    cat_col_id = cat_cols.get('id')
    riv_col_id = riv_cols.get('id')
    riv_col_next_id = riv_cols.get('next_id')

    if shapefile is not None:
        # check `shapefile` dtype
        assert isinstance(shapefile, gpd.GeoDataFrame), "`shapefile` must be of type geopandas.GeoDataFrame"
        # intersection with `cat`
        upstream_ids = cat.overlay(shapefile,
                                   how='intersection',
                                   keep_geom_type=True)[cat_col_id]

    elif outlet_id is not None:
        # if it is a single ID, make a set out of it
        outlet_id = set(outlet_id)

        # check if outlet_id is included in the `riv` IDs
        assert riv[riv_col_id].isin(outlet_id).any(), "`outlet_id` must be chosen from segments included in `riv`"

        # find upstream segments
        upstream_ids = set()
        for element in outlet_id:
            upstream_ids.update(find_upstream(gdf=riv,
                                              target_id=element,
                                              main_id=riv_col_id,
                                              ds_main_id=riv_col_next_id))
    else:
        raise NotImplementedError("Either `shapefile` or `outlet_id` must be specified")

    # index `cat` and `riv` and return `cat_clipped` and `riv_clipped`
    cat_clipped = cat.loc[cat[cat_col_id].isin(upstream_ids), :].copy()
    riv_clipped = riv.loc[riv[riv_col_id].isin(upstream_ids), :].copy()

    # assign `outlet_val` to the segments without downstream segments
    riv_clipped.loc[~riv_clipped[riv_col_next_id].isin(riv_clipped[riv_col_id]), riv_col_next_id] = outlet_val
    
    # reset the index
    cat_clipped.reset_index(drop=True, inplace=True)
    riv_clipped.reset_index(drop=True, inplace=True)
    
    # return
    return cat_clipped, riv_clipped


def box_to_geometry(
    lon_min: float = None,
    lat_min: float = None,
    lon_max: float = None,
    lat_max: float = None,
    lims_box: list = [],
    crs: int = 4326,
    **kwargs,
) -> gpd.GeoDataFrame:
    '''Converts a latitude and longitude bounds to a geopandas
    GeoDataFrame object

    Parameters
    ----------
    lon_min: float, optional [defaults to None]
        Minimum longitude value
    lat_min: float, optional [defaults to None]
        Minimum latitude value
    lon_max: float, optional [defaults to None]
        Maximum longitude value
    lon_min: float, optional [defaults to None]
        Maximum latitude value
    box: list of floats, optional [defaults to [] ]
        List of spatial bounds in the ``[lon_min, lat_min, lon_max, lat_max]``
        form
    crs: int, optional [defaults to 4326]
        Coordinate Reference System EPSG numeber, defaults to
        WGS84 [EPSG:4326]

    Returns
    -------
    geom: geopandas.GeoDataFrame
        a GeoDataFrame object containing the geometry of the
        specified spatial bounds
    '''

    if lon_min and lat_min and lon_max and lat_max is not None:
        geometry = box(lon_min, lat_min, lon_max, lat_max, ccw=True)
    elif any(lims_box):
        geometry = box(*lims_box)
    else:
        raise ValueError("Either the spatial limits or a list of them must be provided")

    return gpd.GeoDataFrame(index=[0], crs=f'epsg:{crs}', geometry=[geometry])