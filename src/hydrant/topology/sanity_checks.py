"""
Common sanity checks on river and subbasin geometries
"""

import geopandas as gpd
import numpy as np
import pandas as np

from typing import ( 
    List,
    Dict,
)

ID_type = Union[int, str, ...]

def spatial_conn(
    gdf: gpd.GeoDataFrame,
    main_id: str=None,
    ds_main_id: str=None,
    print_err: bool=True,
) -> Tuple[Dict, Dict]:
    
    '''Validates downstream IDs and checks spatial connectivity 
    of river elements
    
    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        a GeoPandas dataframe containing details of a river network
        with minimum details of element ID, downstream element ID, and
        geometries
    main_id: str, defaults to `None`
        String defining the column of element IDs in the input geopandas
        dataframe
    ds_main_id: str, defaults to `None`
        String defining the column of downstream element IDs in the
        input geopandas dataframe
    print_err: bool, defaults to `True` (optional)
        Argument to specify either print the detailed errors within the
        input `gdf`
        
    Returns
    -------
    connections: Dict
        A dictionary of connections with keys as `main_id` and
        values of `ds_main_id` if no issues were found
    wrong_conns: Dict
        A dictionary of problematic connections with keys as `main_id`
        and values defining the nature of the problem, available error
        types: 1) spatial disconnection (`spdis`), 2) missing downstream
        element (misds)
        
    '''
    
    # necessary initializations
    connections = {} # key: `main_id`, value: `ds_main_id`
    wrong_conns = {} # key: `main_id`, value: error type
    
    
    if main_id is None:
        raise ValueError("`main_id` cannot be `None`")
        
    if ds_main_id is None:
        raise ValueError("`ds_main_id` cannot be `None`")

    # Loop through each feature in the shapefile
    for idx, feature in gdf.iterrows():

        # Skip over features with DS_Main_ID = 0
        if feature[ds_main_id] == 0:
            continue

        # Find the corresponding descendant feature
        descendant = gdf.loc[gdf[main_id] == feature[ds_main_id]]

        # Check if the descendant DataFrame is empty
        if descendant.empty:
            # find possible downstream IDs of 'DS_Main_ID' does not exist
            possible_ds = gdf.geometry.intersects(gdf.loc[idx].geometry)
            possible_ds_list = gdf[possible_ds][main_id].tolist()
            
            # add the `Main_ID` to the `wrong_conns` list
            wrong_conns[feature[main_id]] = 'misds'

            # print details
            if print_err:
                print(f"Warning: No descendant feature found for `main_id`={feature[main_id]} with `ds_main_id`={feature[ds_main_id]}")
                print(f"  Possible downstream IDs are: {possible_ds_list}")

            continue

        # Check if the feature's geometry intersects with the descendant's geometry
        if not feature.geometry.intersects(descendant.geometry.iloc[0]):

            if print_err:
                print(f"Feature with main_id {feature[main_id]} is not spatially connected to its descendant with ds_main_id {feature[ds_main_id]}")
                
            # add the `Main_ID` to the `wrong_conns` dictionary
            wrong_conns[feature[main_id]] = 'spdis'
        else:
            # add the `Main_ID` to the `connections` dictionary
            connections[feature[main_id]] = feature[ds_main_id]
    
    return connections, wrong_conns

def find_cycles(
    gdf: gpd.GeoDataFrame,
    main_id: str,
    ds_main_id: str,
    l: int=None,
    print_err: bool=True,
) -> list[...]:
    '''Returns the cycles of length `l` of a river network given 
    as a geopandas GeoDataFrame
    
    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        GeoDataFrame of river segments including at least three pieces
        of information: 1) geometries of segments, 2) segment IDs, and
        3) downstream segment IDs
    main_id: str, defaults to `None`
        String defining the column of element IDs in the input geopandas
        dataframe
    ds_main_id: str, defaults to `None`
        String defining the column of downstream element IDs in the
        input geopandas dataframe
    l: int, defaults to None (optional)
        length of cycles to be found, leave as None to find
        cycles of any length
    print_err: bool, defaults to `True` (optional)
        print cycles if stated
    
    Returns
    -------
    cycles: list
        List of `main_id`s indicating cycles
    
    '''
    
    # creating a DiGraph out of `gdf` object
    riv_graph = nx.from_pandas_edgelist(gdf, source=main_id, target=ds_main_id, create_using=nx.DiGraph)

    # extracting cycles
    cycles = list(nx.simple_cycles(riv_graph, l))
    
    if print_err:
        print(cycles)
    
    return cycles

def sanitize_connectivity(
    gdf: gpd.GeoDataFrame,
    main_id: str,
    ds_main_id: str,
    id_dict: Dict[List[ID_type], ID_type],
) -> gpd.GeoDataFrame:
    '''Corrects river network connectivity after union aggregations
    implemented on a few segments. Modified GeoDataFrame where the
    segments with end to elements in `id_list` that are not `target`
    have now a `ds_main_id` of `target`
    
    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        GeoDataFrame of river segments including at least three pieces
        of information: 1) geometries of segments, 2) segment IDs, and
        3) downstream segment IDs
    main_id: str
        String defining the column of element IDs in the input geopandas
        dataframe
    ds_main_id: str
        String defining the column of downstream element IDs in the
        input geopandas dataframe
    id_dict: dict
        A dictionary with keys of list of IDs of elements the union
        of which was added and values of the corresponding ID of the
        element for which the union was reported for  
    
    Returns
    -------
    modified_gdf: geopandas.GeoDataFrame
        Modified GeoDataFrame where the segments with end to elements in `id_list`
        that are not `target` have now a `ds_main_id` of `target`
    '''
    
    modified_gdf = gdf
    
    # loop over elements of `id_dict`
    for id_list, target in id_dict.items():
    
        # correcting ids in the `id_list`
        for segment in id_list:
            # skip if it is the target itself
            if segment == target:
                continue
            # change the downstream value to `target
            else:
                modified_gdf.loc[modified_gdf[ds_main_id] == segment, ds_main_id] = target
    # reseting index
    modified_gdf.reset_index(drop=True, inplace=True)
    
    return modified_gdf

