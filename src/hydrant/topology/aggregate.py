"""
Aggregation methods for river and subbasin geometries
"""
from __future__ import annotations
import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx

import pandas as pd
import numpy as np
import networkx as nx
import hydrant.topology.geom as gm
import matplotlib.pyplot as plt
from typing import (
    Union,
    List
)


import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import hydrant.topology.geom as gm
import subprocess
import os
from   shapely.geometry import Point
import warnings

from typing import (
    Union,
    List,
    Dict,
    Callable,
)

ID_type = Union[int, str]

def union_seg(
    gdf: gpd.GeoDataFrame,
    by: str,
    id_dict: Dict[List[ID_type], ID_type],
    agg_dict: Dict[str, Callable] = None,
) -> gpd.GeoDataFrame:
    '''Return modified `gdf` after substituting `ids` with their union
    and applying `aggfunc`'s functions to their attributes

    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        GeoDataFrame containing various geometries with at least one attribute
    by: str
        Attribute name of the `gdf` for which `ids` are given
    id_dict: dict
        A dictionary with keys of list of IDs of elements the union
        of which is needed and values of the corresponding ID of the
        element for which the union is reported for
    agg_dict: Dict
        A dictionary consisting of keys corresponding to the attributes
        of `gdf` and values corresponding to a function to be applied
        to the attribute values of `ids`

    Returns
    -------
    modified_gdf: geopandas.GeoDataFrame
        Modified GeoDataFrame where the union of `ids` replaces the
        `ids` themselves.
    '''

    # getting list of columns from `gdf`
    cols = gdf.columns.tolist() # getting list of columns
    # removing column names that have specific aggregation
    # method(s) defined for them in the `agg_dict`
    for i in agg_dict.keys(): # removing those having specific aggregation methods
        cols.remove(i)
    # also removing `geometry` column from `gdf` since it is a
    # geopandas.GeoDataFrame object
    cols.remove('geometry')

    # iterate over id_dict elements, `id_list` is a tuple
    for id_list, target in id_dict.items():

        # extracting `ids` from `gdf`
        selected_gdf = gdf.loc[gdf[by].isin(id_list), :].copy()

        # union of `ids`
        union = selected_gdf.unary_union

        # implementing aggregations defined for each column in `agg_dict`;
        # first, preparing the `agg_gdf`
        agg_gdf = selected_gdf.drop(columns='geometry') # gpd.GeoDataFrame to pd.DataFrame
        agg_gdf.set_index(keys=[by], drop=True, inplace=True) # setting index to `by`

        # defining a lambda function which returns "target" (`t`) values
        # for a given "column" (`x`) if it is in the list of columns
        # defined above `cols` or `c` here; note that `c` contains column
        # names that have not been defnined in the `agg_dict`
        f = lambda x, t, c: x.loc[t] if x.name in c else agg_dict[x.name](x)

        # implementing aggregations
        agg_gdf = agg_gdf.agg(func=f, t=target, c=cols)

        # making a geopandas.GeoDataFrame out of created geometry in `union`
        # as it is needed for the concatenation step below
        agg_gdf_geom = gpd.GeoDataFrame(agg_gdf.to_frame().T, geometry=gpd.GeoSeries(union))
        # the `by` itself was set as an index, so the value is assigned here
        agg_gdf_geom[by] = target
        # set the .crs of `agg_gdf_geom` same as gdf
        agg_gdf_geom.set_crs(gdf.crs, inplace=True)

        # removing `id_list` from gdf (`rem_gdf`), concatenating `agg_gdf_geom`,
        # and returning `modified_gdf`
        gdf.drop(gdf.index[gdf[seg_str].isin(id_list)], inplace=True)
        gdf = gpd.GeoDataFrame(pd.concat([gdf, agg_gdf_geom]))

    return gdf



# Suppress SettingWithCopyWarning
warnings.filterwarnings('ignore')

def find_upstream(
    gdf: gpd.GeoDataFrame,
    target_id: Union[str, int, ...],
    main_id: str,
    ds_main_id: str,
) -> set[...]:
    '''Find "ancestors" or upstream segments in a river network given
    in the from of a geopandas.GeoDataFrame `gdf`
    
    Parameters
    ----------
    gdf: geopandas.GeoDataFrame
        GeoDataFrame of river segments including at least three pieces
        of information: 1) geometries of segments, 2) segment IDs, and
        3) downstream segment IDs
    target_id: str, int, or any other data type as included in `gdf`
        Indicating the target ID anscestor or upstream of which is
        desired
    main_id: str
        String defining the column of element IDs in the input geopandas
        dataframe
    ds_main_id: str
        String defining the column of downstream element IDs in the
        input geopandas dataframe
    
    Returns
    -------
    nodes: list
        IDs of nodes being upstream or anscestor of the `target_id`
    
    '''
    # creating a DiGraph out of `gdf` object
    riv_graph = nx.from_pandas_edgelist(gdf,
                                        source=main_id,
                                        target=ds_main_id,
                                        create_using=nx.DiGraph)

    # return nodes in a list
    nodes = nx.ancestors(riv_graph, target_id)

    # adding `target_id` as the last node of the branch
    nodes.add(target_id)

    return nodes

def main_branch(df, df_info):
    
    # get the values
    id_name     = df_info.get('id')
    next_name   = df_info.get('next')
    uparea_name = df_info.get('uparea')
    area_name   = df_info.get('area')
    
    # Create a directed graph
    G = nx.DiGraph()

    # Add edges and assign upa as edge weight
    for _, row in df.iterrows():
        G.add_edge(row[next_name], row[id_name], weight=row[uparea_name])
        
    # get the longest distance weighted based on up area
    longest_path = nx.dag_longest_path(G, weight='weight')

    #print("Longest distance:", longest_path)
    
    # Set flag to 1 where id is in the array_to_check
    df.loc[:,'main_branch'] = 0
    df.loc[df[id_name].isin(longest_path), 'main_branch'] = 1
    
    return df

def pfaf_one_round(df, df_info):
    
    # get the values
    id_name     = df_info.get('id')
    next_name   = df_info.get('next')
    uparea_name = df_info.get('uparea')
    area_name   = df_info.get('area')
    
    # add the pfaf_temp and get the main branch
    df.loc[:,'pfaf_temp'] = 1
    df = main_branch(df, df_info)
    
    # Separate DataFrame based on flag value
    df_main      = df[df['main_branch'] == 1].sort_values(by=uparea_name, ignore_index=True)
    
    if len(df) != len(df_main):
        
        # identify 4 largest upstream segments to main branch
        df_none_main = df[df['main_branch'] == 0].sort_values(by=uparea_name, ignore_index=True, ascending=False)
        max_4_up = df_none_main.loc[df_none_main[next_name].isin(df_main[id_name])].head(4)
        max_4_up = max_4_up.sort_values(by=uparea_name, ignore_index=True)

        # attach the uparea of the next down ID, which are on main ID, to the max_4_up dataframe
        max_4_up['next_up_area'] = 0
        max_4_up['up_confluence_main_id'] = 0
        for index, row in max_4_up.iterrows():
            # get the up area of df_main that is downstream of max_4_up
            max_4_up.loc[index,'next_up_area'] = df_main[uparea_name].loc[df_main[id_name]==row[next_name]].values
            index_temp = df_main.loc[df_main[id_name]==row[next_name]].index-1
            index_temp = np.array(index_temp).item()
            max_4_up.loc[index,'up_confluence_main_id'] = df_main[id_name].loc[index_temp] 

        max_4_up = max_4_up.sort_values(by='next_up_area', ignore_index=True)
        
        #print(max_4_up)
        
        if max_4_up.empty:
            raise ValueError("Error: max_4_up is empty")

        # get the len of max_4_up unique elements
        if len(np.unique(max_4_up[next_name])) == 4:
            odd_pfafs = [3,5,7,9]
        elif len(np.unique(max_4_up[next_name])) == 3:
            odd_pfafs = [3,5,9]
        elif len(np.unique(max_4_up[next_name])) == 2:
            odd_pfafs = [3,9]
        elif len(np.unique(max_4_up[next_name])) == 1:
            odd_pfafs = [9]

        # get the len of max_4_up elements
        if len(max_4_up[next_name]) == 4:
            even_pfafs = [2,4,6,8]
        elif len(max_4_up[next_name]) == 3:
            even_pfafs = [2,4,6]
        elif len(max_4_up[next_name]) == 2:
            even_pfafs = [2,4]
        elif len(max_4_up[next_name]) == 1:
            even_pfafs = [2]

        # add the segment and station into a data frame
        seg_ids = np.array([])
        pfaf_codes = np.array([])

        # Add values to the array for main branch
        seg_ids = np.append(seg_ids, max_4_up['up_confluence_main_id'])
        pfaf_codes = np.append(pfaf_codes, np.flip(odd_pfafs, axis=0))

        # add values to array for tributaries
        seg_ids = np.append(seg_ids, max_4_up[id_name])
        pfaf_codes = np.append(pfaf_codes, np.flip(even_pfafs, axis=0))
        zipped = zip(seg_ids, pfaf_codes)
        sorted_zipped = sorted(zipped, key=lambda x: x[1])

        # get the upstream and assign the pfaf from smaller values to largest values
        df.loc[:,'pfaf_temp'] = 1

        for seg_id, pfaf_code in sorted_zipped:
            ids_selected = find_upstream(df, seg_id, id_name, next_name)
            # replace the ids in df
            indices = df[df[id_name].isin(ids_selected)].index
            # Update 'flag' column to 1 for the rows with matching indices
            df.loc[indices, 'pfaf_temp'] = pfaf_code

    return df


def find_separated_graphs(df, df_info):
    
    # get the values
    id_name     = df_info.get('id')
    next_name   = df_info.get('next')
    
    # Create a graph from the dataframe
    G = nx.from_pandas_edgelist(df, id_name, next_name)

    # Find connected components in the graph
    connected_components = list(nx.connected_components(G))

    # Map each node to its connected component
    node_to_component = {}
    for i, component in enumerate(connected_components, 1):
        for node in component:
            node_to_component[node] = i

    # Create a new column in the dataframe indicating the component number
    df['graph_number'] = df[id_name].map(node_to_component)

    return df


def transform_array(arr):
    # Find the maximum value in the array
    max_value = np.max(arr)
    
    # Determine the number of digits in the maximum value
    num_digits = len(str(max_value))
    
    # Calculate the multiplier
    multiplier = 10 ** num_digits
    
    # Transform the array
    transformed_arr = arr * multiplier
    
    return transformed_arr


def pfaf_sub(df, df_info = {'id':'COMID', 'next':'NextDownID', 'uparea': 'uparea', 'area': 'unitarea'}, depth = 10):
    
    # get the values
    id_name     = df_info.get('id')
    next_name   = df_info.get('next')
    uparea_name = df_info.get('uparea')
    area_name   = df_info.get('area')
    
    # get the current order based on id to sort later
    order_id = df[id_name].values
    
    # initial pfaf set up for the df
    df = pfaf_one_round (df, df_info=df_info)
    df['pfaf'] = df['pfaf_temp']
    
    for i in np.arange(2, depth):
        if len(df) != len(np.unique(df ['pfaf'])):
            df_slice_total = pd.DataFrame()
            for m in np.unique(df ['pfaf']):
                df_slice = df[df['pfaf']==m]
                df_slice = pfaf_one_round(df_slice, df_info)
                df_slice_total = pd.concat([df_slice_total,df_slice], ignore_index=True)
            df = df_slice_total.copy()
            #print(df)
        else:
            df['pfaf_temp'] = 0

        # Convert values in column1 and column2 to strings
        df['pfaf'] = df['pfaf'].astype(str)
        df['pfaf_temp'] = df['pfaf_temp'].astype(str)

        # Concatenate the strings in column1 and column2
        df['pfaf'] = df['pfaf'] + df['pfaf_temp']

        # Convert the concatenated string back to integers
        #df['pfaf'] = df['pfaf'].astype(int)
        df.loc[:,'pfaf'] = df.loc[:,'pfaf'].astype(int)
        
        # set back to the original order
        df = df.set_index(id_name).loc[order_id].reset_index()
        
    return df



def pfaf(df, df_info = {'id':'COMID', 'next':'NextDownID', 'uparea': 'uparea', 'area': 'unitarea'}, depth = 10):
    
    # get the values
    id_name     = df_info.get('id')
    next_name   = df_info.get('next')
    uparea_name = df_info.get('uparea')
    area_name   = df_info.get('area')
    
    # get the current order based on id to sort later
    order_id = df[id_name].values
    
    # replace the possible 0 and negative values for nextid for networkx
    arr = np.array(df[next_name].values)
    max_value = arr.max() + 1
    for i in range(len(arr)):
        if arr[i] <= 0:
            arr[i] = max_value
            max_value += 1
    df[next_name] = arr
    
    # get separated regions
    df = find_separated_graphs(df, df_info)
    df ['pfaf_region'] = transform_array(df ['graph_number'].values)
    df = df.drop(columns= ['graph_number'])
    print(np.unique(df['pfaf_region'].values))
    
    # loop over separated regions
    df_all = pd.DataFrame()
    for i in np.unique(df['pfaf_region'].values):
        # slice the shapefile
        df_slice = df[df['pfaf_region']==i]
        df_slice = pfaf_sub(df_slice, df_info, depth)
        df_all = pd.concat([df_slice,df_all], ignore_index=True)
        
    # get the df in order
    df = df_all.set_index(id_name).loc[order_id].reset_index()    
    
    # drop the extra columns
    df = df.drop(columns=['pfaf_temp','main_branch'])
        
    return df

def pfaf_agg (df, df_info={'pfaf': 'pfaf', 'pfaf_region': 'pfaf_region'}, depth_agg = 2):
    
    # get the values
    pfaf        = df_info.get('pfaf')
    pfaf_region = df_info.get('pfaf_region')
    
    #
    def extract_first_n_digits(x,n=depth_agg):
        return str(x)[:n]
    
    # Apply the function to create a new column
    df['pfaf_temp'] = df[pfaf].apply(extract_first_n_digits)

    # Dissolve based on the new column
    df_agg = df.dissolve(by=['pfaf_temp',pfaf_region])
    
    # return
    return df_agg