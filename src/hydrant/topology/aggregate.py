"""
Aggregation methods for river and subbasin geometries
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import networkx as nx

from typing import (
    Union,
    List,
    Dict,
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

