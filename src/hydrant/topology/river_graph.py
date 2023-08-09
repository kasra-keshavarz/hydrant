"""
Tools for river network analysis in Python using efficient algorithms
"""

from __future__ import annotations

import networkx as nx
import geopandas as gpd

from typing import (
    Union,
)

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
    riv_graph = nx.from_pandas_edgelist(gdf, source=main_id, target=ds_main_id, create_using=nx.DiGraph)
    
    # return nodes in a list
    nodes = nx.ancestors(riv_graph, target_id)
    
    return nodes


