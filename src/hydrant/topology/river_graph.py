"""
Tools for river network analysis in Python using efficient algorithms
"""

from __future__ import annotations

import networkx as nx
import geopandas as gpd

from typing import (
    Union,
    List
)


def _return_nx_graph(
    gdf: gpd.GeoDataFrame,
    main_id: str,
    ds_main_id: str,
) -> nx.DiGraph:
    '''Return networkx DiGraph of the input river network
    '''
    return nx.from_pandas_edgelist(gdf,
                                   source=main_id,
                                   target=ds_main_id,
                                   create_using=nx.DiGraph)


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
    riv_graph = _return_nx_graph(gdf,
                                 main_id,
                                 ds_main_id)

    # return nodes in a list
    nodes = nx.ancestors(riv_graph, target_id)

    # adding `target_id` as the last node of the branch
    nodes.add(target_id)

    return nodes


def find_downstream(
    gdf: gpd.GeoDataFrame,
    target_id: Union[str, int, ...],
    main_id: str,
    ds_main_id: str,
) -> set[...]:
    '''
    Returns "descendants" or downstream segments in a river network given
    in the from of a geopandas.GeoDataFrame (`gdf`)

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
    riv_graph = _return_nx_graph(gdf,
                                 main_id,
                                 ds_main_id)

    # return nodes in a list
    nodes = nx.descendants(riv_graph, target_id)

    return nodes


def longest_branch(
    riv: gpd.GeoDataFrame,
    main_id: Union[str, int] = None,
    ds_main_id: Union[str, int] = None,
) -> List[str, int]:
    """Returns nodes of the longest branch of a river network

    Parameters
    ----------
    riv : gpd.GeoDataFrame
        a geopandas.GeoDataFrame object containing the geometry,
        `main_id`, and `ds_main_id` of a river network
    main_id : str or int, defaults to `None`
        column label within `riv` corresponding to ID values of river
        segments
    ds_main_id : str or int, defaults to `None`
        column label within `riv` corresponding to downstream segments of
        each river

    Returns
    -------
    nodes : set
        a set object containing river segment `main_id` calues of the
        longest branch found in `riv`
    """
    riv_graph = nx.from_pandas_edgelist(riv, main_id, ds_main_id)

    nodes = nx.dag_longest_path(riv_graph)

    return nodes
