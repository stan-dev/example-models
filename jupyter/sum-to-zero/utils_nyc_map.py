import numpy as np
import pandas as pd
import geopandas as gpd
from libpysal.weights import Queen, W
from typing import (Dict, List)

def disconnect_nbs(
        nbs: W, target_indices: List[int], remove_indices: List[int]
        ) -> None:
    """
    Modify a neighbors graph by removing indices from the list of neighbors
    and resetting neighbors, weights lists accordingly.
    - param nbs: neighbor graph
    - param target_indices: list of node ids to process 
    - param remove_indices: list of non-neighbor node ids
    """
    remove_set = set(remove_indices)  # Convert to set for O(1) lookups
    for node in target_indices:
        # Filter neighbors and weights in a single pass
        clean_neighbors = []
        clean_weights = []
        for neighbor, weight in zip(nbs.neighbors[node], nbs.weights[node]):
            if neighbor not in remove_set:
                clean_neighbors.append(neighbor)
                clean_weights.append(weight)
        nbs.neighbors[node] = clean_neighbors
        nbs.weights[node] = clean_weights

def nyc_cleanup(nbs: W, gdf: gpd.GeoDataFrame) -> W:
    """
    Modify neighbor graph of NYC to remove neighbor pairs between
    Manhattan and other boroughs (Brooklyn, Queens).
    - param nbs: neighbor graph
    - param gdf : geopandas.GeoDataFrame of NYC Census Tract data    
    """
    # Get indices for each borough
    manhattan_indices = gdf[gdf['BoroName'] == 'Manhattan'].index
    brooklyn_indices = gdf[gdf['BoroName'] == 'Brooklyn'].index
    queens_indices = gdf[gdf['BoroName'] == 'Queens'].index
    brooklyn_and_queens = list(brooklyn_indices) + list(queens_indices)
    # Disconnect Manhattan from Brooklyn/Queens and vice versa
    disconnect_nbs(nbs, manhattan_indices, brooklyn_and_queens)
    disconnect_nbs(nbs, brooklyn_indices, manhattan_indices)
    disconnect_nbs(nbs, queens_indices, manhattan_indices)

    return W(nbs.neighbors, nbs.weights)


def nyc_sort_by_comp_size(nyc_gdf: gpd.GeoDataFrame
                              ) -> tuple[W, gpd.GeoDataFrame, List[int]]:
    """
    Process NYC geodataframe - sort by component size descending
    - param nyc_gdf : geopandas.GeoDataFrame
    - return: tuple containing: neighbors, sorted gdf,
              nodes_per_component, edges_per_component
    """
    # Compute initial neighborhood graph
    nyc_nbs = Queen.from_dataframe(nyc_gdf, geom_col='geometry')
    # Clean borough connections and get components
    nyc_nbs_tmp  = nyc_cleanup(nyc_nbs, nyc_gdf)
    # Add component info to dataframe
    nyc_gdf['comp_id'] = nyc_nbs_tmp.component_labels
    sizes = nyc_gdf['comp_id'].value_counts()
    nyc_gdf['comp_size'] = nyc_gdf['comp_id'].map(sizes)
    # Sort by component size and reset index
    nyc_gdf_sorted = (nyc_gdf.sort_values(by='comp_size', ascending=False)
                             .reset_index(drop=True)
                             .set_geometry('geometry'))
    # Recompute neighborhood graph, update gdf, get component sizes
    nyc_nbs_sorted = Queen.from_dataframe(nyc_gdf_sorted, geom_col='geometry')
    nyc_nbs_clean = nyc_cleanup(nyc_nbs_sorted, nyc_gdf_sorted)
    nyc_gdf_sorted['comp_id'] = nyc_nbs_clean.component_labels
    component_sizes = list(sizes.sort_values(ascending=False))

    return nyc_nbs_clean, nyc_gdf_sorted, component_sizes
