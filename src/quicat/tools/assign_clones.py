from typing import Dict, List, Union

import anndata
import numpy as np


def assign_clones(
    adata: anndata.AnnData,
    clone_mapping: Dict[str, Union[str, List[str]]],
    threshold: float = 0.5,
    key_added: str = "clone",
) -> anndata.AnnData:
    """
    Assigns clones to cells in the AnnData object based on barcode match ratio.

    For each cell, the function computes the ratio of belonging to each clone
    by calculating the fraction of barcodes matching between the cell and the clone.
    It assigns the clone with the highest ratio to the cell if the ratio
    exceeds the specified threshold.

    Parameters:
        adata (anndata.AnnData): The AnnData object containing single-cell data.
        clone_mapping (Dict[str, Union[str, List[str]]]): A dictionary mapping clone names to barcode(s).
            Barcodes can be a single string or a list of strings.
        threshold (float): The minimum match ratio required to assign a cell to a clone.
            Should be between 0.0 and 1.0. Default is 0.5.
        key_added (str): The key to use for the 'obs' column containing the assigned clones.
    Returns:
        anndata.AnnData: The modified AnnData object with the 'clone' column added to 'obs'.
    """
    clone_barcodes_sets: Dict[str, set] = {}
    for clone_name, clone_barcodes in clone_mapping.items():
        if isinstance(clone_barcodes, str):
            clone_barcodes_sets[clone_name] = {clone_barcodes}
        elif isinstance(clone_barcodes, list):
            clone_barcodes_sets[clone_name] = set(clone_barcodes)
        else:
            e = f"Barcodes for clone '{clone_name}' must be a string or list of strings."
            raise TypeError(e)

    # Now 'barcodes' is a set of all unique barcodes
    barcodes = set.union(*clone_barcodes_sets.values())
    clones: List[str] = []

    for cell in adata.obs_names:
        cell_expression = adata[cell].X
        if hasattr(cell_expression, "toarray"):
            cell_expression = cell_expression.toarray().flatten()
        else:
            cell_expression = np.array(cell_expression).flatten()
        cell_barcodes = adata.var_names[cell_expression > 0].values
        cell_barcodes = set(cell_barcodes) & barcodes

        highest_overlap: float = 0.0
        assigned_clone: str = "Unassigned"

        for clone_name, clone_barcodes_set in clone_barcodes_sets.items():
            intersection = cell_barcodes & clone_barcodes_set
            match_overlap = len(intersection) / len(clone_barcodes_set)
            if (match_overlap > highest_overlap) and (len(cell_barcodes) <= len(clone_barcodes_set)):
                highest_overlap = match_overlap
                assigned_clone_candidate = clone_name

        if highest_overlap >= threshold and highest_overlap > 0.0:
            assigned_clone = assigned_clone_candidate

        clones.append(assigned_clone)

    adata.obs[key_added] = clones

    return adata
