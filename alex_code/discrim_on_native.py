# TODO: drop outlier graphs
# TODO: fix off-by-one error
""" 
Note: For CPAC, because datasets were run with different ANT settings, 
"""
#%%
import os
import re
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import networkx
from sklearn.utils import check_X_y
from sklearn.metrics import euclidean_distances

import graspy
from graspy.utils import pass_to_ranks


def path_files_list(path: str):
    """
    From a path containing ssv edgelist files, 
    return a list of adjacency matrices, 
    sorted alphabetically.
    
    Parameters
    ----------
    path_list : list of ssv files.
    
    Returns
    -------
    List of paths.
    """
    path = Path(path).absolute()
    return sorted([str(fl) for fl in path.iterdir()])


def get_X(arrays: List, PTR=True):
    """
    Ravel every array in list_of_arrays, 
    return single array with rows corresponding to each raveled array.

    Parameters
    ----------
    list_of_arrays : list of arrays, each the same size.

    Returns
    -------
    X : Concatenated X matrix.
    """
    # TODO: check to make sure all of the arrays are the same length
    # TODO: make sure nodes aren't being sorted differently for each array

    # pass to ranks first
    if PTR:
        arrays = [pass_to_ranks(array) for array in arrays]

    # stack 'em up
    return np.stack(list(map(np.ravel, arrays)), axis=0)


# folder = "/Users/alex/Dropbox/NeuroData/ndmg-paper/data/graphs/native_graphs_HNU1"
# path_list = path_files_list(folder)
# arrays = graspy.utils.import_edgelist(path_list, extension="ssv", delimiter=" ")

#%%
def get_Y(list_of_ndmg_paths):
    # TODO: get target vector from ndmg graphs directory
    """
    From a list of paths to graph files, return a list of the subjects those graphs belong to.

    Parameters
    ----------
    list_of_ndmg_paths : list of paths to individual graphs.

    Returns
    -------
    Y : 1d array, shape (n_samples), type(str)
    """
    subrgx = re.compile(r"(sub-)([0-9]*)(_)")

    # used to be a list comprehension, but this is more clear
    subjects_list = []
    for path in list_of_ndmg_paths:
        # to make sure there aren't any random files
        if Path(path).suffix == ".ssv":
            subject_re = re.search(subrgx, path)
            subject = subject_re.group(2)
            subjects_list.append(subject)

    return np.array(subjects_list)


def discr_stat(
    X, Y, dissimilarity="euclidean", remove_isolates=True, return_rdfs=False
):
    """
    Computes the discriminability statistic.

    Parameters
    ----------
    X : array, shape (n_samples, n_features) or (n_samples, n_samples)
        Input data. If dissimilarity=='precomputed', the input should be the dissimilarity matrix.

    Y : 1d-array, shape (n_samples)
        Input labels.

    dissimilarity : str, {"euclidean" (default), "precomputed"}
        Dissimilarity measure to use:

        - 'euclidean':
            Pairwise Euclidean distances between points in the dataset.

        - 'precomputed':
            Pre-computed dissimilarities.

    remove_isolates : bool, optional, default=True
        Whether to remove data that have single label.

    return_rdfs : bool, optional, default=False
        Whether to return rdf for all data points.

    Returns
    -------
    stat : float
        Discriminability statistic. 

    rdfs : array, shape (n_samples, max{len(id)})
        Rdfs for each sample. Only returned if ``return_rdfs==True``.
    """
    check_X_y(X, Y, accept_sparse=True)

    uniques, counts = np.unique(Y, return_counts=True)
    if (counts != 1).sum() <= 1:
        msg = "You have passed a vector containing only a single unique sample id."
        raise ValueError(msg)
    if remove_isolates:
        idx = np.isin(Y, uniques[counts != 1])
        labels = Y[idx]

        if dissimilarity == "euclidean":
            X = X[idx]
        else:
            X = X[np.ix_(idx, idx)]
    else:
        labels = Y

    if dissimilarity == "euclidean":
        dissimilarities = euclidean_distances(X)
    else:
        dissimilarities = X

    rdfs = _discr_rdf(dissimilarities, labels)
    stat = np.nanmean(rdfs)

    if return_rdfs:
        return stat, rdfs
    else:
        return stat


def _discr_rdf(dissimilarities, labels):
    """
    A function for computing the reliability density function of a dataset.

    Parameters
    ----------
    dissimilarities : array, shape (n_samples, n_features) or (n_samples, n_samples)
        Input data. If dissimilarity=='precomputed', the input should be the 
        dissimilarity matrix.

    labels : 1d-array, shape (n_samples)
        Input labels.

    Returns
    -------
    out : array, shape (n_samples, max{len(id)})
        Rdfs for each sample. Only returned if ``return_rdfs==True``.
    """
    check_X_y(dissimilarities, labels, accept_sparse=True)

    rdfs = []
    for i, label in enumerate(labels):
        di = dissimilarities[i]

        # All other samples except its own label
        idx = labels == label
        Dij = di[~idx]

        # All samples except itself
        idx[i] = False
        Dii = di[idx]

        rdf = [1 - ((Dij < d).sum() + 0.5 * (Dij == d).sum()) / Dij.size for d in Dii]
        rdfs.append(rdf)

    out = np.full((len(rdfs), max(map(len, rdfs))), np.nan)
    for i, rdf in enumerate(rdfs):
        out[i, : len(rdf)] = rdf

    return out


# def main(folder):
#     path_list = path_files_list(folder)
#     Y = get_Y(path_list)
#     arrays_list = graspy.utils.import_edgelist(
#         path_list, extension="ssv", delimiter=" "
#     )
#     X = get_X(arrays_list)
#     return discr_stat(X, Y)


#%%
def main(flist):
    for i, val in enumerate(flist):
        print("Dataset: {}".format(Path(flist[i]).name))
        path_list = path_files_list(val)
        Y = get_Y(path_list)
        arrays_list = graspy.utils.import_edgelist(
            path_list, extension="ssv", delimiter=" "
        )
        X = get_X(arrays_list)
        print("Number of observations from target vector: {}".format(len(Y)))
        print("Number of observations from X: {}".format(X.shape[0]))
        print("Unique number of observations from Y: {}".format(len(np.unique(Y))))
        print("Unique subjects considered:", np.unique(Y))
        print("\n")
        print("First (5,5) subset of X: \n {}".format(X[:5, :5]))
        print("Last column of X: {}".format(X[:, -1]))
        print("Discriminability: {}".format(discr_stat(X, Y)))
        print("\n")
        print("-----------------------------------")
        print("\n")


#%%
# testing
folders = Path("/Users/alex/Dropbox/NeuroData/ndmg-paper/data/graphs")
flist = [str(folder) for folder in folders.iterdir()]
val = flist[0]
path_list = path_files_list(val)
Y = get_Y(path_list)
arrays_list = graspy.utils.import_edgelist(path_list, extension="ssv", delimiter=" ")
X = get_X(arrays_list)
# main(flist)
