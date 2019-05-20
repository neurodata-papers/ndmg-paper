#%%
r"""Get all cpac graphs for (HNU1, BNU1, SWU4, KKI2009, NKI1) | graph in cpac,
for every cpac graph, concatenate it with its corresponding native-space graph,
put them all in an output folder to compute discriminability.
"""
# TODO: vectorize the disparate graphs into a single matrix,
#       and get a Y-vector as well.
#%%
import os
from pathlib import Path
from collections import OrderedDict
from collections import Iterable
from functools import partial
from functools import reduce
import warnings

import pandas as pd
import numpy as np
import networkx as nx
import graspy

# my own package with functions for making working with graphs easier
from graphutils.graph_io import get_X
from graphutils.graph_io import get_Y

if not Path(os.getcwd()).name == "alex_code":
    os.chdir("alex_code")

#%%
# steal graspy's `import_edgelist` function and hack it a bit
# so that it's better suited for the situation I'm in
def import_edgelist(
    path,
    extension="edgelist",
    delimiter=None,
    nodetype=int,
    return_vertices=False,
    output_dict=False,
):
    """
    Function for reading a single or multiple edgelists. When importing multiple 
    edgelists, the union of vertices from all graphs is computed so that each output
    graph have matched vertex set. The order of nodes are sorted by node values.
    Parameters
    ----------
    path : str, Path object, or iterable
        If ``path`` is a directory, then the importing order will be sorted in 
        alphabetical order.
    extension : str, optional
        If ``path`` is a directory, then the function will convert all files
        with matching extension. 
    delimiter : str or None, default=None, optional
        Delimiter of edgelist. If None, the delimiter is whitespace.
    nodetype : int (default), float, str, Python type, optional
       Convert node data from strings to specified type.
    return_vertices : bool, default=False, optional
        Returns the union of all ind
    output_dict : bool, default=False, optional
        If `path` is a file, does nothing.
        If True, returns dict of {'filename': np.array} key/value pairs.
    Returns
    -------
    out : list of array-like, or array-like, shape (n_vertices, n_vertices) or dict.
        If ``path`` is a directory, a list of arrays is returned. If ``path`` is both a directory and output_dict is True, a dictionary of {'filename': np.array} is returned. If ``path`` is a file, an array is returned.
    vertices : array-like, shape (n_vertices, )
        If ``return_vertices`` == True, then returns an array of all vertices that were 
        included in the output graphs. 
    """
    # p = Path(path)
    if not isinstance(path, (str, Path, Iterable)):
        msg = "path must be a string or Iterable, not {}".format(type(path))
        raise TypeError(msg)

    # get a list of files to import
    if isinstance(path, (str, Path)):
        p = Path(path)
        if p.is_dir():
            files = sorted(p.glob("*" + extension))
        elif p.is_file():
            files = [p]
        else:
            raise ValueError("No graphs founds to import.")
    else:  # path is an iterable
        files = [Path(f) for f in path]

    if len(files) == 0:
        msg = "No files found with '{}' extension found.".format(extension)
        raise ValueError(msg)

    graphs = [
        nx.read_weighted_edgelist(f, nodetype=nodetype, delimiter=delimiter)
        for f in files
    ]

    if all(len(G.nodes) == 0 for G in graphs):
        msg = (
            "All graphs have 0 vertices. Please double check if proper "
            + "'delimiter' is given."
        )
        warnings.warn(msg, UserWarning)

    # Compute union of all vertices
    vertices = np.sort(reduce(np.union1d, [G.nodes for G in graphs]))
    out = [nx.to_numpy_array(G, nodelist=vertices, dtype=np.float) for G in graphs]

    # only return adjacency matrix if input is only 1 graph
    if len(out) == 1:
        out = out[0]

    if return_vertices:
        return out, vertices
    elif output_dict and p.is_dir():
        filenames = [filepath.stem for filepath in files]
        dict_out = dict(zip(filenames, out))
        return dict_out

    else:
        return out


def concatenate_cpac_and_native(cpac: "np.array", native: "np.array"):
    """
    For each subject,
    and each session within each subject,
    concatenates the graph for cpac and native space.
    """
    # TODO: might need to do some checks on array size, etc
    if cpac.shape[0] == native.shape[0]:
        return np.hstack((cpac, native))
    else:
        print(
            f"current cpac graph is shape {cpac.shape} and current native graph is shape {native.shape}. Skipping."
        )
        pass


def final_concat(native_graphs, cpac_graphs, dataset):
    # Currently incredibly non-general,
    # and specific to the exact name formatting for
    # the particular set of native-space and cpac graphs I have.
    # TODO: is there a way to turn this into a general-form function?
    """Concatenate two dicts of graphs into one dict of graphs.
    
    for every dataset in cpac_graphs and native_graphs,
    for every array in that dataset,
    check if the subject is the same,
    and the session is the same.
    if both are, concatenate the graphs.

    native_graphs, cpac_graphs : dict, {filename: nxd array}
    dataset : str, name of dataset

    Returns : dict, {filename: nxp array}
    """
    concatenated_graphs = []
    for fgraph in native_graphs:
        # turn graph into the equivalent thing in cpac_graphs_KKI
        parts = fgraph.split("_")[0:2]  # e.g., ['sub-849', 'ses-2']
        pieces = [i.split("-")[1] for i in parts]  # e.g., [849, 2]
        sub, ses = pieces
        corresponding_cpac = f"{dataset}_{sub}_session_{ses}_correlation"

        # grab the cpac graph and native graph
        native_graph = native_graphs[fgraph]
        try:
            cpac_graph = cpac_graphs[corresponding_cpac]
        except KeyError:
            # TODO: make this more elegant
            print(
                f"parts: {parts}, pieces: {pieces}, corresponding_cpac: {corresponding_cpac}"
            )
            continue

        # gonna need to turn this into a dictionary appending thing
        graph = concatenate_cpac_and_native(cpac_graph, native_graph)
        key = f"sub_{sub}_session_{ses}"
        concatenated_graphs.append((key, graph))
    return concatenated_graphs


def return_X_and_Y(input_dataset):
    arrays = list(OrderedDict(input_dataset).values())
    names = list(OrderedDict(input_dataset).keys())
    names = [name.split("_")[1] for name in names]

    X = get_X(arrays, PTR=False)
    Y = np.array(names)
    return (X, Y)


def return_sorted_graph(graph_file: str, n_nodes: int):
    """From a graph file with 70 nodes, return a sorted adjacency matrix.

    graph_file : str, filepath to edgelist.
    n_nodes : int, number of nodes the numpy array should have.
    """
    graph = nx.read_weighted_edgelist(graph_file, nodetype=int, delimiter=" ")
    vertices = np.arange(1, n_nodes + 1)
    out = nx.to_numpy_array(graph, nodelist=vertices, dtype=np.float)
    return out


#%%
# Below: Make sure that the functions defined above
# are sorting graph nodes properly.
test_path = "/Users/alex/Dropbox/NeuroData/ndmg-paper/data/graphs/native_graphs_HNU1/sub-0025427_ses-1_dwi_desikan_space-MNI152NLin6_res-2x2x2_measure-spatial-ds_adj.ssv"
dataset = "/Users/alex/Dropbox/NeuroData/ndmg-paper/data/graphs/native_graphs_HNU1"
import_graph = partial(
    import_edgelist, extension="ssv", delimiter=" ", output_dict=True
)
dataset_dict = import_edgelist(
    dataset,
    extension="ssv",
    delimiter=" ",
    nodetype=int,
    return_vertices=False,
    output_dict=True,
)
np.all(
    dataset_dict[
        "sub-0025427_ses-1_dwi_desikan_space-MNI152NLin6_res-2x2x2_measure-spatial-ds_adj"
    ]
    == return_sorted_graph(test_path, 70)
)  # return True

#%%
# grab separate graphs
p = Path.cwd()
p = p.parent / "data" / "graphs"
fcpac_graphs = [folder for folder in p.iterdir() if "cpac" in folder.name]
fnative_graphs = [folder for folder in p.iterdir() if "native" in folder.name]
#%%
# {dataset.stem: {subj_ses: np.array}}
# everything is now correctly the right shape
cpac_graphs = {dataset.stem: import_graph(dataset) for dataset in fcpac_graphs}
native_graphs = {dataset.stem: import_graph(dataset) for dataset in fnative_graphs}
#%%

# how do I concatenate, now that I have the data separately organized nicely?
# I want to:
# for each dataset in cpac_graphs and native_graphs,
# pair up the graphs that are concurrent,
# and output {'sub-#_ses-#': graph}

cpac_graphs_KKI = cpac_graphs["cpac_graphs_KKI"]
cpac_graphs_HNU1 = cpac_graphs["cpac_graphs_HNU1"]
cpac_graphs_BNU1 = cpac_graphs["cpac_graphs_BNU1"]
cpac_graphs_SWU4 = cpac_graphs["cpac_graphs_SWU4"]

native_graphs_KKI = native_graphs["native_graphs_KKI"]
native_graphs_HNU1 = native_graphs["native_graphs_HNU1"]
native_graphs_BNU1 = native_graphs["native_graphs_BNU1"]
native_graphs_SWU4 = native_graphs["native_graphs_SWU4"]
#%%

# Concatenate everything together
graphs_HNU1 = final_concat(native_graphs_HNU1, cpac_graphs_HNU1, "HNU1")
graphs_SWU4 = final_concat(native_graphs_SWU4, cpac_graphs_SWU4, "SWU4")
graphs_BNU1 = final_concat(native_graphs_BNU1, cpac_graphs_BNU1, "BNU1")
# graphs_NKI = final_concat(native_graphs_NKI, cpac_graphs_NKI, "NKI")  # TODO
# graphs_KKI = final_concat(native_graphs_KKI, cpac_graphs_KKI, "KKI")  # TODO

# all graphs turned into (X, Y) matrix + target vector
# out_HNU1 = return_X_and_Y(graphs_HNU1)
# out_SWU4 = return_X_and_Y(graphs_SWU4)
# out_HNU1 = return_X_and_Y(graphs_BNU1)

#%%
out_HNU1 = return_X_and_Y(graphs_HNU1)
out_SWU4 = return_X_and_Y(graphs_SWU4)
out_BNU1 = return_X_and_Y(graphs_BNU1)

# do a check to make sure the order didn't get messed up
allgraphs = [graphs_BNU1, graphs_HNU1, graphs_SWU4]
alldata = [out_BNU1, out_HNU1, out_SWU4]

# graphs_SWU4[0]
return_X_and_Y(graphs_SWU4)[1][0]

# make sure all names in target vector match names in dictionary
assert np.all(
    out_SWU4[1]
    == list([name.split("_")[1] for name in list(OrderedDict(graphs_SWU4).keys())])
)
assert np.all(
    out_HNU1[1]
    == list([name.split("_")[1] for name in list(OrderedDict(graphs_HNU1).keys())])
)
assert np.all(
    out_BNU1[1]
    == list([name.split("_")[1] for name in list(OrderedDict(graphs_BNU1).keys())])
)
