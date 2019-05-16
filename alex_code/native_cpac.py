#%%
r"""Get all cpac graphs for (HNU1, BNU1, SWU4, KKI2009, NKI1) | graph in cpac,
for every cpac graph, concatenate it with its corresponding native-space graph,
put them all in an output folder to compute discriminability.
"""
# TODO: figure out why the fuck the shape of the native-space graphs are all over the place
# TODO: ok figured out the above. Now: edit `make_dict_of_graphs` to return graphs with the right number of nodes.

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

if not Path(os.getcwd()).name == "alex_code":
    os.chdir("alex_code")
#%%
cpac = pd.read_csv("discr_fmri_results.csv")
cpac
#%%
cpac[
    (cpac["Reg"] == "F")
    & (cpac["FF"] == "N")
    & (cpac["Scr"] == "N")
    & (cpac["GSR"] == "G")
]
#%%
# figure out which index has the best overall parcellation
cpac = cpac[cpac.Parcellation == "D"][["Dataset", "discr"]].sort_values(
    by="discr", ascending=False
)
#%%
cpac.reset_index(inplace=True)

#%%
# for each dataset,
# sort by discriminability,
# then see what the mode index is.
best_discrim = cpac.sort_values(["Dataset", "discr"], ascending=False)
best_discrim

# ok that won't work out, just use HNU1
# use that I guess
cpac[cpac["Dataset"] == "HNU1"].sort_values("discr", ascending=False)  # it's FSL11192

#%%
cpac = pd.read_csv("discr_fmri_results.csv")
cpac.reset_index(inplace=True)
cpac.rename(columns={"index": "param"}, inplace=True)
#%%
hnu1_cpac = cpac[cpac.Dataset == "HNU1"]

# Reg: F
# FF: N
# Scr: S
# GSR: G
# xfm: P
hnu1_cpac[hnu1_cpac["Parcellation"] == "D"].sort_values("discr", ascending=False)

#%%
files = [
    "ANT_frf_nsc_gsr_aal",
    "ANT_frf_nsc_ngs_des",
    "ANT_frf_scr_ngs_aal",
    "ANT_nff_nsc_gsr_des",
    "ANT_nff_scr_gsr_aal",
    "ANT_nff_scr_ngs_des",
    "FSL_frf_nsc_ngs_aal",
    "FSL_frf_scr_gsr_des",
    "FSL_nff_nsc_gsr_aal",
    "FSL_nff_nsc_ngs_des",
    "FSL_nff_scr_ngs_aal",
    "ANT_frf_nsc_gsr_cc2",
    "ANT_frf_nsc_ngs_hox",
    "ANT_frf_scr_ngs_cc2",
    "ANT_nff_nsc_gsr_hox",
    "ANT_nff_scr_gsr_cc2",
    "ANT_nff_scr_ngs_hox",
    "FSL_frf_nsc_ngs_cc2",
    "FSL_frf_scr_gsr_hox",
    "FSL_nff_nsc_gsr_cc2",
    "FSL_nff_nsc_ngs_hox",
    "FSL_nff_scr_ngs_cc2",
    "ANT_frf_nsc_gsr_des",
    "ANT_frf_scr_gsr_aal",
    "ANT_frf_scr_ngs_des",
    "ANT_nff_nsc_ngs_aal",
    "ANT_nff_scr_gsr_des",
    "FSL_frf_nsc_gsr_aal",
    "FSL_frf_nsc_ngs_des",
    "FSL_frf_scr_ngs_aal",
    "FSL_nff_nsc_gsr_des",
    "FSL_nff_scr_gsr_aal",
    "FSL_nff_scr_ngs_des",
    "ANT_frf_nsc_gsr_hox",
    "ANT_frf_scr_gsr_cc2",
    "ANT_frf_scr_ngs_hox",
    "ANT_nff_nsc_ngs_cc2",
    "ANT_nff_scr_gsr_hox",
    "FSL_frf_nsc_gsr_cc2",
    "FSL_frf_nsc_ngs_hox",
    "FSL_frf_scr_ngs_cc2",
    "FSL_nff_nsc_gsr_hox",
    "FSL_nff_scr_gsr_cc2",
    "FSL_nff_scr_ngs_hox",
    "ANT_frf_nsc_ngs_aal",
    "ANT_frf_scr_gsr_des",
    "ANT_nff_nsc_gsr_aal",
    "ANT_nff_nsc_ngs_des",
    "ANT_nff_scr_ngs_aal",
    "FSL_frf_nsc_gsr_des",
    "FSL_frf_scr_gsr_aal",
    "FSL_frf_scr_ngs_des",
    "FSL_nff_nsc_ngs_aal",
    "FSL_nff_scr_gsr_des",
    "ANT_frf_nsc_ngs_cc2",
    "ANT_frf_scr_gsr_hox",
    "ANT_nff_nsc_gsr_cc2",
    "ANT_nff_nsc_ngs_hox",
    "ANT_nff_scr_ngs_cc2",
    "FSL_frf_nsc_gsr_hox",
    "FSL_frf_scr_gsr_cc2",
    "FSL_frf_scr_ngs_hox",
    "FSL_nff_nsc_ngs_cc2",
    "FSL_nff_scr_gsr_hox",
]

fileslist = [fil.split("_") for fil in files]
df_fileslist = pd.DataFrame(np.array(fileslist))
df_fileslist.head()

#%%
hnu1_cpac[hnu1_cpac.Parcellation == "D"].head().sort_values("discr", ascending=False)

#%%
possibilities = df_fileslist[df_fileslist[4] == "des"].head().reset_index(drop=True)
possibilities

#%%
# Reg: F
# FF: N
# Scr: S
# GSR: G
# xfm: P
possibilities
possibilities.iloc[[1, 2, 4]]

"_".join(list(possibilities.iloc[2]))

#%%
best = dict(hnu1_cpac.sort_values(by="discr", ascending=False).iloc[0])
best


#%%
df_fileslist = pd.DataFrame(fileslist)
df_fileslist

#%%
df_fileslist = df_fileslist[df_fileslist[4] == "des"]
df_fileslist

#%%
# find the value in `df_fileslist` that corresponds to `best`
# {'param': 'FSL11182',
#  'Dataset': 'HNU1',
#  'Reg': 'F',
#  'FF': 'N',
#  'Scr': 'S',
#  'GSR': 'G',
#  'Parcellation': 'C',
#  'xfm': 'P',
#  'nsub': 30,
#  'nses': 10,
#  'nscans': 300,
#  'nroi': 201,
#  'discr': 0.991390804597701}

# file corresponding to `best`
bestfile_list = df_fileslist[
    (df_fileslist[0] == "FSL")
    & (df_fileslist[1] == "nff")
    & (df_fileslist[2] == "scr")
    & (df_fileslist[3] == "gsr")
]
"_".join(list(bestfile_list.iloc[0, :]))
# '_'.join((str(i) for i in bestfile_list))
#%%
""" 
Conclusion: we want `FSL_nff_scr_gsr_des`
"""

# now that I have the graphs,
# define a function that,
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


# grab a function that computes discriminability on all these graphs.
# then, run this on each dataset.

#%%
# get DFXSG pipeline graphs
# NKI1 isn't in CPAC
map_to_files = OrderedDict({"F": "FSL", "X": "nff", "S": "scr", "G": "gsr", "D": "des"})
"_".join([map_to_files[i] for i in list("FXSGD")])

# Concatenation
# ------------------------------------------------------------ #
#%%
# grab separate graphs
p = Path.cwd()
p = p.parent / "data" / "graphs"
fcpac_graphs = [folder for folder in p.iterdir() if "cpac" in folder.name]
fnative_graphs = [folder for folder in p.iterdir() if "native" in folder.name]
#%%
# partial function to make graph importing easier, then
# make an [array of 'thing.name': np.array for array in cpac_graphs]
def make_dict_of_graphs(dataset: "Path"):  # TODO: depracated.
    """
    From a Path of a directory containing a bunch of graph .ssv files,
    return: dict, {subj_ses: np.array}
    """
    import_graph = partial(
        import_edgelist, extension="ssv", delimiter=" ", output_dict=True
    )

    list_of_graphs = import_graph(dataset)

    # more clear than a dict comprehension
    output = dict()
    for pgraph in dataset.iterdir():
        if pgraph.suffix == ".ssv":
            graph = import_graph(pgraph)
            if np.any(graph):  # exclude graphs that are all 0s
                output[pgraph.stem] = graph
    return output


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
        # TODO: make sure these are sorted properly.
        filenames = [filepath.stem for filepath in files]
        dict_out = dict(zip(filenames, out))
        return dict_out

    else:
        return out


dataset = "/Users/alex/Dropbox/NeuroData/ndmg-paper/data/graphs/native_graphs_HNU1"
test_path = "/Users/alex/Dropbox/NeuroData/ndmg-paper/data/graphs/native_graphs_HNU1/sub-0025427_ses-1_dwi_desikan_space-MNI152NLin6_res-2x2x2_measure-spatial-ds_adj.ssv"
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

#%%
# test to make sure dataset_dict sorts graphs properly
def return_sorted_graph(graph_file):
    graph = nx.read_weighted_edgelist(test_path, nodetype=int, delimiter=" ")
    vertices = np.arange(1, 71)
    out = nx.to_numpy_array(graph, nodelist=vertices, dtype=np.float)
    return out


np.all(
    dataset_dict[
        "sub-0025427_ses-1_dwi_desikan_space-MNI152NLin6_res-2x2x2_measure-spatial-ds_adj"
    ]
    == return_sorted_graph(test_path)
)  # return True


#%%
# {dataset.stem: {subj_ses: np.array}}
# everything is now correctly the right shape
cpac_graphs = {dataset.stem: import_graph(dataset) for dataset in fcpac_graphs}
native_graphs = {dataset.stem: import_graph(dataset) for dataset in fnative_graphs}
#%%
# for every dataset in cpac_graphs and native_graphs,
# for every array in that dataset,
# check if the subject is the same,
# and the session is the same.
# if both are, concatenate the graphs.

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
# graph = (
#     "sub-0025427_ses-1_dwi_desikan_space-MNI152NLin6_res-2x2x2_measure-spatial-ds_adj"
# )
# parts = graph.split("_")[0:2]
# print(parts)  # ['sub-#', 'ses-#']
# pieces = [i.split("-")[1] for i in parts]
# print(pieces)  # e.g., [849, 2]
# sub, ses = pieces
# corresponding_cpac = f"HNU1_{sub}_session_{ses}_correlation"
# corresponding_cpac
# graph = concatenate_cpac_and_native(
#     cpac_graphs_HNU1[corresponding_cpac], native_graphs_HNU1[graph]
# )
# np.shape(graph)

concatenated_graphs_HNU1 = []
dataset = "HNU1"
native_graphs_HNU1

# TODO: maybe turn this into a function? could be useful


def final_concat(native_graphs, cpac_graphs, dataset):
    concatenated_graphs = []
    for fgraph in native_graphs:
        # turn graph into the equivalent thing in cpac_graphs_KKI
        parts = fgraph.split("_")[0:2]  # ['sub-#', 'ses-#']
        pieces = [i.split("-")[1] for i in parts]  # e.g., [849, 2]
        sub, ses = pieces
        corresponding_cpac = f"{dataset}_{sub}_session_{ses}_correlation"

        # gonna need to turn this into a dictionary appending thing
        graph = concatenate_cpac_and_native(
            cpac_graphs[corresponding_cpac], native_graphs[fgraph]
        )
        key = f"sub_{sub}_session_{ses}"
        try:
            concatenated_graphs.append((key, graph))
        except:
            print(
                f"parts: {parts}, pieces: {pieces}, corresponding_cpac: {corresponding_cpac}"
            )
            pass
        return dict(concatenated_graphs)


#%%
# Whoa what the fuck, it works
graphs_HNU1 = final_concat(native_graphs_HNU1, cpac_graphs_HNU1, "HNU1")
graphs_SWU4 = final_concat(native_graphs_SWU4, cpac_graphs_SWU4, "SWU4")
graphs_BNU1 = final_concat(native_graphs_BNU1, cpac_graphs_BNU1, "BNU1")
# graphs_NKI = final_concat(native_graphs_NKI, cpac_graphs_NKI, "NKI")  # TODO
# graphs_KKI = final_concat(native_graphs_KKI, cpac_graphs_KKI, "KKI")  # TODO

