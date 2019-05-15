#%%
r"""Get all cpac graphs for (HNU1, BNU1, SWU4, KKI2009, NKI1) | graph in cpac,
for every cpac graph, concatenate it with its corresponding native-space graph,
put them all in an output folder to compute discriminability.
"""
# TODO: figure out why the fuck the shape of the native-space graphs are all over the place

import os
from pathlib import Path
from collections import OrderedDict
from functools import partial

import pandas as pd
import numpy as np
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
def make_dict_of_graphs(dataset: "Path"):
    """
    From a Path of a directory containing a bunch of graph .ssv files,
    return: dict, {subj_ses: np.array}
    """
    import_graph = partial(graspy.utils.import_edgelist, extension="ssv", delimiter=" ")

    # more clear than a dict comprehension
    output = dict()
    for pgraph in dataset.iterdir():
        if pgraph.suffix == ".ssv":
            graph = import_graph(pgraph)
            if np.any(graph):  # exclude graphs that are all 0s
                output[pgraph.stem] = graph
    return output


# {dataset.stem: {subj_ses: np.array}}
cpac_graphs = {dataset.stem: make_dict_of_graphs(dataset) for dataset in fcpac_graphs}
native_graphs = {
    dataset.stem: make_dict_of_graphs(dataset) for dataset in fnative_graphs
}
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


def combined_graph_dict(cpac, native):
    """
    takes two dicts containing sub# and ses#: array,
    and for each matching key,
    concatenates the values using `concatenate_cpac_and_native`
    """
    pass


#%%
graph = "sub-906_ses-2_dwi_desikan_space-MNI152NLin6_res-2x2x2_measure-spatial-ds_adj"
parts = graph.split("_")[0:2]  # ['sub-#', 'ses-#']
pieces = [int(i.split("-")[1]) for i in parts]  # e.g., [849, 2]
sub, ses = pieces
corresponding_cpac = f"{sub}_session_{ses}"
concatenate_cpac_and_native(
    cpac_graphs_KKI[corresponding_cpac], native_graphs_KKI[graph]
)

# for graph in native_graphs_KKI:
#     # turn graph into the equivalent thing in cpac_graphs_KKI
#     parts = graph.split('_')[0:2]  # ['sub-#', 'ses-#']
#     pieces = [int(i.split('-')[1]) for i in parts]  # e.g., [849, 2]
#     sub, ses = pieces
#     corresponding_cpac = f"{sub}_session_{ses}"

#     # gonna need to turn this into a dictionary appending thing
#     concatenate_cpac_and_native(cpac[corresponding_cpac], native_graphs_KKI[graph])

#%%
# check lengths of graphs

# import_graph = partial(graspy.utils.import_edgelist, extension="ssv", delimiter=" ")
# print([np.shape(graph) for graph in import_graph(fnative_graphs[0])])
# print([np.shape(graph) for graph in import_graph(fnative_graphs[1])])
# print([np.shape(graph) for graph in import_graph(fnative_graphs[2])])
# print([np.shape(graph) for graph in import_graph(fnative_graphs[3])])
# print([np.shape(graph) for graph in import_graph(fnative_graphs[4])])

print([np.shape(graph) for graph in import_graph(fcpac_graphs[0])])

