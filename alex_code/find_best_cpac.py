import os
from pathlib import Path
from collections import OrderedDict
from collections import Iterable
from functools import partial
from functools import reduce
import warnings

import pandas as pd
import numpy as np

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
