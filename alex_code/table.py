# TODO: sort within each category highest to lowest
# TODO: KKI and NKI native+cpac numbers
# TODO: proper sig figs
#%%
import os
import shutil
from pathlib import Path
from collections import OrderedDict

import pandas as pd
from numpy import NaN
import rpy2.robjects as robjects

if not Path(os.getcwd()).name == "alex_code":
    os.chdir("alex_code")
#%%
def latexify(df):
    df.to_latex(buf="table_latex.txt", na_rep="-")


# all datasets we've run with diffusion + func
dif_func = ["HNU1", "BNU1", "SWU4", "BNU3"]

# all datasets run with only diffusion
dif_only = ["KKI2009", "NKI1", "Temp114", "NKI-ENH", "Temp255", "PING", "MRN1313"]

# all datasets run with only functional
func_only = [
    "IPCAS6",
    "IPCAS8",
    "SWU1",
    "IPCAS5",
    "SWU3",
    "XHCUMS",
    "UWM",
    "NYU1",
    "SWU2",
    "IPCAS1",
    "IPCAS2",
    "IBATRT",
    "MRN1",
    "BNU2",
]

# all datasets
dsets = dif_func + dif_only + func_only
dsets = {i: [] for i in dsets}

# make table skeleton
table = pd.DataFrame(
    columns=(
        "Study",
        "dMRI",
        "Native-Space",
        "fMRI",
        "CPAC",
        "dMRI + fMRI",
        "Native-Space + CPAC",
    )
)
table["Study"] = list(dsets.keys())
table.set_index("Study", inplace=True)
table
#%%
# list of discriminability values, with [0] being IPCAS6
fMRI_discrim_list = [
    0.994,
    0.958,
    0.935,
    0.819,
    0.986,
    0.823,
    0.849,
    0.931,
    0.874,
    0.893,
    0.868,
    0.974,
    0.859,
    0.863,
]

# pleased with how efficient this is
table["fMRI"].loc["IPCAS6":"BNU2"] = fMRI_discrim_list


# manually add values for which it isn't easy to add programatically
# dMRI, Native-Space, fMRI, CPAC, dMRI + fMRI, Native-Space + CPAC
table.loc["HNU1"] = [0.993, 0.996, 0.956, NaN, 0.993, NaN]
table.loc["BNU1"] = [0.984, 0.723, 0.906, NaN, 0.990, NaN]
table.loc["SWU4"] = [0.884, 0.718, 0.891, NaN, 0.903, NaN]
table.loc["NKI1"] = [0.984, 0.984, NaN, NaN, NaN, NaN]
table.loc["KKI2009"] = [1.00, 1.00, NaN, NaN, NaN, NaN]

table
#%%
# grab cpac data using rpy2 and stick it into a pandas dataframe
rdf = robjects.r("dat <- readRDS('discr_fmri_results.rds')")
rdf.to_csvfile("discr_fmri_results.csv")
df_cpac = pd.read_csv("discr_fmri_results.csv")

#%%
# subset cpac numbers to be what we want for the table
Reg = "F"
FF = "N"
Scr = "S"
Gsr = "G"
Parcellation = "D"
xfm = "P"
df_cpac = df_cpac[
    (df_cpac["Reg"] == Reg)
    & (df_cpac["FF"] == FF)
    & (df_cpac["Scr"] == Scr)
    & (df_cpac["GSR"] == Gsr)
    & (df_cpac["Parcellation"] == Parcellation)
    & (df_cpac["xfm"] == xfm)
]
df_cpac
#%%
# Take the average discriminability grouped by dataset
df_cpac.set_index("Dataset", drop=True, inplace=True)
df_cpac = df_cpac.discr
df_cpac
#%%
# Rename in preparation for merge, merge, then clean
df_cpac.index.name = "Study"
df_cpac.name = "CPAC"
table = table.merge(df_cpac, how="left", right_index=True, left_index=True)
table["CPAC_x"], table["CPAC_y"] = table["CPAC_y"], table["CPAC_x"]
table = table.rename(columns={"CPAC_x": "CPAC"}).drop("CPAC_y", axis=1)
table.rename(
    columns={"dMRI": "ndmg-d", "fMRI": "ndmg-f", "dMRI + fMRI": "ndmg-d + ndmg-f"},
    inplace=True,
)
table
#%%
new_native_space_numbers = {"NKI1": 0.918, "HNU1": 0.989, "BNU1": 0.999, "SWU4": 0.815}
table["Native-Space"].update(pd.Series(new_native_space_numbers))

#%%
table.loc["KKI2009", "Native-Space"] = NaN


#%%
table.loc["BNU1", "Native-Space"] = 0.989
table.loc["HNU1", "Native-Space"] = 0.989
table.loc["NKI", "Native-Space"] = 0.918
table.loc["SWU4", "Native-Space"] = 0.815
table.loc["KKI2009", "Native-Space"] = np.NaN
table
#%%
# add native+cpac numbers
table.loc["BNU1", "Native-Space + CPAC"] = 0.866
table.loc["HNU1", "Native-Space + CPAC"] = 0.831
table.loc["SWU4", "Native-Space + CPAC"] = 0.568

# drop weird extra NKI row
table.drop("NKI", axis=0, inplace=True)

#%%
# move all of NaN rows to the bottom
nans = table.isnull().all(axis=1)
notnans = ~nans
table = table[notnans].append(table[nans])
table
#%%
# TODO: sort within each category highest to lowest

# 
#%%
# latexify(table)
