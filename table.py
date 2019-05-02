#%%
import pandas as pd
from numpy import NaN

#%%

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
