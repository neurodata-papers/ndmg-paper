# TODO: account for the fact that HNU1 is now `ssv` files

#%%
import os
import re
import csv
import pickle
import os.path as op
from copy import deepcopy as dc
from collections import OrderedDict
from pathlib import Path
import shutil

import networkx as nx
import scipy as sp
import numpy as np
from sklearn.preprocessing import normalize
from sklearn.neighbors import DistanceMetric
from scipy.linalg import svd
from scipy.linalg import norm
from scipy.stats import gaussian_kde

# from plotly.offline import download_plotlyjs, init_notebook_mode, iplot, plot
# from plotly import tools
# from plotly.graph_objs import *
import colorlover as cl

# import plotly_helper

# init_notebook_mode()
np.random.seed(12345678)  # for reproducibility, set random seed


#%%

#  dsets : ['NKI1', 'KKI2009', 'HNU1', 'BNU3', 'BNU1', 'SWU4']
dpath = Path("/Users/alex/Dropbox/NeuroData/ndmg-paper/data/multisite/")
dsets = [dset.name for dset in dpath.iterdir()]

# datadict: Dict of {'dset1': [path11, path12, ...], 'dset2': [path21, path22, ...], ...}
# datadict_strs: same thing, but values are strings rather than Paths
datadict = {
    dset.name: list(dset.iterdir())
    for dset in dpath.iterdir()
    if not dset.name == "avg"
}
datadict_strs = {dset: list(map(str, datadict[dset])) for dset in datadict}

# graph labels
labels = [
    "Betweenness Centrality",
    "Clustering Coefficient",
    "Normalized Degree",
    "Normalized Edge Weight",
    "Eigenvalue",
    "Locality Statistic-1",
    "Density",
]

# nsubs: number of unique subjects in datadict
# totalsubs: total unique subjects across datasets
# (that was a lot of code for that list of 5 numbers)
#        = [20, 57, 30, 227, 21]
nsubs = []
subrgx = re.compile(r"(sub-)([\d]*)")
for dataset in datadict:
    nsubs.append(
        len(
            set(
                [
                    re.search(subrgx, str(datadict[dataset][graph])).group(2)
                    for graph, _ in enumerate(datadict[dataset])
                ]
            )
        )
    )
totalsubs = sum(nsubs)
print("Datasets: " + "{}, ".format(list(zip(dsets, nsubs))))
#%% [markdown]
# Only run the below if you want to re-generate the averages for the data

#%%
N = 70  # desikan atlas


def avg_data(basepath, fs):

    op = "{}/avg".format(basepath)
    Path(op).mkdir(exist_ok=True, parents=True)

    # if originally, stats == `labels`, rut roh
    # I think fs was originally {'dataset1': ['something.betweenness_centrality.pkl', ...,]}
    stats = [
        "_".join(key.split(".")[0].split("_")[1:]) for key in fs[list(fs.keys())[0]]
    ]

    for stat in stats:
        print("Analyzing: {}".format(stat))
        if stat == "edge_weight":
            print(
                "Evaluating by proxy of mean connectome -- cannot average unequal lists"
            )
            continue
        # create empty dict
        average_dict = OrderedDict()
        for dset in fs:
            # load the data
            # note: for every dset, this currently only grabs the first file in it.
            fil = [fil for fil in fs[dset] if stat in fil][0]
            with open(fil, "rb") as f:
                dat = pickle.load(f)
            # f = open(fil)
            # dat = pickle.load(f)[stat]
            # f.close()

            # average it
            if stat == "degree_distribution":
                ipsi = np.zeros((len(dat["ipso_deg"].keys()), N))
                contra = np.zeros((len(dat["contra_deg"].keys()), N))
                total = np.zeros((len(dat["total_deg"].keys()), N))
                for idx, d in enumerate(dat["ipso_deg"]):
                    ipsi[idx, :] = dat["ipso_deg"][d]
                    contra[idx, :] = dat["contra_deg"][d]
                    total[idx, :] = dat["total_deg"][d]
                ipsi_avg = np.mean(ipsi, axis=0)
                contra_avg = np.mean(contra, axis=0)
                total_avg = np.mean(total, axis=0)
                avg = {
                    "ipso_deg": ipsi_avg,
                    "contra_deg": contra_avg,
                    "total_deg": total_avg,
                }
            elif stat == "study_mean_connectome":
                dat_reduced = np.array(
                    [
                        dat[i, j]
                        for i in np.arange(0, dat.shape[0])
                        for j in np.arange(i, dat.shape[1])
                    ]
                )
                avg = dat_reduced[np.where(dat_reduced > 0)]
            elif stat == "number_non_zeros":
                avg = np.array(dat.values())
            else:
                array = np.zeros((len(dat.keys()), N))
                for idx, d in enumerate(dat):
                    array[idx, :] = dat[d]
                avg = np.mean(array, axis=0)

            # place in dict in same format that we grabbed it
            average_dict[dset] = avg

        # save new dict of averages
        if stat == "degree_distribution":
            # reorganize
            tmp = average_dict.keys()[0]
            new = OrderedDict()
            for key in average_dict[tmp].keys():
                new[key] = {d: average_dict[d][key] for d in average_dict.keys()}
            average_dict = new
        elif stat == "study_mean_connectome":
            stat = "edge_weight"
        f = open("{}/avg_{}.pickle".format(op, stat), "wb")
        pickle.dump({stat: average_dict}, f)
        f.close()


#%%


# What the hell does this function do
avg_data(str(dpath), datadict_strs)
#%%
from IPython.display import HTML

# categorical 10
# cols = cl.scales['11']['qual']['Paired']
# cols = {d:cols[idx] for idx, d in enumerate(dsets)}
# HTML(cl.to_html(cols.values()))

# sequential orange-red 10
cols = dc(cl.scales["9"]["seq"]["OrRd"])
cols = cols[1:]
cols = cols + ["rgb(90,0,0)", "rgb(30,0,0)"]
# print(len(cols))
print(cols)
cols2 = OrderedDict()
for idx, d in enumerate(dsets):
    cols2[d] = cols[idx]
cols = cols2
print(cols.values())
HTML(cl.to_html(cols.values()))

# my categorical 10
# cols = ['rgba(228,26,28,{})', 'rgba(55,126,184,{})', 'rgba(77,175,74,{})',
#         'rgba(152,78,163,{})', 'rgba(166,86,40,{})', 'rgba(247,129,191,{})',
#         'rgba(255,204,0,{})', 'rgba(136,136,136,{})', 'rgba(55,224,169,{})',
#         'rgba(0,85,85,{})']
# cols = {d:cols[idx] for idx, d in enumerate(dsets)}

# sequential greys 13
# cols = ['#dedede', '#cdcdcd', '#bcbcbc', '#ababab', '#9a9a9a',
#         '#898989', '#787878', '#676767', '#565656', '#454545',
#         '#343434', '#232323', '#121212']

# sequential blue-green 13
# cols = cl.scales['9']['seq']['GnBu']
# cols = cols[2:4] + ['rgb(141,211,184)'] + cols[4:5] + ['rgb(101,182,192)'] + cols[5:6] +\
#        ['rgb(52,152,200)'] + cols[6:7] + ['rgb(21,121,181)'] + cols[7:] + ['rgb(6,55,100)', 'rgb(5,42,82)']
# HTML(cl.to_html(cols))


#%%
cols


#%%
normfactor = np.min(nsubs)
relsize = 1.4 * np.log2(np.array(nsubs)) / np.log10(totalsubs)
print(normfactor)
print(relsize)


#%%
fnames = [
    name
    for name in os.listdir(dpath + "/avg")
    if os.path.splitext(name)[1] == ".pickle" and "degree" in name
]
fnames = sorted(fnames)
paths = [os.path.join(dpath + "/avg", item) for item in fnames]
keys = ["_".join(n.split(".")[0].split("_")[1:]) for n in fnames]
for idx, curr in enumerate(paths):
    f = open(curr)
    dat = pickle.load(f)[keys[idx]]
    f.close()
    break

from itertools import chain, izip


d = dat["total_deg"]["BNU1"]
h1 = np.array(sorted(range(len(d[:35])), key=lambda k: d[k]))
h2 = h1[:] + 35
ordering = np.concatenate((h1, h2))
ordering = np.array(list(chain.from_iterable(izip(h1, h2))))
print(len(ordering))
print(ordering)


with open("./desikan.txt") as fil:
    rois = fil.read().split("\n")
rois
# dat = [0, 1, 2, 3, 8, 5, 6, 7]
# print sorted(range(len(dat)), key=lambda k: dat[k])


#%%
# ordering = np.arange(70)
# ordering2 = sorted(ordering, reverse=True)
# ordering[ordering2]


#%%
def plot_degrees(dats, name=None, ylab=None, xlab=None, hemi=False):
    data = list()
    if hemi:
        main = dats["ipso_deg"]
        contra = dats["contra_deg"]
    else:
        main = dats["total_deg"]
    al = 4.0 / len(main.keys())

    for idx, key in enumerate(dsets):
        lgth = len(main[key])
        data += [
            Scatter(
                x=np.linspace(1, lgth, lgth),
                y=main[key][ordering],
                line=Line(
                    color=cols[key],
                    #                                    width=relsize[idx],
                ),
                hoverinfo="x",
                name=key,
                showlegend=True,
                legendgroup=key,
            )
        ]
        if hemi:
            data += [
                Scatter(
                    x=np.linspace(1, lgth, lgth),
                    y=contra[key][ordering],
                    line=Line(
                        color=cols[key],
                        #                                        width=relsize[idx],
                        dash="dash",
                    ),
                    hoverinfo="x",
                    name=key,
                    showlegend=False,
                    legendgroup=key,
                )
            ]
    fig = Figure(data=data)
    return fig


def plot_series(
    stat, name=None, ylab=None, xlab=None, sort=False, leg=False, reverse=False
):
    data = list()
    for idx, key in enumerate(dsets):
        ys = stat[key]
        if sort:
            hov = "x"
            ys = np.array(sorted(ys, reverse=reverse))
        else:
            ys = ys[ordering]
            if idx == 0:
                hov = "text"
            else:
                hov = "none"
        data += [
            Scatter(
                x=np.linspace(1, len(ys), len(ys)),
                y=ys,
                line=Line(
                    color=cols[key],
                    #                                    width=relsize[idx],
                ),
                hoverinfo=hov,
                text=[rois[ordi] for ordi in ordering],
                # text= ordering+1,
                name=key,
                showlegend=False,
                legendgroup=key,
            )
        ]
    fig = Figure(data=data)
    return fig


def plot_jitter_scatter(stat, name=None, ylab=None, xlab=None):
    data = list()
    ys = [stat[key] for key in dsets]
    xs = np.linspace(-0.5, 0.5, len(ys))
    for idx, y in enumerate(ys):
        data += [
            Scatter(
                x=0.05 * np.random.rand(len(ys[idx])) + xs[idx],
                y=ys[idx],
                mode="markers",
                marker=Marker(
                    color=cols[dsets[idx]],
                    size=5,
                    #                             size=relsize[idx],
                    #                             opacity=0.2,
                ),
                hoverinfo="text",
                text=dsets[idx],
                name="{} ({})".format(dsets[idx], nsubs[idx]),
                showlegend=True,
                legendgroup=dsets[idx],
            )
        ]
    fig = Figure(data=data)
    return fig


# def plot_rugdensity(stat, name=None, ylab=None, xlab=None):
#     series = [stat[dset] for dset in dsets]
#     dens = gaussian_kde(series)
#     x = np.linspace(np.min(series), np.max(series), 100)
#     y = dens.evaluate(x)*np.max(series)

#     d_rug = Scatter(
#                 x=series,
#                 y=[0]*len(series),
#                 mode='markers',
#                 marker=Marker(
#                          color=[cols[dset] for dset in dsets],
#                          size=10,
#                        ),
#                 name=name,
#                 text=dsets,
#                 hoverinfo='text',
#                 showlegend=False
#           )

#     d_dens = Scatter(
#                 x=x,
#                 y=y,
#                 line=Line(
#                        color='rgba(0,0,0,0.6)'
#                      ),
#                 hoverinfo='x',
#                 name=name,
#                 showlegend=False
#            )
#     data = [d_dens, d_rug]
#     fig = Figure(data=data)
#     return fig


def make_panel_plot(basepath, dataset=None, atlas=None, log=True, hemispheres=False):
    fnames = [
        name for name in os.listdir(basepath) if os.path.splitext(name)[1] == ".pickle"
    ]
    fnames = sorted(fnames)
    paths = [os.path.join(basepath, item) for item in fnames]
    keys = ["_".join(n.split(".")[0].split("_")[1:]) for n in fnames]
    labels = [
        "Betweenness Centrality",
        "Clustering Coefficient",
        "Degree",
        "Edge Weight",
        "Eigenvalue",
        "Locality Statistic-1",
        "Number of Non-zeros",
        "Mean Connectome",
    ]

    traces = list(())
    for idx, curr in enumerate(paths):
        f = open(curr)
        dat = pickle.load(f)[keys[idx]]
        f.close()
        if keys[idx] == "number_non_zeros":
            #             fig = plot_rugdensity(dat)
            fig = plot_jitter_scatter(dat)
        elif keys[idx] == "edge_weight":
            edges = np.max([len(dat[i]) for i in dat.keys()])
            fig = plot_series(dat, sort=True)
        elif keys[idx] == "degree_distribution":
            if not hemispheres:
                tmp = dat["total_deg"]
                fig = plot_series(tmp, leg=True)
            else:
                fig = plot_degrees(dat, hemi=hemispheres)
                anno = [
                    dict(
                        x=dims / 3,
                        y=4 * dims / 7,
                        xref="x3",
                        yref="y3",
                        text="ipsilateral",
                        showarrow=False,
                        font=dict(color="rgb(0,0,0)", size=14),
                    ),
                    dict(
                        x=dims / 3,
                        y=3.7 * dims / 7,
                        xref="x3",
                        yref="y3",
                        text="contralateral",
                        showarrow=False,
                        font=dict(color="rgb(0,0,0)", size=14),
                    ),
                ]
        else:
            dims = len(dat.values()[0])
            if keys[idx] == "eigen_sequence":
                fig = plot_series(dat, sort=True, reverse=True)
            else:
                fig = plot_series(dat)
        traces += [pp.fig_to_trace(fig)]

    multi = pp.traces_to_panels(traces)
    for idx, curr in enumerate(paths):
        key = "axis%d" % (idx + 1)
        d = multi.layout["x" + key]["domain"]
        multi.layout["x" + key]["domain"] = [d[0], d[1] - 0.0125]
        multi.layout["x" + key]["zeroline"] = False
        multi.layout["y" + key]["zeroline"] = False
        multi.layout["y" + key]["title"] = labels[idx]
        multi.layout["x" + key]["title"] = "Node"
        multi.layout["x" + key]["nticks"] = 3
        multi.layout["y" + key]["nticks"] = 3
        if idx in [0, 1, 2, 3, 5]:
            multi.layout["x" + key]["range"] = [1, dims]
            multi.layout["x" + key]["tickvals"] = [1, dims / 2, dims]
            if idx in [2]:
                if hemispheres:
                    multi.layout["annotations"] = anno
            elif log:
                multi.layout["y" + key]["type"] = "log"
                multi.layout["y" + key]["title"] += " (log scale)"
        if idx in [3]:
            multi.layout["x" + key]["range"] = [1, edges]
            multi.layout["x" + key]["tickvals"] = [1, edges / 2, edges]
            multi.layout["x" + key]["title"] = "Edge"
        if idx in [4]:
            multi.layout["x" + key]["range"] = [1, dims]
            multi.layout["x" + key]["tickvals"] = [1, dims / 2, dims]
            multi.layout["x" + key]["title"] = "Dimension"
        if idx in [6]:
            multi.layout["x" + key]["title"] = "Dataset"
            multi.layout["y" + key]["range"] = [0, 1500]
            multi.layout["x" + key]["tickvals"] = [0]
            multi.layout["x" + key]["ticktext"] = [""]
        if idx in [7]:
            multi.layout["y" + key]["title"] = None
            multi.layout["x" + key]["title"] = labels[idx]
            multi.layout["y" + key]["autorange"] = "reversed"
            multi.layout["x" + key]["tickvals"] = [0, dims / 2 - 1, dims - 1]
            multi.layout["y" + key]["tickvals"] = [0, dims / 2 - 1, dims - 1]
            multi.layout["x" + key]["ticktext"] = [1, dims / 2, dims]
            multi.layout["y" + key]["ticktext"] = [1, dims / 2, dims]
            if log:
                multi.layout["x" + key]["title"] += " (log10)"

    if dataset is not None and atlas is not None:
        if atlas == "desikan":
            atlas = atlas.capitalize()
        tit = dataset + " (" + atlas + " parcellation)"
    else:
        tit = None

    multi.layout["showlegend"] = True
    multi.layout["legend"]["orientation"] = "v"
    multi.layout["legend"]["xanchor"] = "x7"
    multi.layout["legend"]["yanchor"] = "y8"
    multi.layout["legend"]["x"] = 0.78
    multi.layout["legend"]["y"] = -0.17
    multi.layout["legend"]["tracegroupgap"] = 0
    multi.layout["legend"]["font"] = {"size": 14}

    anno = [
        dict(
            x=0.97,
            y=0.46,
            xref="paper",
            yref="paper",
            text="Dataset (Number of Subjects)",
            showarrow=False,
            font=dict(color="#28282e", size=14),
        )
    ]
    multi.layout.annotations = anno
    multi.layout["title"] = tit
    return multi


#%%
atlas = "Desikan"
multi = make_panel_plot(
    dpath + "/avg", hemispheres=False, atlas=atlas, dataset="Multiple Datasets"
)

# iplot(multi)
plot(
    multi,
    validate=False,
    filename="/Users/gkiar/code/ocp/ndmg-paper/code/multisite_graphs/multisite_sorted_dset.html",
)


#%%

