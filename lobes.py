#!/usr/bin/env python
# coding: utf-8
import numpy as np

from matplotlib.pyplot import cm, colorbar
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['figure.figsize'] = ((3+3/8), (3+3/8))
plt.rc('text.latex', preamble=r'\usepackage{bm,braket}')

import matplotlib.pyplot as plt
import pandas as pd

import sys, os, re
"""
    detect files in folder pwd_path that match find_str
"""

df = pd.read_csv("observables.csv")
obs = np.sort(np.unique(df["obs"]))
evs = np.sort(np.unique(df["ev"]))
ks = np.sort(df.keys())

print(obs)
lbl = ["\\braket{\\hat a}","\\braket{\\hat a^\\dag}","\\braket{\\hat n}"]
for (id,o) in enumerate(obs):
    for ev in evs:
        for v in ["val_re", "val_im"]:
            data = df[df["obs"] == o]
            data = data[data["ev"] == ev]
            data = data.groupby(["u", "t", "mu"]).mean().reset_index()

            pivotted = np.asarray(data.pivot("mu","t",v))
            extent = [np.min(data["t"]),np.max(data["t"]),np.min(data["mu"]),np.max(data["mu"])]
            fig, ax = plt.subplots(1,1)
            imag = ax.imshow(pivotted, extent=extent, cmap='RdBu_r', interpolation='bicubic', aspect='auto', origin='lower')
            # levels = [0.97,0.995,1.005,1.03,1.97,1.995,2.005,2.03]
            # if ob == 'N':
            #     ax.contour(pivotted, levels=levels, corner_mask=False, colors=['gray','black','black','gray'], extent=extent, origin='lower')
            fig.colorbar(imag, extend='both', pad=0.05, label=f'${lbl[id]}$')
            ax.set_xlabel("$-t/U$")
            ax.set_ylabel("$-\mu/U$")
            plt.tight_layout()
            # plt.show()
            plt.savefig(f"mott_lobes_{o}_{ev}_{v}.jpg", bbox_inches='tight', dpi=600, pad_inches=0.0)
