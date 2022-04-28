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

import os, sys

def E0(n,x):
    return -x*n + 0.5*n*(n-1)

def n0(x):
    return np.floor(x) + 1

xs = np.linspace(-1,1.9,100)
ns = [n0(x) for x in xs]
E0s = [E0(n,x) for (n,x) in zip(ns,xs)]
plt.plot(xs, ns, )
plt.plot(xs, E0s)
plt.show()