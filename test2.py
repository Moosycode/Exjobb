import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import datetime
import json
import os
import re

root = '/Users/niwi9751/Dropbox/Nils_files/Srim_Results/Zr300keV_in_UN_range.txt'
depth, height = np.loadtxt(root,usecols=(0,1),unpack=True,encoding='cp437')
rho = 14.05
M_a = 252
N_a = 6.022e23
fluence = 1e16
height = height*fluence
N = N_a*rho/M_a

conc = height/(height + N)

plt.plot(depth/10000,conc)
plt.show()