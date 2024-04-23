import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import datetime
import json
import os
import re

roots = ['/Users/niwi9751/Dropbox/Nils_files/Srim_Results/Zr330keV_in_UN_range.txt']
Names = ['Zirconium']


i = 0
for root in roots:
    depth, height = np.loadtxt(root,usecols=(0,1),unpack=True,encoding='cp437')
    rho = 14.05
    M_a = 252
    N_a = 6.022e23
    fluence = 3.73e16
    height = height*fluence
    N = N_a*rho/M_a
    conc = height/(height + N)
    plt.plot(depth/10000,conc,label=Names[i])
    i = i+1
plt.grid()
plt.legend()
plt.title('Implanted concentration')
plt.xlabel('Depth [micrometer]')
plt.ylabel('Concentration [at. fraction]')


roots2 = ['/Users/niwi9751/Dropbox/Nils_files/Srim_Results/Zr330keV_in_UN_Vacancies.txt']

plt.figure()
i = 0
for root in roots2:
    depth, ionvac, recoilvac = np.loadtxt(root,usecols=(0,1,2),unpack=True,encoding='cp437')
    totvac = [sum(x) for x in zip(ionvac,recoilvac)]
    totvac = [x*100000000 for x in totvac]
    rho = 5.68
    M_a = 123
    N_a = 6.022e23
    fluence = 1e17
    N = N_a*rho/M_a
    dpa = [x*fluence/(100*N) for x in totvac]
    plt.plot(depth/10000,dpa,'-.',label=Names[i])
    i = i+1

plt.grid()
plt.legend()
plt.title('Implanted damage')
plt.xlabel('Depth [micrometer]')
plt.ylabel('dpa')
plt.show()
