import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

file = 'Data/SRIM/Xe300keV_in_ZrO2_Vacancies.txt'
file = 'Data/SRIM/Kr300keV_in_ZrO2_Vacancies.txt'
rho = 6.025
M_a = 123.218
N_a = 6.022e23

def D_d(x,a,b):
    damage_factor = (np.exp(-a*(x-b)**2)+1)
    return damage_factor

depth, ionvac, recoilvac = np.loadtxt(file,usecols=(0,1,2),unpack=True,encoding='cp437')
totvac = [sum(x) for x in zip(ionvac,recoilvac)]
totvac = [x*100000000 for x in totvac]
fluence = 1e17

N = N_a*rho/M_a
dpa = [x*fluence/(N) for x in totvac]
dpa = dpa[:30]
np.insert(depth, 0, 0)
norm = [i/sum(dpa) for i in dpa]
depth = depth[:30]
np.insert(depth, 0, 0)
depth = [d/10 for d in depth]
dx = 0.05
Nx =30
x = np.linspace(0,298,Nx)
Diff = np.array([D_d(i * dx,5,0.1) for i in range(Nx)])

Diff = [d/sum(Diff) for d in Diff]

popt,pcov = curve_fit(D_d,x,norm, p0=[20,0.1])
print(popt)

plt.plot(depth,norm,'--')
plt.plot(x,Diff)
# plt.plot(x,Diff)

# fig.legend(plts, labels= ['Xe-concentration','Kr-concentration','Xe-damage','Kr-damage'],loc = 'upper center',bbox_to_anchor = (0.67,0.94))
plt.xlabel('Depth [nm]')
plt.ylabel('Damage [dpa]')
plt.tight_layout()
plt.show()
