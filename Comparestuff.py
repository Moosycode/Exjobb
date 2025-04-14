import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin
import pandas as pd

def rebin(data,x_res):
    data = [(data[1::2][i] + data[::2][i])/2 for i in range(len(data[1::2]))]
    x_res = x_res[1::2]
    return data, x_res

def rebin_array(data, n_new_bins, method='mean'):
    data = np.array(data)
    old_len = len(data)
    trim_len = (old_len // n_new_bins) * n_new_bins
    trimmed = data[:trim_len]
    
    factor = trim_len // n_new_bins
    reshaped = trimmed.reshape(n_new_bins, factor)

    if method == 'mean':
        return reshaped.mean(axis=1)
    elif method == 'sum':
        return reshaped.sum(axis=1)
    else:
        raise ValueError("Method must be 'mean' or 'sum'")


# Diffusion coefficient dependant on temperature
def D(D0, Ea, Temp):
    D = D0*np.exp(-Ea/(k*Temp))

    return D
    
def histog(data, length):
    L = length*10000 #Make sure unit is correct
    n, bins = np.histogram(data, bins = 100,range=(0,L))
    width = bins[1]-bins[0]
    y = n/sum(n) #Normalize
    return width,y

def hist_integral(n, width):
    n = [item*width for item in n]
    return sum(n)#Definition of integrals :))

def V(x):
    return a
#Constants----------------------------
k = 8.6e-5 # boltzmann constant [ev/K]
Na = 6.022e23 # avogadros number [atoms/mole]
#-------------------------------------

#GENERAL FILEPATHS-----------------------------------------------------------------
pre_data_path = 'Data/ToFERDA/Data/pre-anneal_Xe-imp.csv'
post_data_path = 'Data/ToFERDA/Data/post-anneal_Xe-imp.csv'
# pre_data_path = 'Data/ToFERDA/Data/pre-anneal_Kr-imp_corrected.csv'
# post_data_path = 'Data/ToFERDA/Data/post-anneal_Kr-imp_corrected.csv'
srim_data = 'Data/SRIM/Xe300keV_in_ZrO2_vacancies.txt'
#---------------------------------------------------------------------------------
depth, ionvac, recoilvac = np.loadtxt(srim_data,usecols=(0,1,2),unpack=True,encoding='cp437')
totvac = [sum(x) for x in zip(ionvac,recoilvac)]

#Global Parameters-----------------------------------------------------------------------
Times_in = 9
Temp = 1473.15 #Target emperature [K] 
fluence = 1e17# Input fluence of implantation [atoms/cm^2]
Integrate = False
Concentrations = []#Result list
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=16)
#---------------------------------------------------------------------------------
T = Times_in*3600 #Convert to seconds
#---------------------------------------------------------------------------------

#Dictionary with needed values of each element
elementdict = {
    'Fe_ZrO2':{'D0':2.26e-6,'Ea':2.3, 'rho':6.025, 'Ma': 123.218}, #Data from Springer
    'Kr_ZrO2':{'D0':8.11e-7,'Ea':2.53, 'rho':6.025, 'Ma': 123.218},
    'Xe_ZrO2':{'D0':1.83e-6,'Ea':2.91, 'rho':6.025, 'Ma': 123.218},
    'Zr_UN':{'D0':6.9e-7,'Ea':2.7, 'rho':14.05, 'Ma': 252.036}, 
    'Kr_UN':{'D0':8.11e-7,'Ea':2.53, 'rho':14.05, 'Ma': 252.036},
    'Xe_UN':{'D0':1.83e-6,'Ea':2.91, 'rho':14.05, 'Ma': 252.036},
    'Zr_UO2':{'D0':6.9e-7,'Ea':2.7, 'rho':10.6, 'Ma': 270.02},
    'Kr_UO2':{'D0':8.11e-7,'Ea':2.53, 'rho':10.6, 'Ma': 270.02}
}



vals = np.linspace(1,1000,5)
vals = [1,10,100,500]
plt.figure(figsize=(8,6))
for v in vals:
    # Parameters------------------------------------------------------
    element = 'Xe_ZrO2'
    sample = 'Xe'
    df_pre = pd.read_csv(pre_data_path)
    df_pre = df_pre[df_pre['depth'] <= 3000]
    c_pre = df_pre[sample]
    err_pre = df_pre[f'{sample}_sd']
    tot_err_pre = (sum([e**2 for e in err_pre]))**0.5
    print('Avg conc pre: ' + str(sum(c_pre)/len(c_pre)) + ' error: ' + str(tot_err_pre/len(err_pre)))
    x_pre = df_pre['depth']
    df_post = pd.read_csv(post_data_path)
    df_post = df_post[df_post['depth'] <= 3000]
    c_post = df_post[sample]
    err_post = df_post[f'{sample}_sd']
    tot_err_post = (sum([e**2 for e in err_post]))**0.5
    print('Avg conc post: ' + str(sum(c_post)/len(c_post)) + ' error: ' + str(tot_err_post/len(err_post)) )
    x_post = df_post['depth']
  
    optivals = [84.43171634990631, 2.4644511371221833, 5.700699145455733e-05] #Kr
    optivals= [567.2290114003844, 3.107999491728278, 5.728181804784536e-05] #Xe
    #Constants--------------------------------------------------------
    
    D0 = optivals[0]
    Ea = optivals[1]
    a = optivals[2]
    rho = elementdict[element]['rho']# density of target [g/cm^3]
    m_a = elementdict[element]['Ma']# atomic mass of target in [g/mole]
    n_atoms = rho*Na/m_a #atomic density of target [atoms/cm^3]
    Nx = len(x_pre)  # Number of spatial points per micrometer
    Nt = int(T/60)*2 # Number of time steps, can be anything really, code finds this for you but do not start lower than this.
    x_pre = [3*1e21*x/(n_atoms) for x in x_pre] #Convert to micrometer
    
    Extendby = 2
    L = int(x_pre[-1])*Extendby
    studyL = int(x_pre[-1])
    dx = L / (L*Nx - 1) # Spatial step size
    dt = T / Nt # Time step size
     #integration--------------------------
    
    binwidth = (x_pre[1]-x_pre[0])*1e-7
    Integrate=True
    if Integrate:
        #Integrate over total length
        I_inital = hist_integral(c_pre,binwidth)
        I_diffused = hist_integral(c_post, binwidth)
        ratio = I_inital/I_diffused
        print(f'Total integral of inital concentration: ')
        print(I_inital*n_atoms)
        print(f'Total integral of diffused concentration: ')
        print(I_diffused*n_atoms)
        print(f'Ratio between total integrals: ')
        print(ratio)
        print()
    #-----------------------------------------------------------------
    # Create spatial grid
    x = np.linspace(0, L, Nx*Extendby)
    # Initialize solution matrix
    C = np.zeros((Nt, Nx))
    
    # Apply initial condition 
    C[0, :] = c_pre
    C = np.hstack((C, np.zeros((Nt,(Extendby-1)*Nx)))) #Add zeros to desired length, SRIM length is only 1 micron usually
    rebinned = rebin_array(totvac, Nx*Extendby, 'sum')
    d = v
    rebinned = [d*x/sum(rebinned) + 1  for x in rebinned]
    D_x = np.array([D(D0, Ea, Temp)*rebinned[i] for i in range(Nx*Extendby)])
    V_x = np.array([V(i * dx) for i in range(Nx*Extendby)]) 
    # Time-stepping loop
    for n in range(0, Nt - 1):         
        for i in range(1, Extendby*Nx - 1): #Update interior points 
            D_left = D_x[i-1]
            D_right = D_x[i+1]
            D_center = D_x[i]
            diffusion_term = (D_right * (C[n, i+1] - C[n, i]) - D_left * (C[n, i] - C[n, i-1])) / dx**2
            C[n+1, i] = C[n, i] + dt * (diffusion_term)
        C[n+1,:] -= C[n,:]*V_x*dt
        C[n+1, -1] = 0
        C[n+1, 0] = C[n+1, 1]
    C[-1,:] = [c*100 for c in C[-1,:]]
    plt.plot(x,C[-1,:],label =f'd = {v}')


# Plot the results
# x = [x*1e3 for x in x] #Convert to nm
C[0,:] = [c*100 for c in C[0,:]]
c_pre = [c*100 for c in c_pre]
err_pre = [e*100 for e in err_pre]
c_post= [c*100 for c in c_post]
x_post = [3*1e21*x/(n_atoms) for x in x_post]
err_post = [e*100 for e in err_post]

plt.errorbar(x_pre,c_pre,yerr= err_pre, ecolor='blue', capsize=3, label = 'As impanted', color ='blue')
plt.errorbar(x_post,c_post,yerr= err_post, ecolor='orange', capsize=3, label = 'Post-annealing', color ='orange')

plt.xlabel('Depth [nm]')
plt.ylabel(r'Concentration [at.\%]')
plt.grid(True)
plt.tight_layout()
plt.xlim([0,300])
plt.ylim([0,5.5])
plt.title(r'\textbf{Xe}')
plt.tight_layout()
plt.legend()

plt.show()
