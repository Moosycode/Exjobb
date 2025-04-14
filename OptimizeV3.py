import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd

# ----------------- Physical Constants -----------------
k = 8.6e-5         # Boltzmann constant [eV/K]
Na = 6.022e23      # Avogadro's number [atoms/mol]
# Plot settings
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=16)

def D(D0, Ea, Temp):
    D = D0*np.exp(-Ea/(k*Temp))
    return D  

def hist_integral(n, width):

    return sum(count * width for count in n)

class SDTrimSP_output:

    def __init__( self, filename ):
        
        d_marker = 0  # data-type marker 
        h_marker = 0  # history-step marker
        d_data = []   # list to depth data
        s_data = []   # list to hold statistical data
        
        for line in reversed(list(open(filename))):
            
            if d_marker == 0 and line[0] == ' ':   # if stat. data hit while dm=0:
                d_marker +=1                       # mark for dep. data
                
            if d_marker == 1 and line[0] != ' ':   # if dep. data hit while dm=0:
                d_marker = 0                       # mark for stat. data
                h_marker +=1                       # mark for next hist. step
                
            if h_marker >= 1:                      # if 2nd hist. step is hit:
                break                              # end the file read --->>> ONLY LOOKS AT LAST HISTORY OUTPUT RIGHT NOW
            
            if d_marker == 0: d_data.append([float(i) for i in line.rstrip().split()])
            
            if d_marker == 1: s_data.append(line.rstrip().split())
            
        d_data = np.flipud(np.array(d_data)).T
        
        self.s = s_data
        
        self.nhist   =  float(s_data[-1][0])                    # number of histories
        self.srem    =  float(s_data[-1][1])/10                 # surface removal (nm)
        self.areald  = [1e16*float(i) for i in s_data[-4][:2]]  # areal densities (atoms/cm^2)
        self.npart   =  float(s_data[-5][0])                    # number of particles
        self.ibackp  = [float(i) for i in s_data[-6][:2]]       # backscattered particles
        self.ibackr  = [float(i) for i in s_data[-10][:2]]      # backscattered target
        
        self.depth   = d_data[0]/10                             # depth (nm)
        self.density = d_data[1]*1e24                           # atomic density (atoms/cm^3)
        self.afracs  = d_data[2:-2]
        self.Pdam    = d_data[-2]                               # relative probability of damage (max=1)
        
        self.Pdam_dn = self.Pdam/sum(self.Pdam)                 # relative probability of damage (area=1)
        self.binw_nm = self.depth[1]-self.depth[0]              # width of each depth-bin (nm)
        self.binw_cm = self.binw_nm * 1e-7                      # width of each depth-bin (cm)

        self.something    = d_data[-1]                          # NOT SURE WHAT THIS IS !!!!!!!!!

class data_tabel():

    def __init__( self, filename, delimiter='\t' ):
        with open(filename, 'r') as file:
            header_line = file.readline().strip()    # Read the first line and strip any trailing newline characters
        self.headers = header_line.split(delimiter)  # Split the header line into a list
        self.ncols = len(self.headers)
        self.data = np.loadtxt( filename, delimiter=delimiter, skiprows=1, usecols=np.arange(self.ncols)).T
        
    def get_data( self, header ):
        for i in range(len(self.headers)):
            if self.headers[i] == header:
                break
        return self.data[i]

def get_trim_res( path ):
    filename = 'E0_31_target.dat'
    d = SDTrimSP_output( path+filename )
    print('-------------------------------------------')
    print()
    print(f'number of histories: {int(d.nhist)} ')
    print(f'number of particles: {int(d.npart)} ')
    print(f'surface removal:     {np.round(d.srem,1)} nm')
    print(f'areal density of Kr from SDTrimSP: {d.areald[0]:.2e} atoms/cm^2')
    print(f'total area of Pdam = {sum(d.Pdam)}')
    print()
    return d


# ---------------- Element Database --------------------
#Inital params
elementdict = {
    'Fe_ZrO2':{'D0':2.26e-6,'Ea':2.3, 'rho':6.08, 'Ma': 123.218}, #Data from Springer
    'Kr_ZrO2':{'D0':8.11e-7,'Ea':2.53, 'rho':6.08, 'Ma': 123.218},
    'Xe_ZrO2':{'D0':1.83e-6,'Ea':2.91, 'rho':6.08, 'Ma': 123.218},
    'Zr_UN':{'D0':6.9e-7,'Ea':2.7, 'rho':14.05, 'Ma': 252.036}, 
    'Kr_UN':{'D0':8.11e-7,'Ea':2.53, 'rho':14.05, 'Ma': 252.036},
    'Xe_UN':{'D0':1.83e-6,'Ea':2.91, 'rho':14.05, 'Ma': 252.036},
    'Zr_UO2':{'D0':6.9e-7,'Ea':2.7, 'rho':10.6, 'Ma': 270.02},
    'Kr_UO2':{'D0':8.11e-7,'Ea':2.53, 'rho':10.6, 'Ma': 270.02}
}
# ---------------- Simulation Setup --------------------

# File paths

pre_data_path = 'Data/ToFERDA/Data/pre-anneal_Xe-imp.csv'
post_data_path = 'Data/ToFERDA/Data/post-anneal_Xe-imp.csv'
# pre_data_path = 'Data/ToFERDA/Data/pre-anneal_Kr-imp_corrected.csv'
# post_data_path = 'Data/ToFERDA/Data/post-anneal_Kr-imp_corrected.csv'


# Experiment conditions
sample = 'Xe'
target = 'ZrO2'
element = sample + '_' + target
T = 9 * 3600                                     # in hours
Temp = 1473.15                                   # in Kelvin
fluence = 1e17                                   # at/cm2
Integrate=True
rho = elementdict[element]['rho']               #density of target [g/cm^3]
m_a = elementdict[element]['Ma']                # atomic mass of target in [g/mole]
n_atoms = rho*Na/m_a                            #atomic density of target [atoms/cm^3]
scores =[]
Concentrations = []
d = get_trim_res('Data/SDTrim/' + sample + '/')
interp_func = interp1d(d.depth, d.Pdam, kind='linear', fill_value="extrapolate")
df_pre = pd.read_csv(pre_data_path)
df_pre = df_pre[df_pre['depth'] <= 3000]
df_post = pd.read_csv(post_data_path)
df_post = df_post[df_post['depth'] <= 3000]


x_pre = df_pre['depth'].values
x_post = df_post['depth'].values

#Constants--------------------------------------------------------
rho = elementdict[element]['rho']# density of target [g/cm^3]
m_a = elementdict[element]['Ma']# atomic mass of target in [g/mole]
n_atoms = rho*Na/m_a #atomic density of target [atoms/cm^3]
Nx = len(x_pre)  # Number of spatial points per micrometer
Nt = int(T/60) # Number of time steps, can be anything really, code finds this for you but do not start lower than this.
x_pre = [3*1e21*x/(n_atoms) for x in x_pre] #Convert to micrometer
Extendby = 2                                                                                            
L = int(x_pre[-1])*Extendby
studyL = int(x_pre[-1])
dx = L / (L*Nx - 1) # Spatial step size
dt = T / Nt # Time step size
x_post = [3*1e21*x/(n_atoms) for x in x_post]
#-----------------------------------------------------------------


def optifunc(vars,plot = False):
    c_pre = df_pre[sample].values
    c_post = df_post[sample].values
    D0, Ea, a , d= vars
    # Create spatial grid
    x = np.linspace(0, L, Nx*Extendby)
    # Initialize solution matrix
    C = np.zeros((Nt, Nx))
    # Apply initial condition 
    C[0, :] = c_pre
    C = np.hstack((C, np.zeros((Nt,(Extendby-1)*Nx)))) #Add zeros to desired length
    rebinned = interp_func(x)
    rebinned = [d*x/sum(rebinned) + 1  for x in rebinned]
    D_x = np.array([D(D0, Ea, Temp)*rebinned[i] for i in range(Nx*Extendby)])

    for n in range(Nt - 1):                         # Main loop
        Cn = C[n]               
        diffusion_term = (D_x[2:] * (Cn[2:] - Cn[1:-1]) - D_x[:-2] * (Cn[1:-1] - Cn[:-2])) / dx**2
        C[n+1, 1:-1] = Cn[1:-1] + dt * diffusion_term
        C[n+1] -= C[n] * a * dt
        C[n+1, -1] = 0
        C[n+1, 0] = C[n+1, 1]
    Concentrations.append(C[-1,:])

    score = sum([(c-c_)**2 for c,c_ in zip(C[-1,:],c_post)])
    scores.append(score)
    
    if plot:
        plt.figure(figsize =(8,6))
        c_pre = [c*100 for c in c_pre]
        Cres = [c*100 for c in C[-1,:]]
        c_post = [c*100 for c in c_post]
        plt.plot(x_pre,c_pre, label = 'Pre annealing')
        plt.plot(x,Cres, label = 'Optimized')
        plt.plot(x_post,c_post, label = 'Post annealing')
        plt.xlabel('Depth [nm]')
        plt.ylabel(r'Concentration [at. \%]')
        plt.grid(True)
        plt.xlim((0,300))
        plt.tight_layout()
        plt.legend()
        plt.show()
    return score

initial_guess = [elementdict[element]['D0']*1e8, elementdict[element]['Ea'], 5e-5, 500]
bounds = [(1,None),(1,4),(1e-5,1e-4), (1e-15,None)]

result = minimize(optifunc, initial_guess, method="Nelder-Mead", bounds=bounds, tol=1e-2)
optivals = []
for val in result['x']:
    optivals.append(val)
print(optivals)
optifunc(result['x'],True)
