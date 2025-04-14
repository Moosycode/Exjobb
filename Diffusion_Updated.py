import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd

# ----------------- Physical Constants -----------------
k = 8.6e-5         # Boltzmann constant [eV/K]
Na = 6.022e23      # Avogadro's number [atoms/mol]

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
#Optimzed params
opt_dict = {
    'Kr': [84.98484603683865, 2.778815397425465, 6.001900117458542e-05, 128.6670625458783],
    'Xe': [34.719190870006244, 3.1622358926720224, 5.724240280220646e-05, 367.1035590595884]
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
T = 9 * 3600                                     # in hours
Temp = 1473.15     # in Kelvin
fluence = 1e17
element = sample + '_' + target
Integrate=True
rho = elementdict[element]['rho']               #density of target [g/cm^3]
m_a = elementdict[element]['Ma']                # atomic mass of target in [g/mole]
n_atoms = rho*Na/m_a                            #atomic density of target [atoms/cm^3]
Concentrations = []
d = get_trim_res('Data/SDTrim/' + sample + '/')
interp_func = interp1d(d.depth, d.Pdam, kind='linear', fill_value="extrapolate")
# Plot settings
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=16)

#Select sample
optivals = opt_dict[sample]
param_sets = [
    {'label': 'Initial guess', 'D0': elementdict[element]['D0']*1e8, 'Ea': elementdict[element]['Ea'], 'a': 1e-5,       'd': 1},
    {'label': 'Optimized',     'D0': optivals[0],                    'Ea': optivals[1],                'a':optivals[2], 'd':optivals[3]}
]

for params in param_sets:
    print(f"\n--- Running diffusion with: {params['label']} ---")

    D0, Ea, a, d = params['D0'], params['Ea'], params['a'], params['d']
    
    # Parameters------------------------------------------------------
    df_pre = pd.read_csv(pre_data_path)
    df_pre = df_pre[df_pre['depth'] <= 3000]
    df_post = pd.read_csv(post_data_path)
    df_post = df_post[df_post['depth'] <= 3000]

    c_pre = df_pre[sample].values
    c_post = df_post[sample].values
    x_pre = df_pre['depth'].values
    x_post = df_post['depth'].values

    binwidth = (x_pre[1]-x_pre[0])*1e15 #at /cm2
    if Integrate:
        I_inital = hist_integral(c_pre,binwidth)
        I_diffused = hist_integral(c_post, binwidth)
        ratio = I_inital/I_diffused
        print(f'Total integral of inital concentration: ')
        print(I_inital)
        print(f'Total integral of diffused concentration: ')
        print(I_diffused)
        print(f'Ratio between total integrals: ')
        print(ratio)

    # --- Physical Parameters ---
    Nx = len(x_pre)                                 # Number of spatial points per micrometer
    Nt = int(T/60)*2                                # Number of time steps, can be anything really, code finds this for you but do not start lower than this.
    x_pre = [3*1e21*x/(n_atoms) for x in x_pre]     # Convert to nm
    x_post = [3*1e21*x/(n_atoms) for x in x_post]   # Convert to nm
    
    Extendby = 2                                    # Extend domain by X many times
    L = int(x_pre[-1])*Extendby                     # Length of domain
    dx = L / (L*Nx - 1)                             # Spatial step size
    dt = T / Nt                                     # Time step size
    x = np.linspace(0, L, Nx*Extendby)              # Create spatial grid
    C = np.zeros((Nt, Nx * Extendby))               # Initialize solution matrix
    C[0, :Nx] = c_pre                               # Initial condition

    Pdam = interp_func(x)                           # Damage factor
    Pdam = [d*x/sum(Pdam) + 1 for x in Pdam]        # Damage factor, normalized + 1

    D_x = np.array([D(D0, Ea, Temp)*Pdam[i] for i in range(Nx*Extendby)])  # Diffusion coefficient

    for n in range(Nt - 1):                         # Main loop
        Cn = C[n]               
        diffusion_term = (D_x[2:] * (Cn[2:] - Cn[1:-1]) - D_x[:-2] * (Cn[1:-1] - Cn[:-2])) / dx**2
        C[n+1, 1:-1] = Cn[1:-1] + dt * diffusion_term
        C[n+1] -= C[n] * a * dt
        C[n+1, -1] = 0
        C[n+1, 0] = C[n+1, 1]

    Concentrations.append(C[-1])

# Integrate post distribution
bin_w = (x[1]-x[0])*1e-7
postint = hist_integral(Concentrations[-1][:21],bin_w)
print(f'Total integral of diffused concentration as predicted by model:')
print(f'{postint*n_atoms*3:.2e}')

# Plot the results
plt.figure(figsize=(8, 6))

scale = 100
plt.errorbar(x_pre, c_pre * scale, yerr=df_pre[f'{sample}_sd'] * scale,
             label='As implanted', capsize=3, color='blue')
plt.errorbar(x_post, c_post * scale, yerr=df_post[f'{sample}_sd'] * scale,
             label='Post-anneal', capsize=3, color='orange')

for i, C_ in enumerate(Concentrations):
    plt.plot(x, C_ * scale, linestyle='--',
             label=f'{i+1}$^{{nd}}$ model' if i else f'{i+1}$^{{st}}$ model', color='green' if i else 'red')

handles, labels = plt.gca().get_legend_handles_labels()
plot_order = [2, 3, 0, 1]
plt.xlabel('Depth [nm]')
plt.ylabel('Concentration [at.\%]')
plt.title(r'\textbf{' + sample + '}')
plt.grid(True)
plt.xlim([0, 300])
plt.ylim([0, 6])
plt.legend([handles[i] for i in plot_order], [labels[i] for i in plot_order] )
plt.tight_layout()



plt.figure(figsize=(8,6))
plt.plot(x, D_x)
plt.title('Diffusion constant')
plt.xlabel('Depth [nm]')
plt.ylabel(r'D(x) [$\mu$ m$^2$/s]')
plt.grid(True)
plt.tight_layout()

plt.show()
