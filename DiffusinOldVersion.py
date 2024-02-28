import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Parameters
L = 1 # Length of sample [micrometer]
T = 300  # Total time [seconds]
Nx = 100  # Number of spatial points
Nt = 300 # Number of time steps
D0 = 5.3e-3   # Diffusion coefficient inital value
Ea = 1.08 #Activation energy for diffusion in kJ/mole or eV depending on choise of k
dx = L / (Nx - 1) # Spatial step size
dt = T / Nt # Time step size
k = 8.6e-5 # boltzmann constant in eV/K
Temp = 2000 #Temperature in K
#k = 8.3e-3 # boltzmann constant [kJ/(mol*K)]
k = 8.6e-5 # boltzmann constant [ev/K]
Temp = 2000 #Temperature [K]
root = '/Users/niwi9751/Srim_Results/Fe_inZrO2_300keV.txt'

rho = 5.68 # density of target [g/cm^3]
m_a = 123.218 # molar mass of target in [g/mole]
Na = 6.022e22 # avogadros number 
n_atoms = rho*Na/m_a*10**-4 # atomic density of target [atoms/cm^3]
fluence = 1e17 # Input fluence of implantation [atoms/cm^2]
#/Users/niwi9751/Srim_Results/Fe_in_ZrO2.txt

def  read_columns(root):
    columns =  []
    with open(root,'r') as depthprofiles:
        lines = depthprofiles.readlines()
        for line in lines:
            column_data = line.split()
            for i,  column_data in enumerate(column_data):
                if len(columns) <= i:
                    columns.append([])
                columns[i].append(float(column_data.strip()))
    return columns
# Diffusion coefficient dependant on temperature
def D(D0, Ea, Temp):
    D = D0*np.exp(-Ea/(k*Temp))
    return D

# Gaussian
def gaussian(x,amp,mu,sigma):
    return amp*np.exp(-(x-mu)**2/(2*sigma**2))
# Fit data to gaussian
def fit_gauss(data, plot = False):
    n, bins, patches = plt.hist(data, bins = 60)
    x = np.linspace(min(data), max(data),60)
    x = x*0.0001 #rescale to micrometer
    y = n/0.11
    
    popt, pcov = curve_fit(gaussian,x,y)
    if plot:
        plt.plot(x,gaussian(x,*popt))
        plt.show()
    return popt

# Create spatial grid
x = np.linspace(0, L, Nx)

# Initialize solution matrix
C = np.zeros((Nt, Nx))

# Apply initial condition
data = read_columns(root)
popt = fit_gauss(data[1])
y = gaussian(x,*popt)
y = y*fluence/(3*n_atoms/100)
C[0, :] = y
print(C[0,:])


# # Time-stepping loop
# for n in range(0, Nt - 1):
#     # Apply direchlet boundary condition to RHS
#     C[n+1, -1] = 0
#     # Update interior points using forward difference in time and central difference in space
#     for i in range(1, Nx - 1):
#         C[n+1, i] = C[n, i] + D(D0, Ea, Temp) * dt / dx**2 * (C[n, i+1] - 2*C[n, i] + C[n, i-1])
#         # Apply Neumann boundary condition to LHS
#         C[n+1, 0] = C[n+1, 1]

# # Plot the results
# plt.figure(figsize=(8, 6))
# # for n in range(0, Nt, Nt//10):
# #     plt.plot(x, C[n, :], label=f"t={n*dt:.2f}")
# plt.plot(x,C[0,:], label = 'Initial distribution')
# plt.plot(x,C[-1,:], label = 'Final distribution')
# plt.title('Diffusion of Concentration')
# plt.xlabel('Position [micrometer]')
# plt.ylabel('Concentration [At. %]')
# plt.ylabel('Concentration [(atoms/cm^3)/(atoms/cm^2)]')
# plt.legend()
# plt.grid(True)
# plt.show()
# plt.show()
