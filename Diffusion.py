import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Parameters
L = 3 # Length in micrometer
T = 500  # Total time in seconds
Nx = 100  # Number of spatial points
Nt = 500 # Number of time steps
D0 = 2e-4  # Diffusion coefficient inital value
Ea = 2.5 #Activation energy for diffusion in eV
dx = L / (Nx - 1) # Spatial step size
dt = T / Nt # Time step size
k = 8.6e-5 # boltzmann constant in eV/K
Temp = 2000 #Temperature in K
root = '/Users/niwi9751/Srim_Results/Fe_in_ZrO2.txt'

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
    return D0

# Gaussian
def gaussian(x,amp,mu,sigma):
    return amp*np.exp(-(x-mu)**2/(2*sigma**2))

# Fit data to gaussian
def fit_gauss(data, plot = False):
    n, bins, patches = plt.hist(data, bins = 60)
    x = np.linspace(min(data), max(data),60)
    x = x*0.0001 #rescale to micrometer
    y = n
    y = 100*y/len(data) #rescale to atomic %
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
C[0, :] = gaussian(x,*popt)

# Time-stepping loop
for n in range(0, Nt - 1):
    # Apply direchlet boundary condition to RHS
    C[n+1, -1] = 0
    # Update interior points using forward difference in time and central difference in space
    for i in range(1, Nx - 1):
        C[n+1, i] = C[n, i] + D(D0, Ea, Temp) * dt / dx**2 * (C[n, i+1] - 2*C[n, i] + C[n, i-1])
        # Apply Neumann boundary condition to LHS
        C[n+1, 0] = C[n+1, 1]

# Plot the results
plt.figure(figsize=(8, 6))
for n in range(0, Nt, Nt//10):
    plt.plot(x, C[n, :], label=f"t={n*dt:.2f}")

plt.title('Diffusion of Concentration')
plt.xlabel('Position [micrometer]')
plt.ylabel('Concentration [At. %]')
plt.legend()
plt.grid(True)
plt.show()