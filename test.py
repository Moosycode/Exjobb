import numpy as np

L = 4# Studied region [micrometer]
studyL = 1 #Region of intrest [micrometer] (HAS TO BE SAME LENGTH AS IN SRIM SIM)
T = 30000 # Total time [seconds]
Nx = 100  # Number of spatial points per micrometer
Nt = 6000 # Number of time steps
D0 = 5.3e-3  # Diffusion coefficient inital value [cm^2/s]
Ea = 1.08 #Activation energy for diffusion in kJ/mole or eV depending on choise of k
dx = L / (L*Nx - 1) # Spatial step size
dt = T / Nt # Time step size
k = 8.6e-5 # boltzmann constant [ev/K]
Temp = 2000 #Temperature [K]


def D(D0, Ea, Temp):
    D = D0*np.exp(-Ea/(k*Temp))
    return D

val = D(D0,Ea, Temp)*dt/(dx**2)
print(val)