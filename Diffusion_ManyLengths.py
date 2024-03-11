import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

x_ = []
C_diffused = []
length = [1,4]
Temps = [25000, 30000, 36000, 43200]

for L in length:

    # Parameters
    #L = 4# Studied region [micrometer]
    studyL = 1 #Region of intrest [micrometer] (HAS TO BE SAME LENGTH AS IN SRIM SIM)
    T = 30000 # Total time [seconds]
    Nx = 100  # Number of spatial points per micrometer
    Nt = 30000 # Number of time steps
    D0 = 5.3e-3  # Diffusion coefficient inital value [cm^2/s]
    Ea = 1.08 #Activation energy for diffusion in kJ/mole or eV depending on choise of k
    dx = L / (L*Nx - 1) # Spatial step size
    dt = T / Nt # Time step size
    #k = 8.3e-3 # boltzmann constant [kJ/(mol*K)]
    k = 8.6e-5 # boltzmann constant [ev/K]
    Temp = 2000 #Temperature [K]
    root = '/Users/niwi9751/Srim_Results/Fe_inZrO2_300keV.txt'

    rho = 5.68 # density of target [g/cm^3]
    m_a = 123.218 # molar mass of target in [g/mole]
    Na = 6.022e22 # avogadros number [atoms/mole]
    n_atoms = rho*Na/m_a*10**-4 # atomic density of target [atoms/cm^3]
    fluence = 1e17 # Input fluence of implantation [atoms/cm^2]

    #/Users/niwi9751/Srim_Results/Fe_in_ZrO2.txt
    def  read_columns(root): #Use np.loadtxt instead.
        columns =  []
        with open(root,'r') as file:
            lines = file.readlines()
            for line in lines:
                column_data = line.split()
                for i,  column_data in enumerate(column_data):
                    if len(columns) <= i:
                        columns.append([])
                    columns[i].append(float(column_data.strip()))
        return columns

    def find_start(filename):
        try:
            with open(filename, 'r', encoding='utf-8') as file:
                for i, line in enumerate(file, 1):
                    if '--' in line:
                        # Check if the line contains a hyphen
                        return i  # Return the index (line number) if found
        except Exception as e:
            print("An error occurred while reading the file:", e)
            return None

    # Diffusion coefficient dependant on temperature
    def D(D0, Ea, Temp):
        D = D0*np.exp(-Ea/(k*Temp))
        return D

    # Gaussian
    def gaussian(x,amp,mu,sigma):
        return amp*np.exp(-(x-mu)**2/(2*sigma**2))


    def histog(data, length):
        L = length*10000 #Make sure unit is correct
        n, bins = np.histogram(data, bins = 100,range=(0,L))
        width = bins[1]-bins[0]
        x = np.linspace(min(data), L,100)
        x = x*0.0001 #rescale to micrometer, if file is in Å
        y = n/sum(n) #Normalize
        return width,y

    # Fit data to gaussian
    def fit_gauss(data, plot = False):
        n, bins, patches = plt.hist(data, bins = 60)
        x = np.linspace(min(data), max(data),60)
        x = x*0.0001 #rescale to micrometer
        y = n/0.11
        popt, pcov = curve_fit(gaussian,x,y)
        if plot:
            plt.figure('Initial data and gaussian fit')
            plt.plot(x,gaussian(x,*popt))
            plt.plot(x,y)
            plt.show()
        return popt

    def hist_integral(n, width):
        return sum(n*width)

    # Create spatial grid
    x = np.linspace(0, L, L*Nx)

    # Initialize solution matrix
    C = np.zeros((Nt, Nx))

    # Initialize initial condition
    data = read_columns(root)
    binwidth,y_hist = histog(data[1], studyL)
    y = y_hist*fluence/(3*n_atoms/100)

    # Apply initial condition 
    C[0, :] = y
    C = np.hstack((C, np.zeros((Nt,(L-studyL)*Nx)))) #Add zeros to desired length

    # Time-stepping loop
    for n in range(0, Nt - 1):
        # Update interior points using forward difference in time and central difference in space
        for i in range(1, L*Nx - 1):
            C[n+1, i] = C[n, i] + D(D0, Ea, Temp) * dt / dx**2 * (C[n, i+1] - 2*C[n, i] + C[n, i-1])
            # Apply Neumann boundary condition to boundaries
            C[n+1, 0] = C[n+1, 1]
            C[n+1, -1] = C[n+1,-2]




    x_.append(x)
    C_original = C[ 0,:]
    C_diffused.append(C[-1,:])
    
    I1 = hist_integral(C[0,:],binwidth)
    I2 = hist_integral(C[-1,:], binwidth)
    # I3 = hist_integral(C[-1,0:100], binwidth)
    # I4 = hist_integral(C[-1,100:200], binwidth)
    # I5 = hist_integral(C[-1,200:300], binwidth)
    print('-------------------')
    print(f'SRIM intagrel for length {L} = {I1}')
    print(f'diff intagrel for length {L} = {I2}')
    # print(I3)
    # print(I4)
    # print(I5)
    print(I1/I2)
    print('\n')




# Plot the results
plt.figure(figsize=(8, 6))
# for n in range(0, Nt, Nt//10):
#     plt.plot(x, C[n, :], label=f"t={n*dt:.2f}")
#plt.plot(x,C[0,:], label = 'Initial distribution')
#plt.plot(x,C[-1,:], label = 'Final distribution')
plt.plot(x,C_original, label = 'Initial distribution')
plt.plot(x_[0],C_diffused[0], label = f'Final distribution, length = {length[0]}')
plt.plot(x_[1],C_diffused[1], label = f'Final distribution, length = {length[1]}')
plt.title('Diffusion of Concentration')
plt.xlabel('Position [micrometer]')
plt.ylabel('Concentration [at. %]')
plt.legend()
plt.grid(True)
#plt.yscale('log')
plt.show()



