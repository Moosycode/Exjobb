import matplotlib.pyplot as plt
import numpy as np
import datetime

def D(D0, Ea, Temp):
    D = D0*np.exp(-Ea/(k*Temp))
    return D

def excel_time_to_hours(excel_time):
    # Extracting the integer part as days
    days = int(excel_time)
    # Extracting the fractional part as seconds
    seconds = (excel_time - days) * 24 * 3600
    # Converting to timedelta
    time_delta = datetime.timedelta(seconds=seconds)
    # Extracting hours, minutes, and seconds
    hours = int(time_delta.total_seconds() // 3600)
    minutes = int((time_delta.total_seconds() % 3600) // 60)
    seconds = int(time_delta.total_seconds() % 60)
    return hours*60 + minutes + seconds/60
    
# Read the CSV file
filename = '/Users/nilsw/Dropbox/Nils_files/SRIM_Results/test4.txt'  # Replace 'your_file.csv' with the path to your file
with open(filename, 'r') as file:
     lines = file.readlines()

# Extract data from the first two columns
column1 = []
column2 = []
for line in lines:
    parts = line.strip().split(',')
    if len(parts) >= 2:  # Ensure there are at least two columns of data
        column1.append(float(parts[0]))
        column2.append(float(parts[1]))
        
        
Times_in = [2,10]#Times in hours
Times = [T*3600 for T in Times_in]#Convert to seconds
Concentrations = []#Result list
#root = '/Users/niwi9751/Srim_Results/Fe_inZrO2_300keV.txt'
root = '/Users/nilsw/Dropbox/Nils_Files/Srim_Results/B300keV_in_Si.txt'
furnace_root = '/Users/niwi9751/furnanceData.txt'
#potku_path = '/Users/niwi9751/potku/requests/20240304-KrXe-In-ZrO2.potku'
elementdict = {
    'Fe_ZrO2':{'D0':1.13e-7*1e8,'Ea':2.7, 'rho':5.68, 'Ma': 123.218, 'N_at':3},
    'Zr_UN':{'D0':2.69e-4,'Ea':0.521, 'rho':13.9, 'Ma': 518.078, 'N_at':2},
    'Kr_ZrO2':{'D0':8.11e-7,'Ea':3.04, 'rho':5.68, 'Ma': 123.218, 'N_at':3},
    'B_Si':{'D0':0.22*1e8,'Ea':3.46, 'rho':2.33, 'Ma': 28.09, 'N_at':1}
}

for T in Times:
    # Parameters
    L = 2# Studied region [micrometer]
    studyL = 2 #Region of intrest, where SRIM starts/ends [micrometer] (HAS TO BE SAME LENGTH AS IN SRIM SIM)
    Nx = 100  # Number of spatial points per micrometer
    Nt = 1 # Number of time steps, can be anything really, code finds this for you, START LOW
    dx = L / (L*Nx - 1) # Spatial step size
    dt = T / Nt # Time step size
    k = 8.6e-5 # boltzmann constant [ev/K]
    Temp_fin = 1273.15 #Target emperature [K] 
    Temp = 1273.15#Initial temperature [K]
    Na = 6.022e23 # avogadros number [atoms/mole]
    fluence = 1e17# Input fluence of implantation [atoms/cm^2]
    temps = []
    cooltemps = []
    
    element = 'B_Si'
    D0 = elementdict[element]['D0']# Diffusion coefficient inital value [cm^2/s]
    Ea = elementdict[element]['Ea']# Activation energy for diffusion in kJ/mole or eV depending on choise of k
    rho = elementdict[element]['rho']# density of target [g/cm^3]
    m_a = elementdict[element]['Ma']# molar mass of target in [g/mole]
    at_in_mol = elementdict[element]['N_at']# Number of atoms in each molecule 
    
    n_atoms = rho*Na/m_a #atomic density of target [atoms/cm^3]
    heat = False
    cool = False
    cool_time = 563 #Cooling time in minutes
    
    #Check stability
    stability_cond = D(D0,Ea, Temp_fin)*dt/(dx**2)
    if stability_cond > 0.5:
        print('Not stable, increasing number of timesteps')
        while stability_cond > 0.5:
            Nt = Nt+1
            dt = T/Nt
            stability_cond = D(D0,Ea, Temp_fin)*dt/(dx**2)
        print(f'Continuing with new number of timesteps: {Nt}')

    # Create spatial grid
    x = np.linspace(0, L, Nx)

    # Initialize solution matrix
    C = np.zeros((Nt, Nx))

    # # Initialize initial condition
    C[0,:] = column2
    #/(at_in_mol*n_atoms*1e-6) # Atoms/cm^2 DIVIDED BY atoms/cm^3 * thickness (USE THICKNESS OF BIN)

    # Apply initial condition 
    
    # Time-stepping loop
    for n in range(0, Nt - 1):
        # Update interior points using forward difference in time and central difference in space
        Diff = D(D0, Ea, Temp)
        for i in range(1, Nx - 1):
            C[n+1, i] = C[n, i] + Diff * dt / dx**2 * (C[n, i+1] - 2*C[n, i] + C[n, i-1])
            # Apply Neumann boundary condition to boundaries
            C[n+1, 0] = C[n+1, 1]
            C[n+1, -1] = C[n+1,-2]
          
    Concentrations.append(C[-1,:])
            
# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(x,C[0,:], label = 'Nils code, 0 hours')
i = 0
for C_ in Concentrations:        
    plt.plot(x,C_, label = f'Nils code, {Times_in[i]} hours')
    i = i + 1

#Measured plot
# potku_data = Initialize_Profile(potku_path)
# x_pot = potku_data['Samples']['Fe-Imp']['Fe']['x']
# x_pot = [at_in_mol*1e18*x/(n_atoms) for x in x_pot] #Convert to micrometer
# c_pot = potku_data['Samples']['Fe-Imp']['Fe']['C']
# plt.plot(x_pot,c_pot, label = 'Measured distribution')
plt.title('Diffusion of Concentration')
plt.xlabel('Position [micrometer]')
plt.ylabel('Concentration [at. fraction]')
plt.grid(True)
plt.legend()
# Read the CSV file
filenames = ['/Users/nilsw/Dropbox/Nils_Files/SRIM_Results/test.txt','/Users/nilsw/Dropbox/Nils_Files/SRIM_Results/test2.txt','/Users/nilsw/Dropbox/Nils_Files/SRIM_Results/test3.txt']  # Replace 'your_file.csv' with the path to your file
i = 0
times = [0,2,10]
for filename in filenames:
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Extract data from the first two columns
    column1 = []
    column2 = []
    for line in lines:
        parts = line.strip().split(',')
        if len(parts) >= 2:  # Ensure there are at least two columns of data
            column1.append(float(parts[0]))
            column2.append(float(parts[1]))

    # Plot the data
    plt.plot(column1, column2, label = f'Cleanroom code, {times[i]} hours')
    i = i +1


plt.grid(True)
plt.legend()
plt.show()
