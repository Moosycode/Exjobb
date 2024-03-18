import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import datetime

Times_in = [18, 24]#Times in hours
Times = [T*3600 for T in Times_in]#Convert to seconds
Concentrations = []#Result list
root = '/Users/niwi9751/Srim_Results/Fe_inZrO2_300keV.txt'
furnace_root = '/Users/niwi9751/furnanceData.txt'
# D0 = 2.69e-4 #For Zr in U
# Ea = 0.521 # For Zr in U   
for T in Times:
    # Parameters
    L = 4# Studied region [micrometer]
    studyL = 1 #Region of intrest, where SRIM starts/ends [micrometer] (HAS TO BE SAME LENGTH AS IN SRIM SIM)
    Nx = 100  # Number of spatial points per micrometer
    Nt = 1 # Number of time steps, can be anything really, code finds this for you, START LOW
    D0 = 5.3e-3  # Diffusion coefficient inital value [cm^2/s]
    Ea = 1.08# Activation energy for diffusion in kJ/mole or eV depending on choise of k
    dx = L / (L*Nx - 1) # Spatial step size
    dt = T / Nt # Time step size
    k = 8.6e-5 # boltzmann constant [ev/K]
    Temp_fin = 1873 #Temperature [K]
    Temp_in = 300 #Initial temperature [K]
    Temp = 300
    rho = 5.68 # density of target [g/cm^3]
    m_a = 123.218# molar mass of target in [g/mole]
    Na = 6.022e22 # avogadros number [atoms/mole]
    n_atoms = rho*Na/m_a #atomic density of target [atoms/cm^3]
    at_in_mol = 3 # Number of atoms in each molecule 
    fluence = 1e17 # Input fluence of implantation [atoms/cm^2]
    heat = True
    cool = False
    cool_time = 563 #Cooling time in minutes
    
    def find_start(filename):
        try:
            with open(filename, 'r', encoding='cp437') as file:
                for i, line in enumerate(file, 1):
                    if '-------  ----------- ----------- -----------' in line:
                        # Check if the line contains a hyphen
                        return i  # Return the index (line number) if found
        except Exception as e:
            print("An error occurred while reading the file:", e)
            return 0

    # Diffusion coefficient dependant on temperature
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

    # Gaussian
    def gaussian(x,amp,mu,sigma):
        return amp*np.exp(-(x-mu)**2/(2*sigma**2))

    def histog(data, length):
        L = length*10000 #Make sure unit is correct
        n, bins = np.histogram(data, bins = 100,range=(0,L))
        width = bins[1]-bins[0]
        # x = np.linspace(min(data), L,100)
        # x = x*0.0001 #rescale to micrometer, if file is in Å
        y = n/sum(n) #Normalize
        return width,y

    # Fit data to gaussian, NOT USED IN THIS VERSION
    def fit_gauss(data, plot = False):
        n, bins, patches = plt.hist(data, bins = 60)
        x = np.linspace(min(data), max(data),60)
        x = x*0.0001 #rescale to micrometer
        y = n/0.11 #sketchy normalization, not correct...
        popt, pcov = curve_fit(gaussian,x,y)
        if plot:
            plt.figure('Initial data and gaussian fit')
            plt.plot(x,gaussian(x,*popt))
            plt.plot(x,y)
            plt.show()
        return popt

    def hist_integral(n, width):
        return sum(n*width)#Definition of integrals :))

    #Check stability
    stability_cond = D(D0,Ea, Temp_fin)*dt/(dx**2)
    if stability_cond > 0.5:
        print('Not stable, increasing number of timesteps')
        while stability_cond > 0.33:
            Nt = Nt+1
            dt = T/Nt
            stability_cond = D(D0,Ea, Temp_fin)*dt/(dx**2)
        print(f'Continuing with new number of timesteps: {Nt}')

    # Create spatial grid
    x = np.linspace(0, L, L*Nx)

    # Initialize solution matrix
    C = np.zeros((Nt, Nx))

    # # Initialize initial condition
    start = find_start(root)
    if start == None:
        start = 0
    
    data = np.loadtxt(root,skiprows=start,usecols=1,unpack=True,encoding='cp437')
    binwidth,y_hist = histog(data, studyL)
    y = y_hist*fluence/(at_in_mol*n_atoms*1e-6) # Atoms/cm^2 DIVIDED BY atoms/cm^3 * thickness (USE THICKNESS OF BIN)

    # Apply initial condition 
    C[0, :] = y
    C = np.hstack((C, np.zeros((Nt,(L-studyL)*Nx)))) #Add zeros to desired length
    Current_min = 0
    Fin_min = 0
    
    time, cool_curve = np.loadtxt(furnace_root,usecols=(0,1),unpack=True, dtype= 'U25')
    time = [t.replace(',','.') for t in time]
    cool_curve = [c.replace(',','.') for c in cool_curve]
    time = np.array(time, dtype=float)
    time = [excel_time_to_hours(t) for t in time] 
    time = [int(t - time[0]) for t in time]
    cool_curve = np.array(cool_curve, dtype=float)
    
    # Time-stepping loop
    for n in range(0, Nt - 1):
        # Update interior points using forward difference in time and central difference in space
        second = T/Nt*n
        min = int(second/60)
        hour = min/60
        time_left = T/60-min
        if heat:
            if min != Current_min:  # Check if minute has changed
                #print("A minute has passed at", second, "seconds.")
                Current_min = min  # Update current minute
                if Temp <= 1273:
                    Temp = Temp + 10
                    #print('Temp increasing! ')
                    # print(Temp)
                elif Temp < Temp_fin and Temp > 1273:
                    Temp = Temp + 5
                    Fin_min = Current_min
                    #print('Temp increasing! ')
                    # print(Temp)
                else:
                    heat = False
                
        if (time_left - cool_time == 0):
            cool = True
        if cool:
            if min != Current_min:  # Check if minute has changed
                #print("A minute has passed at", second, "seconds.")
                Current_min = min  # Update current minute
                time_temp =int(cool_time - time_left)
                index = time.index(time_temp)
                Temp = cool_curve[index]+273
        print(Temp)
        for i in range(1, L*Nx - 1):
            C[n+1, i] = C[n, i] + D(D0, Ea, Temp) * dt / dx**2 * (C[n, i+1] - 2*C[n, i] + C[n, i-1])
            # Apply Neumann boundary condition to boundaries
            C[n+1, 0] = C[n+1, 1]
            C[n+1, -1] = C[n+1,-2]

    #Tell time
    print(f'Time to reach max temp: {Fin_min} minutes')
    print(f'Time at maxtemp: {min -Fin_min} minutes')

    #Integrate over intresting areas in steps of 1 micrometer
    I_interval = []
    for i in range(L):
        Integral = hist_integral(C[-1,i*100:(i+1)*100],binwidth)
        I_interval.append(Integral)
        print(f'Integral between {i} and {i+1} micrometer: {Integral}')
        print('----------------------')
        
    #Integrate over total length
    I_inital = hist_integral(C[0,:],binwidth)
    I_diffused = hist_integral(C[-1,:], binwidth)
    ratio = I_inital/I_diffused
    print(f'Total integral of inital concentration at length {L} micrometer: ')
    print(I_inital)
    print(f'Total integral of diffused concentration at length {L} micrometer: ')
    print(I_diffused)
    print(f'Ratio between total integrals: ')
    print(ratio)
    print()
    Int_x = np.linspace(1,L,L)
    Concentrations.append(C[-1,:])

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(x,C[0,:], label = 'Initial distribution')
i = 0
for C in Concentrations:        
    plt.plot(x,C, label = f'Distribution after {Times_in[i]} hours')
    i = i + 1
plt.title('Diffusion of Concentration')
plt.xlabel('Position [micrometer]')
plt.ylabel('Concentration [at. %]')
plt.legend()
plt.grid(True)
plt.show()




