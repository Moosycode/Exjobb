import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import datetime
import json
import os
import re

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

#Finding profiles
def Initialize_Profile(folder_path):
    if '.potku' in folder_path:                     #Check if potku is in the path name
        print('Folder is compatible, proceeding')   
        files = os.listdir(folder_path)
        possible_files = [folder_path] #Make list of folders in the request
    else:                                           
        print('Searching directory for potku files') #Search for potku files
        possible_files = [file for file in os.listdir(folder_path) if '.potku' in file] #make list of possible files
        
    if len(possible_files) == 0:                #if none are found, tell you
        print('Could not find .potku file, please try again')
    
    elif len(possible_files) > 1:
        print('Found compatible files: ')       
        [print(str(i + 1) + '.' + possible_files[i]) for i in range(len(possible_files))] #print list of possible files
        choise = possible_files[int(input('Chose the file you want (1 - ' + str(len(possible_files)) + ')'))-1]
        folder_path = folder_path + choise #chose one of them
    
    dict = {'Beam':{}, 'Samples':{}, 'Settings':{}} #Create dictionary
    
    try:
        beamdata = json.load(open(folder_path +'/Default/Default.measurement')) #Load info from folders
        beamprofile = json.load(open(folder_path+'/Default/Default.profile'))
        dict['Beam']['Ion'] = re.sub(r'\d+', '', beamdata['beam']['ion']) #Assign data
        dict['Beam']['Mass'] = re.sub(r'[a-zA-z]',  '' , beamdata['beam']['ion'])
        dict['Beam']['Energy'] = beamdata['beam']['energy']
        dict['Settings']['Num_step'] = beamprofile['depth_profiles']['number_of_depth_steps']
        dict['Settings']['Stop_step'] = beamprofile['depth_profiles']['depth_step_for_stopping']
        dict['Settings']['Out_step'] = beamprofile['depth_profiles']['depth_step_for_output']
    except:
        print('Corrupt Default.profile or Default.measurement file, please check them')
    
    for root, dirs, files in os.walk(folder_path):#Check all files in folder path
        for file in files:
            if file.endswith('.info'): #If infofile, we are in the right directory, keep looking here!!!
                currentsample = file.removesuffix('.info')
                dict['Samples'][currentsample] = {}  #Make it a sample
            if file.startswith('depth.') and 'total' not  in file: #Find the corresponding depht profiles, and save them.
                newroot = root.replace(os.path.sep,  '/') + '/' + file
                depthprofile = file.removeprefix('depth.')
                columns =  read_columns(newroot)
                if len(columns)== 7:
                    dict['Samples'][currentsample][depthprofile] = {'x': columns[0], 'C': columns[3],'NoNormC': columns[4],'N': columns[6]}
    return dict

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


#GENERAL DATA, INPUT URSELF
root = '/Users/niwi9751/Dropbox/Nils_files/Srim_Results/Zr300keV_in_UN_range.txt'
furnace_root = '/Users/niwi9751/Dropbox/Nils_files/furnanceDataTest.txt'
#potku_path = '/Users/niwi9751/potku/requests/20240304-KrXe-In-ZrO2.potku'
Times_in = [2]#Times in hours
Integrate = False

time, cool_curve = np.loadtxt(furnace_root,usecols=(0,1),unpack=True, dtype= 'U25') #Fix data format
time = [t.replace(',','.') for t in time]
cool_curve = [c.replace(',','.') for c in cool_curve]
time = np.array(time, dtype=float)
time = [excel_time_to_hours(t) for t in time] 
time = [int(t - time[0]) for t in time]
cool_curve = np.array(cool_curve, dtype=float)
Times = [T*3600 for T in Times_in]#Convert to seconds
Concentrations = []#Result list
MaxT_Times = []

#Dictionary with needed values of each element
elementdict = {
    'Fe_ZrO2':{'D0':2.26e-6,'Ea':2.3, 'rho':5.68, 'Ma': 123.218, 'N_at':3}, #Data from Springer
    'Zr_UN':{'D0':2.26e-6,'Ea':2.3, 'rho':14.05, 'Ma': 252.036, 'N_at':2}, 
    'Kr_ZrO2':{'D0':8.11e-7,'Ea':3.04, 'rho':5.68, 'Ma': 123.218, 'N_at':3},
    'B_Si':{'D0':0.76*1e8,'Ea':3.46, 'rho':2.33, 'Ma': 28.09, 'N_at':1}
}

for T in Times:
    # Parameters
    L = 3# Studied region [micrometer]
    studyL = 1 #Region of intrest, where SRIM starts/ends [micrometer] (HAS TO BE SAME LENGTH AS IN SRIM SIM)
    Nx = 100  # Number of spatial points per micrometer
    Nt = 100 # Number of time steps, can be anything really, code finds this for you, START LOW
    dx = L / (L*Nx - 1) # Spatial step size
    dt = T / Nt # Time step size
    k = 8.6e-5 # boltzmann constant [ev/K]
    Temp_fin = 1473.15 #Target emperature [K] 
    Temp = 300#Initial temperature [K]
    Na = 6.022e23 # avogadros number [atoms/mole]
    fluence = 1e16# Input fluence of implantation [atoms/cm^2]
    temps = []
    minutes = []
    
    element = 'Zr_UN'
    D0 = elementdict[element]['D0']*1e8# Diffusion coefficient inital value [um^2/s]
    Ea = elementdict[element]['Ea']# Activation energy for diffusion in kJ/mole or eV depending on choise of k
    rho = elementdict[element]['rho']# density of target [g/cm^3]
    m_a = elementdict[element]['Ma']# molar mass of target in [g/mole]
    at_in_mol = elementdict[element]['N_at']# Number of atoms in each molecule 
    
    n_atoms = rho*Na/m_a #atomic density of target [atoms/cm^3]
    heat = True
    cool = False
    
    #Check stability
    stability_cond = D(D0,Ea, Temp_fin)*dt/(dx**2)
    if stability_cond > 0.5:
        print('Not stable, increasing number of timesteps')
        while stability_cond > 0.4:
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
    
    height = np.loadtxt(root,usecols=1,unpack=True,encoding='cp437') #load height of bins
    height = height*fluence #convert into atoms/cm^3
    conc = height/(height + n_atoms) #Calculate concentration from number density of SRIM

    # Apply initial condition 
    C[0, :] = conc
    C = np.hstack((C, np.zeros((Nt,(L-studyL)*Nx)))) #Add zeros to desired length, SRIM length is only 1 micron usually
    
    
    cool_time = (Temp_fin-1273.15)/10 + time[-1] #Cooling time in minutes
    Current_min = 0
    Fin_min = 0
    # Time-stepping loop
    for n in range(0, Nt - 1):
        # Update interior points using forward difference in time and central difference in space
        second = T/Nt*n
        min = int(second/60)
        hour = min/60
        time_left = T/60-min
        if heat:
            if min != Current_min:  # Check if minute has changed
                minutes.append(min)
                temps.append(Temp)
                #print("A minute has passed at", second, "seconds.")
                Current_min = min  # Update current minute
                if Temp <= 1273:
                    Temp = Temp + 10
                    check = D(D0, Ea, Temp)*dt / dx**2
                    #print('Temp increasing! ')
                    #print(Temp)
                elif Temp < Temp_fin and Temp > 1273:
                    Temp = Temp + 5
                    Fin_min = Current_min
                    check = D(D0, Ea, Temp)*dt / dx**2
                    #print('Temp increasing! ')
                    #print(Temp)
                else:
                    heat = False     
        if (time_left - cool_time == 0):
            cool =  True
        if cool:
            if min != Current_min:  # Check if minute has changed
                minutes.append(min) 
                temps.append(Temp)
                Current_min = min  # Update current minute
                if Temp > 1273:
                    Temp = Temp - 10 #Linear decrease from Final temp to 1000 deg C
                else:
                    time_temp = int(cool_time - time_left) #Find
                    try:
                        index = time.index(time_temp) #Find index of time
                    except:
                        Temp = Temp #Temp is low enough here so fixing this is not worth it to fix this bug
                    Temp = int(cool_curve[index]+273) #Update temperature according to data
                
        Diff = D(D0, Ea, Temp) #calculate diffusion coefficient
        for i in range(1, L*Nx - 1): #Update interior points
            C[n+1, i] = C[n, i] + Diff * dt / dx**2 * (C[n, i+1] - 2*C[n, i] + C[n, i-1])
            # Apply Neumann boundary condition to boundaries
            C[n+1, 0] = C[n+1, 1]
            C[n+1, -1] = C[n+1,-2]

    #Tell time to reach maxtemp
    MaxT_time = min - Fin_min- cool_time
    print(f'Time to reach max temp: {Fin_min} minutes')
    print(f'Time at maxtemp: {MaxT_time} minutes')
    Concentrations.append(C[-1,:])
    MaxT_Times.append(int(MaxT_time/60))
    if Integrate:
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

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(x,C[0,:], label = 'Initial distribution')
i = 0
for C_ in Concentrations:        
    plt.plot(x,C_, label = f'Distribution after {Times_in[i]} hours. Time at maxtemp = {MaxT_Times[i]} hours')
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

plt.figure()
plt.title('Temperature change over time')
plt.xlabel('Minutes')
plt.ylabel('Temperature [K]')
plt.plot(minutes,temps)
plt.show()
