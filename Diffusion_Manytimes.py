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

def rebin(data,x_res):
    data = [(data[1::2][i] + data[::2][i])/2 for i in range(len(data[1::2]))]
    x_res = x_res[1::2]
    return data, x_res

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
                if '-----------  ----------  ------------' in line:
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

def histog(data, length):
    L = length*10000 #Make sure unit is correct
    n, bins = np.histogram(data, bins = 100,range=(0,L))
    width = bins[1]-bins[0]
    y = n/sum(n) #Normalize
    return width,y

def hist_integral(n, width):
    n = [item*width for item in n]
    return sum(n)#Definition of integrals :))


#Constants----------------------------
k = 8.6e-5 # boltzmann constant [ev/K]
Na = 6.022e23 # avogadros number [atoms/mole]
#-------------------------------------

#GENERAL FILEPATHS-----------------------------------------------------------------
# root = '/Users/niwi9751/Dropbox/Nils_files/Srim_Results/Kr300keV_in_ZrO2_range.txt' 
root = '/Users/niwi9751/Dropbox/Nils_files/Srim_Results/Zr330keV_in_UN_Range.txt'
furnace_root = '/Users/niwi9751/Dropbox/Nils_files/furnanceDataTest.txt'
# # potku_path = '/Users/niwi9751/potku/requests/20240410-Zr-in-UN.potku'
# potku_path = '/Users/niwi9751/potku/requests/20240506-UNUO2Samples.potku'
potku_path2 = '/Users/niwi9751/potku/requests/20240304-KrXe-In-ZrO2.potku'
potku_path = '/Users/niwi9751/potku/requests/20240521-PostAnnealZrO2.potku'
#---------------------------------------------------------------------------------

#Global Parameters-----------------------------------------------------------------------
Times_in = [5,10,25,50]#Times in hours
L = 3# Studied region [micrometer]
studyL = 1 #Region of intrest, where SRIM starts/ends [micrometer] (HAS TO BE SAME LENGTH AS IN SRIM SIM)
Temp_fin = 1673.15 #Target emperature [K] 
fluence = 9.7e15# Input fluence of implantation [atoms/cm^2]
Integrate = True
Concentrations = []#Result list
MaxT_Times = []#Honestly do not remember
Mins = []#All minutes globally
Temperatures = []#All temperatures globally
#---------------------------------------------------------------------------------


#Stupid stuff since i am not that good at coding-------------------------------------
time, cool_curve = np.loadtxt(furnace_root,usecols=(0,1),unpack=True, dtype= 'U25') #Fix data format
time = [t.replace(',','.') for t in time]
cool_curve = [c.replace(',','.') for c in cool_curve]
time = np.array(time, dtype=float)
time = [excel_time_to_hours(t) for t in time] 
time = [int(t - time[0]) for t in time]
cool_curve = np.array(cool_curve, dtype=float)
Times = [T*3600 for T in Times_in]#Convert to seconds
#---------------------------------------------------------------------------------

#Dictionary with needed values of each element
elementdict = {
    'Fe_ZrO2':{'D0':2.26e-6,'Ea':2.3, 'rho':6.025, 'Ma': 123.218}, #Data from Springer
    'Kr_ZrO2':{'D0':8.11e-7,'Ea':2.53, 'rho':6.025, 'Ma': 123.218},
    'Xe_ZrO2':{'D0':1.83e-6,'Ea':2.91, 'rho':6.025, 'Ma': 123.218},
    'Zr_UN':{'D0':6.9e-7,'Ea':2.7, 'rho':14.05, 'Ma': 252.036}, 
    'Kr_UN':{'D0':2e-4,'Ea':4.66, 'rho':14.05, 'Ma': 252.036},
    'Xe_UN':{'D0':2e-4,'Ea':4.71, 'rho':14.05, 'Ma': 252.036},
    'Zr_UO2':{'D0':2e-4,'Ea':4.71, 'rho':10.6, 'Ma': 270.02},
}

for T in Times:
    # Parameters------------------------------------------------------
    element = 'Zr_UN'
    Temp = 1473.15#Initial temperature [K]
    #-----------------------------------------------------------------
    
    #Constants--------------------------------------------------------
    Nx = 100  # Number of spatial points per micrometer
    Nt = int(T/60) # Number of time steps, can be anything really, code finds this for you but do not start lower than this.
    dx = L / (L*Nx - 1) # Spatial step size
    dt = T / Nt # Time step size
    temps = [] #Result list for temperatures each minute
    minutes = []#Result list for minutes passed
    D0 = elementdict[element]['D0']*1e8 # Diffusion coefficient inital value [um^2/s]
    Ea = elementdict[element]['Ea']# Activation energy for diffusion [eV] 
    rho = elementdict[element]['rho']# density of target [g/cm^3]
    m_a = elementdict[element]['Ma']# atomic mass of target in [g/mole]
    n_atoms = rho*Na/m_a #atomic density of target [atoms/cm^3]
    heat = False
    cool = False
    #-----------------------------------------------------------------
    
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

    # Initialize initial condition
    start = find_start(root)
    if start == None:
        start = 0
    
    # Read data from SRIM
    depth,height = np.loadtxt(root,usecols=(0,1),unpack=True,encoding='cp437') #load height and width of bins
    height = height*fluence #convert into atoms/cm^3
    binwidth = (depth[1]-depth[0])*1e-8 #define binwidth (in cm)
    conc = height/(height + n_atoms) #Calculate concentration from number density of SRIM
    
    
    # Apply initial condition 
    C[0, :] = conc
    C = np.hstack((C, np.zeros((Nt,(L-studyL)*Nx)))) #Add zeros to desired length, SRIM length is only 1 micron usually
    
    
    cool_time = (Temp_fin-1273.15)/10 + time[-1] #Cooling time in minutes
    print(f'Cooling time will be {cool_time/60} hours!!')
    Current_min = 0 #Some awesome parameters
    Fin_min = 0 #same as above
    # Time-stepping loop
    for n in range(0, Nt - 1):
        # Update interior points using forward difference in time and central difference in space
        second = T/Nt*n
        min = int(second/60)
        hour = min/60
        time_left = T/60-min
        temps.append(Temp)
        minutes.append(min)
        if min != Current_min:  # Check if minute has changed
            Current_min = min  # Update current minute
            if heat:
                if Temp <= 1273:
                    Temp = Temp + 10
                    #print('Temp increasing! ')
                    #print(Temp)
                elif Temp < Temp_fin and Temp > 1273:
                    Temp = Temp + 5
                    Fin_min = Current_min
                    #print('Temp increasing! ')
                    #print(Temp)
                else:
                    heat = False     
            if (time_left - cool_time == 0):
                cool =  False
            if cool:
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
    print(f'Diffusion coefficient at maxtemp: {D(D0,Ea,Temp_fin)*1e-8}')
    Concentrations.append(C[-1,:])
    MaxT_Times.append(int(MaxT_time/60))
    
    if Integrate:
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
    Mins.append(minutes)
    Temperatures.append(temps)

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(x,C[0,:], label = 'SRIM simulation')
i = 0
for C_ in Concentrations:        
    plt.plot(x,C_, label = f'Distribution after {Times_in[i]} hours. ')
    i = i + 1

#Measured plot
potku_data = Initialize_Profile(potku_path)
x_pot = potku_data['Samples']['Xe-Imp']['Xe']['x']
c_pot = potku_data['Samples']['Xe-Imp']['Xe']['C']
    
# c_pot2 = potku_data['Samples']['UN-05']['Ru']['C']
# c_pot = [c1 - c2 for c1,c2 in zip(c_pot,c_pot2)]
# c_pot = [0 if c < 0 else c for c in c_pot]
x_pot = [3*1e18*x/(n_atoms) for x in x_pot] #Convert to micrometer
c_pot,x_pot = rebin(c_pot,x_pot)
c_pot,x_pot = rebin(c_pot,x_pot)
c_pot,x_pot = rebin(c_pot,x_pot)
pot_width = (x_pot[1]-x_pot[0])*1e-4 #Convert to cm
pot_Integral = hist_integral(c_pot,pot_width)
print(f'Fluence put in acc. to SRIM:{fluence*0.95} at/cm^2') #0.965 for Zr in UN
print(f'Fluence put in acc. to measurement: {pot_Integral*n_atoms} at/cm^2')

print(f'Ratio: {pot_Integral*n_atoms/(fluence*0.95)}')
plt.plot(x_pot,c_pot, label = 'ToF-ERDA Measurement')
plt.title(f'Comparison between SRIM and ToF-ERDA Measurement ')
plt.xlabel('Position [micrometer]')
plt.ylabel('Concentration [at. fraction]')
plt.grid(True)
plt.legend()

plt.figure()
for i in range(len(Temperatures)):
    plt.title('Temperature change over time')
    plt.xlabel('Minutes')
    plt.ylabel('Temperature [K]')
    plt.plot(Mins[i],Temperatures[i], label = f'Time: {Times_in[i]} h')
plt.legend()
plt.grid(True)
plt.show()
