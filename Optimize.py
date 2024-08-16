import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from scipy.optimize import minimize
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
root = '/Users/nilsw/Dropbox/Nils_files/Srim_Results/Xe300keV_in_ZrO2_range.txt' 
potku_path = '/Users/nilsw/Dropbox/Nils_Files/Tof_ERDA_Files/requests/20240304-KrXe-In-ZrO2.potku'
# root = '/Users/nilsw/Dropbox/Nils_files/Srim_Results/Xe300keV_in_UN_Range.txt'
# # potku_path = '/Users/niwi9751/potku/requests/20240410-Zr-in-UN.potku'
# potku_path = '/Users/niwi9751/potku/requests/20240506-UNUO2Samples.potku'
# potku_path = '/Users/nilsw/potku/requests/20240304-KrXe-In-ZrO2.potku'
# potku_path = '/Users/niwi9751/potku/requests/20240521-PostAnnealZrO2.potku'
# potku_path = '/Users/nilsw/Potku/requests/20240611UNAuBeam.potku'
#---------------------------------------------------------------------------------

#Global Parameters-----------------------------------------------------------------------
Times_in = 9 #simtime in hrs
Temp = 1473.15 #Target emperature [K] 
fluence = 1e17# Input fluence of implantation [atoms/cm^2]
Concentrations = []#Result list
#---------------------------------------------------------------------------------
T = Times_in*3600 #Convert to seconds
#---------------------------------------------------------------------------------

#Dictionary with needed values of each element
elementdict = {
    'Fe_ZrO2':{'D0':2.26e-6,'Ea':2.3, 'rho':6.025, 'Ma': 123.218}, #Data from Springer
    'Kr_ZrO2':{'D0':8.11e-7,'Ea':2.53, 'rho':6.025, 'Ma': 123.218},
    'Xe_ZrO2':{'D0':1.83e-6,'Ea':2.91, 'rho':6.025, 'Ma': 123.218},
    'Zr_UN':{'D0':6.9e-7,'Ea':2.7, 'rho':14.05, 'Ma': 252.036}, 
    'Kr_UN':{'D0':8.11e-7,'Ea':2.53, 'rho':14.05, 'Ma': 252.036},
    'Xe_UN':{'D0':1.83e-6,'Ea':2.91, 'rho':14.05, 'Ma': 252.036},
    'Zr_UO2':{'D0':6.9e-7,'Ea':2.7, 'rho':10.6, 'Ma': 270.02},
    'Kr_UO2':{'D0':8.11e-7,'Ea':2.53, 'rho':10.6, 'Ma': 270.02}
}
# Parameters------------------------------------------------------
element = 'Xe_ZrO2'
sample = 'Xe'
#-----------------------------------------------------------------
potku_data = Initialize_Profile(potku_path)
x_pot = potku_data['Samples'][f'{sample}-Imp'][sample]['x']
c_pot = potku_data['Samples'][f'{sample}-Imp'][sample]['C']
c_pot,x_pot = rebin(c_pot,x_pot)
c_pot,x_pot = rebin(c_pot,x_pot)
c_pot,x_pot = rebin(c_pot,x_pot)

potku_path2 = '/Users/nilsw/Dropbox/Nils_Files/Tof_ERDA_Files/requests/20240521-PostAnnealZrO2.potku'
data2 = Initialize_Profile(potku_path2)
x2 = data2['Samples'][f'{sample}-Imp'][sample]['x']
c2 = data2['Samples'][f'{sample}-Imp'][sample]['C']
c2,x2 = rebin(c2,x2)
c2,x2 = rebin(c2,x2)
c2,x2 = rebin(c2,x2)

#Constants--------------------------------------------------------
rho = elementdict[element]['rho']# density of target [g/cm^3]
m_a = elementdict[element]['Ma']# atomic mass of target in [g/mole]
n_atoms = rho*Na/m_a #atomic density of target [atoms/cm^3]
Nx = len(x_pot)  # Number of spatial points per micrometer
Nt = int(T/60) # Number of time steps, can be anything really, code finds this for you but do not start lower than this.
x_pot = [3*1e21*x/(n_atoms) for x in x_pot] #Convert to micrometer
Extendby = 1
L = int(x_pot[-1])*Extendby
studyL = int(x_pot[-1])
dx = L / (L*Nx - 1) # Spatial step size
dt = T / Nt # Time step size
x2 = [3*1e21*x/(n_atoms) for x in x2]

def optifunc(vars,plot = False):
    D0, Ea, q = vars
    # Create spatial grid
    x = np.linspace(0, L, Nx*Extendby)
    # Initialize solution matrix
    C = np.zeros((Nt, Nx))
    # Apply initial condition 
    C[0, :] = c_pot
    C = np.hstack((C, np.zeros((Nt,(Extendby-1)*Nx)))) #Add zeros to desired length

    # Time-stepping loop
    for n in range(0, Nt - 1):
        # Update interior points using forward difference in time and central difference in space
        Diff = D(D0, Ea, Temp) #calculate diffusion coefficient
        C[n+1, 0] = C[n, 0] + Diff * dt / dx**2 * (C[n, 1] - 2*C[n, 0])
        for i in range(1, Extendby*Nx - 1): #Update interior points
            C[n+1, i] = C[n, i] + Diff * dt / dx**2 * (C[n, i+1] - 2*C[n, i] + C[n, i-1])
        # Apply Neumann boundary condition to boundaries
        C[n+1, -1] = C[n+1,-2] 

    score = sum([(c-c_)**2 for c,c_ in zip(C[-1,:],c2)])
    
    if plot:
        c_pot2 = [c*100 for c in c_pot]
        Cres = [c*100 for c in C[-1,:]]
        Cmeas = [c*100 for c in c2]
        plt.step(x_pot,c_pot2)
        plt.step(x,Cres)
        plt.step(x,Cmeas)
        plt.xlabel('Depth [nm]', fontsize = 18)
        plt.ylabel('Concentration [at. %]', fontsize = 18)
        plt.grid(True)
        plt.rcParams.update({'font.size':18})
        plt.tight_layout()
        plt.show()
    return score

initial_guess = [elementdict[element]['D0']*1e8, elementdict[element]['Ea'], 0.5]
# bounders = [(None,None),(None,None),(0,None)]

result = minimize(optifunc, initial_guess)

print(result)
optivals = []
for val in result['x']:
    print(val)
    optivals.append(val)
print(optivals)
optifunc(result['x'],True)
