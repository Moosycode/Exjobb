import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import json
import os
import re
Na = 6.022e23
# rho = 6.025
# Ma = 123.218
rho = 14.05
Ma = 252.036
n_atoms = rho*Na/Ma

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

def plot_profiles(data):
    i = 0
    for sample in data['Samples']:
        plt.figure(i)
        plt.title(sample)
        for depth in data['Samples'][sample]:
            x = data['Samples'][sample][depth]['x']
            x = [3*1e21*x/(n_atoms) for x in x]
            C = data['Samples'][sample][depth]['C']
            C,x = rebin(C,x)
            # C,x = rebin(C,x)
            C,x = rebin(C,x)
            C = [c*100 for c in C]
            # plt.yscale('log')
            plt.ylim(0,1)
            plt.xlim([0,400])
            plt.plot(x,C, label = depth, color = color_dict[depth])
            plt.xlabel('depth [nanometer]')
            plt.ylabel('atomic %')
            plt.grid(linestyle='--')
            plt.legend(loc = 'upper right')
        i = i+1

def plot_profiles_noNorm(data):
    i = 0
    for sample in data['Samples']:
        plt.figure(i)
        plt.title(sample)
        for depth in data['Samples'][sample]:
            x = data['Samples'][sample][depth]['x']
            x = [3*1e21*x/(n_atoms) for x in x]
            C = data['Samples'][sample][depth]['NoNormC']
            C,x = rebin(C,x)
            # C,x = rebin(C,x)
            C,x = rebin(C,x)
            C = [c*100 for c in C]
            # plt.yscale('log')
            # plt.ylim(0,1)
            plt.xlim([0,400])
            plt.plot(x,C, label = depth, color = color_dict[depth])
            plt.xlabel('depth [nanometer]')
            plt.ylabel('atomic %')
            plt.grid(linestyle='--')
            plt.legend(loc = 'upper right')
        i = i+1


def normalize_potku(data):
    for sample in data['Samples']:
        summ = [0 for i in range(data['Settings']['Num_step']+10)]
        for element in data['Samples'][sample]:
            C = data['Samples'][sample][element]['C']
            summ = [c1 + c2 for c1,c2 in zip(C,summ)]
        for element in data['Samples'][sample]:
            C = data['Samples'][sample][element]['C']
            C = [c1/c2 if c2 != 0 else 0 for c1,c2 in zip(C,summ) ]
            data['Samples'][sample][element]['C'] = C
    return data

def normalize_potku_2(data):
    for sample in data['Samples']:
        summ = [0 for i in range(data['Settings']['Num_step']+10)]
        for element in data['Samples'][sample]:
            C = data['Samples'][sample][element]['NoNormC']
            summ = [c1 + c2 for c1,c2 in zip(C,summ)]
        for element in data['Samples'][sample]:
            C = data['Samples'][sample][element]['NoNormC']
            C = [c1/c2 if c2 != 0 else 0 for c1,c2 in zip(C,summ) ]
            data['Samples'][sample][element]['NoNormC'] = C
    return data

# color_dict = {'Zr':'r','O':'b', 'Fe':'gray', 'Xe': 'c', 'Kr':'g', 'Hf': 'y', 'Al':'m', 'C':'k', 'Cr': 'silver'}
color_dict = {'U':'r','N':'deepskyblue', 'O': 'b','Zr':'gray', 'Xe': 'c', 'Kr':'g', 'Hf': 'y', 'Al':'m', 'C':'k', 'Br': 'y','Ru':'y'}

# potku_path = '/Users/niwi9751/potku/requests/20240410-Zr-in-UN.potku'
potku_path = '/Users/niwi9751/potku/requests/20240506-UNUO2Samples.potku'
# potku_path = '/Users/niwi9751/potku/requests/20240304-KrXe-In-ZrO2.potku'
# potku_path = '/Users/niwi9751/potku/requests/20230205_KrXe_in_ZrO2.potku'
# potku_path = '/Users/niwi9751/potku/requests/20240521-PostAnnealZrO2.potku'
# potku_path = '/Users/niwi9751/potku/requests/20240319-Fe-In-ZrO2.potku'


root = '/Users/niwi9751/Dropbox/Nils_files/Srim_Results/Zr330keV_in_UN_Range.txt'
fluence = 9.7e15
x_srim = np.linspace(0,1,100)
# Read data from SRIM
depth,height = np.loadtxt(root,usecols=(0,1),unpack=True,encoding='cp437') #load height and width of bins
height = height*fluence #convert into atoms/cm^3
binwidth = (depth[1]-depth[0])*1e-8 #define binwidth (in cm)
conc = height/(height + n_atoms) #Calculate concentration from number density of SRIM
plt.plot(x_srim,conc)
    

data = Initialize_Profile(potku_path)
samples = ['UN-05', 'UN-1-MIT']
elements = ['Kr', 'Zr']
sample = samples[1]
element = elements[1]
x = data['Samples'][f'{sample}'][f'{element}']['x']
x = [3*1e18*x/(n_atoms) for x in x] #Convert to micrometer
C1 = data['Samples'][f'{sample}'][f'{element}']['NoNormC']
C2 = data['Samples'][f'{sample}']['Kr']['NoNormC']
C = [c1 - c2 for c1,c2 in zip(C1,C2)]
data['Samples'][f'{sample}'][f'{element}']['NoNormC'] = C
data = normalize_potku_2(data)
C = data['Samples'][f'{sample}'][f'{element}']['NoNormC']
C, x = rebin(C,x)
C, x = rebin(C,x)
C, x = rebin(C,x)

plt.plot(x,C, label = 'asfasdf')
plt.xlabel('Position [micrometer]')
plt.ylabel('Concentration [at. fraction]')
plt.grid(True)
plt.legend()
plt.show()
