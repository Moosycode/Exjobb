import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import json
import os
import re
Na = 6.022e23
rho = 6.025
Ma = 123.218
# rho = 14.05
# Ma = 252.036
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
            C,x = rebin(C,x)
            # C,x = rebin(C,x)
            C = [c*100 for c in C]
            plt.yscale('log')
            plt.ylim(0.1,100)
            plt.xlim([0,400])
            plt.plot(x,C, label = depth, color = color_dict[depth])
            plt.xlabel('depth [nanometer]')
            plt.ylabel('atomic %')
            plt.grid(linestyle='--')
            plt.legend(loc = 'upper right')
        i = i+1

def hist_integral(n, width):
    n = [item*width for item in n]
    return sum(n)#Definition of integrals :))


Na = 6.022e23
rho = 6.025
Ma = 123.218
# rho = 14.05
# Ma = 252.036
n_atoms = rho*Na/Ma

# potku_path = '/Users/niwi9751/potku/requests/20240410-Zr-in-UN.potku'
# potku_path = '/Users/niwi9751/potku/requests/20240506-UNUO2Samples.potku'
potku_path1 = '/Users/niwi9751/potku/requests/20240304-KrXe-In-ZrO2.potku'
# potku_path = '/Users/niwi9751/potku/requests/20230205_KrXe_in_ZrO2.potku'
potku_path2 = '/Users/niwi9751/potku/requests/20240521-PostAnnealZrO2.potku'
# potku_path2 = '/Users/niwi9751/potku/requests/20240319-Fe-In-ZrO2.potku'

data = Initialize_Profile(potku_path1)
data2 = Initialize_Profile(potku_path2)
elements = ['Xe', 'Kr', 'Fe']
element = elements[1]

x1 = data['Samples'][f'Kr-Imp'][f'{element}']['x']
x1 = [3*1e21*x/(n_atoms) for x in x1]
c1 = data['Samples'][f'Kr-Imp'][f'{element}']['C']
c1,x1 = rebin(c1,x1)
c1,x1 = rebin(c1,x1)
c1 = [c for c in c1]

x2 = data2['Samples'][f'Kr-Imp'][f'{element}']['x']
x2 = [3*1e21*x/(n_atoms) for x in x2]
c2 = data2['Samples'][f'Kr-Imp'][f'{element}']['C']
c2,x2 = rebin(c2,x2)
c2,x2 = rebin(c2,x2)
c2 = [c for c in c2]

binwidth = (x1[1]-x1[0])*1e-7

preInt = hist_integral(c1,binwidth)
postInt = hist_integral(c2,binwidth)

print(preInt*n_atoms)
print(postInt*n_atoms)


plt.xlabel('depth [nanometer]')
plt.ylabel('atomic %')
plt.grid(linestyle='--')
plt.plot(x1,c1, label = f'Implanted {element}, pre annealing')
plt.plot(x2,c2, label = f'Implanted {element}, post annealing')
plt.title('Implanted Xenon, before and after annealing, normalised')
plt.legend()
plt.show()
