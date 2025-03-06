# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import json
import os
import re
Na = 6.022e23
rho = 6.025
Ma = 123.218
# rho = 14.0
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

def rebinn(N):
    N = [(N[1::2][i] + N[::2][i]) for i in range(len(N[1::2]))]
    return N

def plot_profiles(data):
    i = 0
    for sample in data['Samples']:
        plt.figure(i,figsize=(11,4))
        # plt.title(sample)
        for depth in data['Samples'][sample]:
            if depth in plot_elements:
                x = data['Samples'][sample][depth]['x']
                x = [3*1e21*x/(n_atoms) for x in x]
                C = data['Samples'][sample][depth]['C']
                N = data['Samples'][sample][depth]['N']
                N = rebinn(rebinn(N))
                C,x = rebin(C,x)
                C,x = rebin(C,x)
                C = [c*100 for c in C]
                N = [c/n**(1/2) if n != 0 else 0 for c,n in zip(C,N)]
                plt.yscale('log')
                plt.ylim(0.1,100)
                plt.xlim([-50,300])
                plt.errorbar(x,C, yerr = N, fmt='.k', capsize= 2,capthick=1, ecolor = color_dict[depth])
                plt.plot(x,C, label = depth, color = color_dict[depth])
                plt.xlabel('Depth [nm]',fontsize = 18)
                plt.ylabel(r'Concentration [at.\%]',fontsize = 18)
                plt.grid(linestyle='--')
                plt.legend(loc = 'upper right')
                plt.tight_layout()
        i = i+1

def find_closest_index(arr, target):
    # Use the min function with a custom key to find the closest value in the array
    closest_value = min(arr, key=lambda x: abs(x - target))
    # Find the index of the closest value
    index = arr.index(closest_value)
    return index

def normalize_potku(data,x_start=0,x_stop=250):
    for sample in data['Samples']:
        summ = [0 for i in range(data['Settings']['Num_step']+10)]
        for element in data['Samples'][sample]:
            C = data['Samples'][sample][element]['C']
            summ = [c1 + c2 for c1,c2 in zip(C,summ)]
        for element in data['Samples'][sample]:
            x = data['Samples'][sample][element]['x']
            x = [3*1e21*x/(n_atoms) for x in x]
            x_start_ind = find_closest_index(x,x_start)
            x_stop_ind = find_closest_index(x,x_stop)
            C = data['Samples'][sample][element]['C']
            for i in range(x_start_ind,x_stop_ind):
                C[i] = C[i]/summ[i]
            data['Samples'][sample][element]['C'] = C
    return data

def normalize_potku_2(data):
    for sample in data['Samples']:
        summ = [0 for i in range(data['Settings']['Num_step'])]
        for element in data['Samples'][sample]:
            C = data['Samples'][sample][element]['NoNormC']
            summ = [c1 + c2 for c1,c2 in zip(C,summ)]
        for element in data['Samples'][sample]:
            C = data['Samples'][sample][element]['NoNormC']
            C = [c1/c2 if c2 != 0 else 0 for c1,c2 in zip(C,summ) ]
            data['Samples'][sample][element]['NoNormC'] = C
    return data


color_dict = {'Zr':'r','O':'b', 'Fe':'gray', 'Xe': 'c', 'Kr':'g', 'Hf': 'y', 'Al':'m', 'C':'k', 'Cr': 'silver'}
# color_dict = {'U':'r','N':'deepskyblue', 'O': 'b','Zr':'gray', 'Xe': 'orange', 'Kr':'g', 'Hf': 'y', 'Al':'m', 'C':'k', 'Ru': 'y', 'Ba':'g', 'H':'y'}
plot_elements = ['C', 'Al', 'O', 'Kr', 'Xe', 'Zr', 'Hf', 'Cr']

potku_path = 'Data/ToFERDA/20240304-KrXe-In-ZrO2.potku'

data = Initialize_Profile(potku_path)
data = normalize_potku(data)
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=20)
plot_profiles(data)
plt.show()