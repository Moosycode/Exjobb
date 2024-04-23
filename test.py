import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import datetime
import json
import os
import re

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

potku_path = '/Users/niwi9751/potku/requests/20240304-KrXe-In-ZrO2.potku'


potku_data = Initialize_Profile(potku_path)
x_pot = potku_data['Samples']['Fe-Imp']['Fe']['x']
x_pot = [3*1e18*x/(n_atoms) for x in x_pot] #Convert to micrometer
c_pot = potku_data['Samples']['Fe-Imp']['Fe']['C']
plt.plot(x_pot,c_pot, label = 'ToF-ERDA Measurement')
plt.title(f'Diffusion of Concentration ')
plt.xlabel('Position [micrometer]')
plt.ylabel('Concentration [at. fraction]')
plt.grid(True)
plt.legend()
