import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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

def hist_integral(n, width):
    n = [item*width for item in n]
    return sum(n)#Definition of integrals :))

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

#Parameters 
potku_path = '/Users/niwi9751/potku/requests/20240410-Zr-in-UN.potku' #Folder path to potku file (ends with .potku)
x_start = 0 #Start value of integral (in 10^15 at/cm^2)
x_end = 1000 # stop value of integral (in 10^15 at/cm^2)

#-------------------Example usage-------------------------
data = Initialize_Profile(potku_path)

x = data['Samples']['UN-AimedLow']['U']['x']
C = data['Samples']['UN-AimedLow']['U']['C']
binwidth = x[1]-x[0]

plt.step(x,C)
plt.xlabel('depth [10$^{15}$ atoms/cm$^2$]')
plt.ylabel('atomic fraction')
plt.grid(linestyle='--')
plt.show()

diffarray1 = [abs(x - x_start) for x in x]
diffarray1 = np.array(diffarray1)
diffarray2 = [abs(x - x_end) for x in x]
diffarray2 = np.array(diffarray2)
start_index = diffarray1.argmin()
stop_index = diffarray2.argmin()

integral = hist_integral(C[start_index:stop_index],binwidth)

print(f'Integral between {x[start_index]} and {x[stop_index]} is: {integral}')