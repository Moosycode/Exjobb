# Importing necessary library
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import json
import os
import re


# Function to clean headers by removing parentheses and contents within
def clean_header(header):
    return re.sub(r"\(.*?\)", "", header).strip()

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

# Load the data
file_path = '/Users/niwi9751/CONTES/Output_Data/I127_44MeV_pos7_Zr_imp_UN_AIMED_LOW_profiles.txt'
potku_path = '/Users/niwi9751/potku/requests/20240410-Zr-in-UN.potku'
data = pd.read_csv(file_path, delimiter='\t', skiprows=3, header=None)

# Extract headers and data
headers = data.iloc[0].values
values = data.iloc[1:]



# Create the dictionary with cleaned headers as keys and convert each column's values to appropriate numeric types
data_dict = {}
for index, header in enumerate(headers):
    clean_header_name = clean_header(header)
    headers[index] = clean_header_name
    # Attempt to convert column to float, if fails keep as string (this should not happen if all are numbers)
    try:
        data_dict[clean_header_name] = values[index].astype(float).tolist()
    except ValueError:
        data_dict[clean_header_name] = values[index].tolist()

data_dict.popitem()
headers = np.delete(headers, len(headers)-1)
headers = np.delete(headers, 0)
x_cont = data_dict['x']
c_cont = data_dict['Zr']


potku_data = Initialize_Profile(potku_path)
x_pot = potku_data['Samples']['UN-AimedLow']['Zr']['x']
c_pot = potku_data['Samples']['UN-AimedLow']['Zr']['C']
plt.step(x_pot,c_pot,label = 'Potku')
plt.step(x_cont,c_cont,label = 'Contes')
plt.legend()    
plt.show()
