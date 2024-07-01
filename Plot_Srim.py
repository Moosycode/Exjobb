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


# potku_path = '/Users/niwi9751/potku/requests/20240410-Zr-in-UN.potku'
# data = Initialize_Profile(potku_path)

roots = ['/Users/nilsw/Dropbox/Nils_files/Srim_Results/Xe300keV_in_ZrO2_range.txt',
         '/Users/nilsw/Dropbox/Nils_files/Srim_Results/Kr300keV_in_ZrO2_range.txt',
         '/Users/nilsw/Dropbox/Nils_files/Srim_Results/Fe300keV_in_ZrO2_range.txt']
Names = ['Xe', 'Kr', 'Fe']




fluences = [1e17,1e17,1e17]
i = 0
for root in roots:
    depth, height = np.loadtxt(root,usecols=(0,1),unpack=True,encoding='cp437')
    # rho = 14.05
    # M_a = 252
    rho = 6.025
    M_a = 123.218
    N_a = 6.022e23
    fluence = fluences[i]
    height = height*fluence
    N = N_a*rho/M_a
    conc = height/(height + N)
    plt.plot(depth/10000,conc,label=Names[i])
    i = i+1
plt.grid()
plt.legend()
plt.title('Implanted concentration')
plt.xlabel('Depth [micrometer]')
plt.ylabel('Concentration [at. fraction]')


roots2 = ['/Users/nilsw/Dropbox/Nils_files/Srim_Results/Xe300keV_in_ZrO2_Vacancies.txt',
          '/Users/nilsw/Dropbox/Nils_files/Srim_Results/Kr300keV_in_ZrO2_Vacancies.txt',
          '/Users/nilsw/Dropbox/Nils_files/Srim_Results/Fe300keV_in_ZrO2_Vacancies.txt']

plt.figure()
i = 0
for root in roots2:
    depth, ionvac, recoilvac = np.loadtxt(root,usecols=(0,1,2),unpack=True,encoding='cp437', dtype=str)
    depth = [x.replace(',','.') for x in depth]
    depth = [float(x)/10000 for x in depth]
    ionvac = [x.replace(',','.') for x in ionvac]
    ionvac = [float(x) for x in ionvac]
    recoilvac = [x.replace(',','.') for x in recoilvac]
    recoilvac = [float(x) for x in recoilvac]
    totvac = [x1 + x2 for x1,x2 in zip(ionvac,recoilvac)]
    totvac = [x*1e8 for x in totvac]
    # rho = 14.05
    # M_a = 252
    rho = 6.025
    M_a = 123.218
    N_a = 6.022e23
    fluence = fluences[i]
    N = N_a*rho/M_a
    print(N)
    dpa = [x*fluence/(N) for x in totvac]
    plt.plot(depth,dpa,'-.',label=Names[i])
    i = i+1

plt.grid()
plt.legend()
plt.title('Implanted damage')
plt.xlabel('Depth [micrometer]')
plt.ylabel('dpa')
plt.show()
