# %%
"""
Create dictionary from potku with 

Nils Wikström
2024-02-16

"""

import os
import re
#------------------------------------------------------------------------------
import numpy as np                       ## for handling of data
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
import json
#------------------------------------------------------------------------------
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['font.family'] = 'monospace'
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['k','r','g','b']) 
plt.rcParams['figure.autolayout'] = True
#------------------------------------------------------------------------------
#import ImportData as ID
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# --- element dictionary -----------------------------------
#------------------------------------------------------------------------------
element = {
	'N':{'A':14,'Z':7},
	'Si':{'A':28,'Z':14},
	'Cl':{'A':35,'Z':17},
	'Ti':{'A':48,'Z':22},
	'Cu':{'A':63,'Z':29},
	'Br':{'A':79,'Z':35},
	'Mo':{'A':96,'Z':42},
	'Ag':{'A':107,'Z':47},
	'I':{'A':127,'Z':53},
	'Ta':{'A':181,'Z':73},
	'Au':{'A':197,'Z':79}
	}
#------------------------------------------------------------------------------

# %%
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


# %%

#'C:\Users\Nils\potku\requests\WorkFromAms.potku'
folder_path = '/Users/Nils/potku/requests/WorkFromAms.potku/'
def Initialize_Profile(folder_path):
    dict = {'Beam':{}, 'Samples':{}}
    beamdata = json.load(open(folder_path +'/Default/Default.measurement'))
    beamprofile = json.load(open(folder_path+'/Default/Default.profile'))
    dict['Beam']['Ion'] = re.sub(r'\d+', '', beamdata['beam']['ion'])
    dict['Beam']['Mass'] = re.sub(r'[a-zA-z]',  '' , beamdata['beam']['ion'])
    dict['Beam']['Energy'] = beamdata['beam']['energy']
    
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.info'):
                currentsample = file.removesuffix('.info')
                dict['Samples'][currentsample] = {} 
            if file.startswith('depth.') and 'total' not  in file:
                newroot = root.replace(os.path.sep,  '/') + '/' + file
                depthprofile = file.removeprefix('depth.')
                columns =  read_columns(newroot)
                dict['Samples'][currentsample][depthprofile] = {'Depth': columns[0], 'Concentration': columns[3],'sagas': columns[4],'asg': columns[6]}
                
    return dict

# %%
data = Initialize_Profile(folder_path)

dp_list = [data['Samples']['TiN']['Ti']['Depth'],data['Samples']['TiN']['Ti']['Concentration']]
dp_list2 = [data['Samples']['TiN']['N']['Depth'],data['Samples']['TiN']['N']['Concentration']]
dp_list3 = [data['Samples']['TiN']['Si']['Depth'],data['Samples']['TiN']['Si']['Concentration']]


plt.plot(dp_list[0],dp_list[1])
plt.plot(dp_list2[0],dp_list2[1])
plt.plot(dp_list3[0],dp_list3[1])

plt.xlim([-500,   3000])
plt.ylim([    0.0, 1])
plt.xlabel('depth [10$^{15}$ atoms/cm$^2$]')
plt.ylabel('atomic fraction')
plt.grid(linestyle='--')
plt.legend()



