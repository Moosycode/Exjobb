import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import datetime
import json
import os
import re

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
    
root = '/Users/niwi9751/CONTES/Output_Data/I127_44MeV_pos7_Zr_imp_UN_AIMED_HIGH_profiles.txt'
    
data = np.loadtxt(root,skiprows=3, unpack=True,encoding='cp437')
depth = data[0]
print(data)

# elements = ['U', 'Al', 'O', 'N',]
# for i in range(len(data)-1):
#     if i != 0:
#         plt.step(depth,data[i],label = elements[i-1])
# plt.legend()
# plt.show()
