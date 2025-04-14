import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
# #Values for 135-Xe in ZrO2
Y1 = 4.78
Y2 = 3.97

# Values for 84-Kr in ZrO2
Y1 = 3.51
Y2 = 2.92



# Y1 = 2.41 #Sputt yield for O [at/ion]
# Y2 = 2.01 #Sputt yield for Zr [at/ion]
Y_avg = (2*Y1 + Y2)/3
print(Y_avg)

fluence = 1e17
rho = 6.05
Na = 6.022e23
Ma = 123.218
n_atoms = rho*Na/Ma

dist_TFU = (2*Y1 + Y2)/3*fluence

dist = dist_TFU*Ma/(3*rho*Na)

print(dist*1e-2/(1e-9))

root = 'Data/SRIM/Kr300keV_in_ZrO2_Range.txt'

# pre_data_path = 'Data/ToFERDA/Data/pre-anneal_Xe-imp.csv'
# post_data_path = 'Data/ToFERDA/Data/post-anneal_Xe-imp.csv'
pre_data_path = 'Data/ToFERDA/Data/pre-anneal_Kr-imp_corrected.csv'
post_data_path = 'Data/ToFERDA/Data/post-anneal_Kr-imp_corrected.csv'

def hist_integral(n, width):
    n = [item*width for item in n]
    return sum(n[0:1])#Definition of integrals :))

def shift_and_average(data, num_copies=4):
    """
    Copies the dataset `num_copies` times, shifts each copy incrementally to the right,
    sums them, and then averages the values.
    
    Parameters:
    data (np.ndarray): The original dataset.
    num_copies (int): The number of copies to create and shift.
    
    Returns:
    np.ndarray: The averaged dataset.
    """
    total_length = len(data)
    
    # Initialize the sum array with zeros
    summed_data = np.zeros(total_length + num_copies - 1)
    
    for i in range(num_copies):
        # Create a new array for shifted data
        shifted_data = np.zeros(total_length + num_copies - 1)
        # Insert the data shifted by 'i' positions to the right
        shifted_data[i:i + total_length] = data
        
        # Add the shifted dataset to the sum array
        summed_data += shifted_data
    # Calculate the average

    averaged_data = summed_data / num_copies
    
    return averaged_data

def find_closest_index(arr, target):
    # Use the min function with a custom key to find the closest value in the array
    closest_value = min(arr, key=lambda x: abs(x - target))
    # Find the index of the closest value
    index = arr.index(closest_value)
    return index

def normalize_potku(data,x_start=0,x_stop=350):
    for sample in data['Samples']:
        summ = [0 for i in range(85)]
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

def rebinn(N):
    N = [(N[1::2][i] + N[::2][i]) for i in range(len(N[1::2]))]
    return N


def rebin(data,x_res):
    data = [(data[1::2][i] + data[::2][i])/2 for i in range(len(data[1::2]))]
    x_res = x_res[1::2]
    return data, x_res


depth,height = np.loadtxt(root,usecols=(0,1),unpack=True,encoding='cp437') #load height and width of bins
height = height*fluence #convert into atoms/cm^3
binwidth = (depth[1]-depth[0])*1e-8 #define binwidth (in cm)
conc = height/(height + n_atoms) #Calculate concentration from number density of SRIM

conc_shift = shift_and_average(conc)
# conc_shift = [c*100/sum(conc_shift) for c in conc_shift]
conc_shift = [c*100 for c in conc_shift]
x_shift = np.linspace(0,1.03,103)
x_shift = [(x-0.04)*1000 for x in x_shift]
x_shift = np.array(x_shift)
# conc = [c*100/sum(conc) for c in conc]
conc = [c*100 for c in conc]
x = np.linspace(0,1000,100)
plt.figure('Xe',figsize=(7,5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=16)
plt.title(r'\textbf{Kr}')
plt.plot(x_shift,conc_shift,label='Corrected')
plt.plot(x,conc, label = 'SRIM')
plt.axvline(x = 0, label = 'Surface', color = 'g',linestyle = '--')
plt.fill_between(x_shift, conc_shift, where=(x_shift >= -60) & (x_shift <= 4), color='skyblue', label = 'Sputtered')
plt.xlim([-50,300])
Sputt = sum(conc_shift[0:4])
print(Sputt)

sample ='Kr'
df_pre = pd.read_csv(pre_data_path)
df_pre = df_pre[df_pre['depth'] <= 3000]
c_pre = df_pre[sample]
err_pre = df_pre[f'{sample}_sd']
tot_err_pre = (sum([e**2 for e in err_pre]))**0.5
x_pre = df_pre['depth']
df_post = pd.read_csv(post_data_path)
df_post = df_post[df_post['depth'] <= 3000]
c_post = df_post[sample]
err_post = df_post[f'{sample}']
c_pre = [c*100 for c in c_pre]
x_pre = [3*1e21*a/(n_atoms) for a in x_pre]
err_pre = [e*100 for e in err_pre]

plt.errorbar(x_pre,c_pre,yerr= err_pre, ecolor='green', label = 'ToF-ERDA', color ='green')
plt.xlabel('Depth [nm]')
plt.ylabel(r'Concentration [at.\%]')
plt.tight_layout()
plt.legend()
plt.grid(True)
plt.show()