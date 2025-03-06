import numpy as np
import matplotlib.pyplot as plt
import os

#Values for 135-Xe in ZrO2
Y1 = 4.78
Y2 = 3.97

#Values for 84-Kr in ZrO2
# Y1 = 3.51
# Y2 = 2.92



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
    return sum(n[0:1])#Definition of integrals :))

def shift_and_average(data, num_copies=5):
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
x_shift = np.linspace(0,1.04,104)
x_shift = [(x-0.04)*1000 for x in x_shift]
x_shift = np.array(x_shift)
# conc = [c*100/sum(conc) for c in conc]
conc = [c*100 for c in conc]
x = np.linspace(0,1000,100)
plt.figure('Xe',figsize=(7,5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=20)
plt.title('Kr')
plt.plot(x_shift,conc_shift,label='Corrected')
plt.plot(x,conc, label = 'SRIM')
plt.axvline(x = 0, label = 'Surface', color = 'g',linestyle = '--')
plt.fill_between(x_shift, conc_shift, where=(x_shift >= -60) & (x_shift <= 4), color='skyblue', label = 'Sputtered')
plt.xlim([-50,300])
Sputt = sum(conc_shift[0:4])
print(Sputt)

sample ='Kr'
potku_path = 'Data/ToFERDA/20240304-KrXe-In-ZrO2.potku'
potku_data = Initialize_Profile(potku_path)
potku_data = normalize_potku(potku_data)
x_pot = potku_data['Samples'][f'{sample}-Imp'][sample]['x']
c_pot = potku_data['Samples'][f'{sample}-Imp'][sample]['C']
N = potku_data['Samples'][f'{sample}-Imp'][sample]['N']
N = rebinn(N)
c_pot = [c*100 for c in c_pot]
x_pot = [3*1e21*a/(n_atoms) for a in x_pot]
# N = [c/n**(1/2) if n != 0 else 0 for c,n in zip(c_pot,N)]
# print(N)
c_pot,x_pot = rebin(c_pot,x_pot)
# c_pot,x_pot = rebin(c_pot,x_pot)
# plt.errorbar(x_pot,c_pot, yerr = N, fmt='.k', capsize= 2,capthick=1, ecolor = 'g')
plt.plot(x_pot,c_pot,label = 'ToF-ERDA')
plt.xlabel('Depth [nm]')
plt.ylabel(r'Concentration [at.\%]')
plt.tight_layout()
plt.legend()
plt.grid(True)
plt.show()