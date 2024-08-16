import numpy as np
import matplotlib.pyplot as plt

#Values for 135-Xe in ZrO2
# Y1 = 4.78
# Y2 = 3.97

#Values for 84-Kr in ZrO2
Y1 = 3.51
Y2 = 2.92



# Y1 = 2.41 #Sputt yield for O [at/ion]
# Y2 = 2.01 #Sputt yield for Zr [at/ion]
Y_avg = (2*Y1 + Y2)/3

fluence = 1e17
rho = 6.05
Na = 6.022e23
Ma = 123.218
n_atoms = rho*Na/Ma

dist_TFU = (2*Y1 + Y2)/3*fluence

dist = dist_TFU*Ma/(3*rho*Na)

print(dist*1e-2/(1e-9))

root = '/Users/nilsw/Dropbox/Nils_files/Srim_Results/Kr300keV_in_ZrO2_Range.txt'
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

depth,height = np.loadtxt(root,usecols=(0,1),unpack=True,encoding='cp437') #load height and width of bins
height = height*fluence #convert into atoms/cm^3
binwidth = (depth[1]-depth[0])*1e-8 #define binwidth (in cm)
conc = height/(height + n_atoms) #Calculate concentration from number density of SRIM

conc_shift = shift_and_average(conc)
# conc_shift = [c*100/sum(conc_shift) for c in conc_shift]
conc_shift = [c*100 for c in conc_shift]
x_shift = np.linspace(0,1.03,103)
x_shift = [(x-0.037)*1000 for x in x_shift]
x_shift = np.array(x_shift)
# conc = [c*100/sum(conc) for c in conc]
conc = [c*100 for c in conc]
x = np.linspace(0,1000,100)
plt.rcParams.update({'font.size':18})
plt.plot(x_shift,conc_shift,label='SRIM - corrected')
plt.plot(x,conc, label = 'SRIM')
plt.axvline(x = 0, label = 'Surface', color = 'g',linestyle = '--')
plt.fill_between(x_shift, conc_shift, where=(x_shift >= -40) & (x_shift <= 5), color='skyblue', label = 'Sputtered ions')
plt.xlim([-50,450])
Sputt = sum(conc_shift[0:4])
print(Sputt)
plt.xlabel('Depth [nm]')
plt.ylabel('Concentration [at. %]')
plt.tight_layout()
plt.legend()
plt.grid(True)
plt.show()