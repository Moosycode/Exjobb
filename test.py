import numpy as np

def structure_data_with_numpy(filepath):
    # First, open the file to extract headers properly
    with open(filepath, 'r') as file:
        lines = file.readlines()
        headers = None
        skip_rows = 0
        for line in lines:
            if line.startswith('%'):
                skip_rows += 1
                continue
            if headers is None:
                headers = line.strip().split()
                break
            
    print(skip_rows)
    print(headers)
    # Load data, skipping the necessary rows and handling spaces as delimiters
    data = np.loadtxt(filepath, delimiter=None, comments='%', skiprows=skip_rows)
    
    # Create a dictionary where keys are headers and values are columns of data
    data_dict = {header: data[:, i].tolist() for i, header in enumerate(headers)}
    
    return data_dict

# Specify the path to the data file


file_path = '/Users/niwi9751/CONTES/Output_Data/I127_44MeV_pos7_Zr_imp_UN_AIMED_HIGH_profiles.txt'

# Use the function to structure the data
structured_data = structure_data_with_numpy(file_path)
print(structured_data)