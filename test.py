import numpy as np

# Path to the file
file_path = '/Users/Nils/Testing.txt'

# Use np.genfromtxt to read the headers separately
headers = np.genfromtxt(file_path, max_rows=1, dtype=str, skip_header=2)

# Now use np.loadtxt to read the data, skipping the first four lines including the header line
data = np.loadtxt(file_path, skiprows=3)

# Creating a dictionary to map headers to data columns
structured_data_np = {header: data[:, i] for i, header in enumerate(headers)}

# You can print structured_data_np to see the output:
print(structured_data_np['Zr'])
