import numpy as np

def find_row_with_hyphen(filename):
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            for i, line in enumerate(file, 1):
                if '--' in line:
                    # Check if the line contains a hyphen
                    return i  # Return the index (line number) if found
    except Exception as e:
        print("An error occurred while reading the file:", e)
        return None

# Example usage
root = '/Users/niwi9751/Srim_Results/Fe_inZrO2_300keV.txt'
matching_row_index = find_row_with_hyphen(root)
if matching_row_index:
    print("Found row containing a hyphen at index:", matching_row_index)
else:
    print("No row containing a hyphen found.")
    
# num, depth = np.loadtxt(root, usecols=(0,1),unpack=True,skiprows=18)
# print(num)

