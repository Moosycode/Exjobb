import numpy as np
root = '/Users/niwi9751/Srim_Results/wubelibubeli.txt'
def find_start(filename):
    try:
        with open(filename, 'r', encoding='cp437') as file:
            for i, line in enumerate(file, 1):
                if '--' in line:
                    # Check if the line contains a hyphen
                    return i  # Return the index (line number) if found
    except Exception as e:
        print("An error occurred while reading the file:", e)
        return None



start = find_start(root)
data = np.loadtxt(root,skiprows=start,usecols=1,unpack=True,encoding='cp437')
print(data)


