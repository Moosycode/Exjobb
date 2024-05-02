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
    
# Step 1: Read the file and identify headers and data
file_path = '/Users/Nils/Testing.txt'

# Initialize a dictionary to hold the structured data
structured_data = {}

# Read the file and process it
with open(file_path, 'r') as file:
    lines = file.readlines()
    
    # The 4th line contains the headers, which are 1-indexed in a normal description, 0-indexed in Python
    headers = lines[2].split()
    print(headers)
    
    # Initialize the dictionary with headers as keys and empty lists as values
    for header in headers:
        structured_data[header] = []
    
    # Iterate over the remaining lines to populate the dictionary
    for line in lines[3:]:
        if line.strip():  # This checks if the line is not just whitespace
            values = line.split()
            for header, value in zip(headers, values):
                structured_data[header].append(float(value))  # Convert string to float and append to the correct list

# The structured_data dictionary now contains the data categorized under each header.
# For this code to run here, you'd need to uncomment the next line:
print(structured_data['Fe'])

# Let's comment out the print for now as per instruction to optimize the code further if necessary.
# print(structured_data)
