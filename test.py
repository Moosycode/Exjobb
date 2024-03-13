import numpy as np
import matplotlib.pyplot as plt

root = '/Users/niwi9751/verticalfurnaceupto1600.txt'
data = np.loadtxt(root,usecols=1,unpack=True, dtype='U25', )

def remove_before_semi(text_list):
    result = []
    for item in text_list:
        semicolon_index = item.find(';')
        if semicolon_index != -1:  # If a semicolon is found
            number_part = item[semicolon_index + 1:].strip()  # Extract the part after semicolon and remove leading/trailing whitespace
            if number_part.isdigit():  # Check if the remaining part is a valid integer
                result.append(int(number_part))
    return result

def remove_text_until_zero_and_convert_to_int(text_list):
    result = []
    for item in text_list:
        semicolon_index = item.find(';')
        if semicolon_index != -1:  # If a semicolon is found
            text_before_semicolon = item[:semicolon_index]  # Extract text before semicolon
            zero_index = text_before_semicolon.find('0')  # Find index of first '0'
            if zero_index != -1:  # If '0' is found before semicolon
                number_part = text_before_semicolon[:zero_index].strip()  # Extract the part before '0' and remove leading/trailing whitespace
                if number_part.isdigit():  # Check if the remaining part is a valid integer
                    result.append(int(number_part))
    return result

print(remove_before_semi(data))
x_vals = np.linspace(0,len(data),len(data))

print(x_vals)
