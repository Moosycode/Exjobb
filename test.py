import numpy as np
import datetime
import matplotlib.pyplot as plt

root = '/Users/niwi9751/furnanceData.txt'
time, data = np.loadtxt(root,usecols=(0,1),unpack=True, dtype= 'U25')

def excel_time_to_hours(excel_time):
    # Extracting the integer part as days
    days = int(excel_time)
    # Extracting the fractional part as seconds
    seconds = (excel_time - days) * 24 * 3600
    # Converting to timedelta
    time_delta = datetime.timedelta(seconds=seconds)
    # Extracting hours, minutes, and seconds
    hours = int(time_delta.total_seconds() // 3600)
    minutes = int((time_delta.total_seconds() % 3600) // 60)
    seconds = int(time_delta.total_seconds() % 60)
    return hours*60 + minutes + seconds/60

time = [t.replace(',','.') for t in time]
time = np.array(time, dtype=float)
time = [excel_time_to_hours(t) for t in time] 
time = [int(t - time[0]) for t in time]
data = [d.replace(',','.') for d in data]
data = np.array(data, dtype=float)
print(time)
plt.plot(time,data)
plt.show()