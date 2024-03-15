import numpy as np
import datetime
import matplotlib.pyplot as plt

root = '/Users/nilsw/Python/Exjobb/verticalfurnaceupto1600.txt'
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
    return hours + minutes/60 + seconds/3600

time = [t.replace(',','.') for t in time]
time = np.array(time, dtype=float)
time = [excel_time_to_hours(t) for t in time] 
time = [t - time[0] for t in time]
data = [d.replace(',','.') for d in data]
data = np.array(data, dtype=float)


plt.plot(time,data)
plt.show()
