import numpy as np
import matplotlib.pyplot as plt

# Define the x values and the function (let's use a sine wave as an example)
x = np.linspace(-100, 100, 1000)
y = np.sin(x)

# Create the plot
plt.plot(x, y, label='sin(x)')

# Fill the area between x = -50 and x = 0
plt.fill_between(x, y, where=(x >= -50) & (x <= 0), color='skyblue', alpha=0.5)

# Add labels and a legend
plt.xlabel('x')
plt.ylabel('y')
plt.title('Filling area between x = -50 and x = 0')
plt.legend()

# Show the plot
plt.show()
