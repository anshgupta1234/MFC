import numpy as np
import matplotlib.pyplot as plt

# Function to read coordinates from a text file
def read_coordinates(filename):
    return np.loadtxt(filename)

# Read coordinates from gp.txt and ip.txt
gp_coords = read_coordinates('/Users/anshgupta1234/Desktop/Coding/MFC-copy/examples/rectangle/gp.txt')
ip_coords = read_coordinates('/Users/anshgupta1234/Desktop/Coding/MFC-copy/examples/rectangle/ip.txt')

# Create the plot
plt.figure()

# Plot gp points as blue
plt.scatter(gp_coords[:, 0], gp_coords[:, 1], color='blue', label='GP Points')

# Plot ip points as red
plt.scatter(ip_coords[:, 0], ip_coords[:, 1], color='red', label='IP Points')

# Add labels and legend
plt.xlabel('X')
plt.ylabel('Y')
plt.title('GP and IP Points')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Show the plot
plt.show()
