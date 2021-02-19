# A) Using the planetary data in the accompanying document (excel sheet titled “MR_data”,
# make a plot mass-radius plot (similar to Fig 11.8 in the book) with different lines indicating
# planets with a pure Fe, pure rock, and water-ice mixture. Please make sure to label your axes,
# including units, as well as the 3 lines. Also, label the region above the water-ice line as
# “atmosphere,” between rock and water-ice as “likely surface volatiles,” and between rock
# and pure Fe as “Fe-rich.”

# B) Filter the data to only include rocky exoplanets, with a radius < 1.6 R_Earth. Then choose 3
# rocky planets (Graduates: choose 5). Calculate the density of these planets in g/cm3. Add these
# rocky planets to your plot and label with them with their name and densities.
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['figure.figsize'] = (10, 8)

# Water/Ice data read
waveO1, flux1 = np.loadtxt('DAsmus_Sy1.txt', unpack = True)
# Pure rock data read:
waveO2, flux2 = np.loadtxt('DAsmus_Sy2.txt', unpack = True)
# Pure Fe (Liquid) data read:
waveO3, flux3 = np.loadtxt('DAsmus_Syint.txt', unpack = True)
# Atmosphere boundary data read:
#waveO4, flux4 = np.loadtxt('HW3atmo.txt', unpack = True)
# Rock boundary data read:
#waveO5, flux5 = np.loadtxt('HW3rocklayer.txt', unpack = True)
# Fe-Rich boundary data read:
#waveO6, flux6 = np.loadtxt('HW3Fe_rich.txt', unpack = True)

# plot main sequence lines
plt.scatter(waveO1, flux1, color = 'y', linewidth = 2, label = 'Sy 1')
plt.scatter(waveO2, flux2, color = 'b', linewidth = 2, label = 'Sy 2', marker = 'd')
plt.scatter(waveO3, flux3, color = 'r', linewidth = 2, Label = 'Sy-Intermediate', marker = 's')


plt.xlabel('Log F [OIV] (ergs s^-1 cm^-2)')
plt.ylabel('Distance (Mpc)')
plt.title('Log F [OIV] vs Distance')
plt.legend()
#plt.savefig('dasmus_2016.png')
plt.show()