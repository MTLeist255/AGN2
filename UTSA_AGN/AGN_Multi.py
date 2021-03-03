# Create multi plots of the final redshift correct wavelength and 9/18um silicate strength fit

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# load redshift correct data
wave1, flux1, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_3c273.txt', unpack = True)
wave2, flux2, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_3c321.txt', unpack = True)
wave3, flux3, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_cenA.txt', unpack = True)
wave4, flux4, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_circinus.txt', unpack = True)
wave5, flux5, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_cygnusA.txt', unpack = True)
wave6, flux6, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_eso323g77.txt', unpack = True)

# plot all data on single figure
plt.figure(figsize=(16,10))
plt.subplot(2,3,1)
plt.plot(wave1, flux1, color = 'k')
plt.xlabel('Wavelength[um]')
plt.ylabel('Flux[Jy]')
plt.title('3c273')

plt.subplot(2,3,2)
plt.plot(wave2, flux2, color = 'k')
plt.xlabel('Wavelength[um]')
plt.ylabel('Flux[Jy]')
plt.title('3c321')

plt.subplot(2,3,3)
plt.plot(wave3, flux3, color = 'k')
plt.xlabel('Wavelength[um]')
plt.ylabel('Flux[Jy]')
plt.title('Centaurus A')

plt.subplot(2,3,4)
plt.plot(wave4, flux4, color = 'k')
plt.xlabel('Wavelength[um]')
plt.ylabel('Flux[Jy]')
plt.title('Circinus')

plt.subplot(2,3,5)
plt.plot(wave5, flux5, color = 'k')
plt.xlabel('Wavelength[um]')
plt.ylabel('Flux[Jy]')
plt.title('Cygnus A')

plt.subplot(2,3,6)
plt.plot(wave6, flux6, color = 'k')
plt.xlabel('Wavelength[um]')
plt.ylabel('Flux[Jy]')
plt.title('ESO 323-G77')
#plt.savefig('Current_AGN_fits/Wave_Corrected/test1.png')
plt.show()

# load 9/18um silicate strength data

fig = plt.figure(figsize=(12,18))
imT = Image.open('Current_AGN_fits/silicate/9_18_continuum/test.png') # bottom image
imM = Image.open('Current_AGN_fits/silicate/9_18_continuum/test3.png') # middle image
imB = Image.open('Current_AGN_fits/silicate/9_18_continuum/test4.png') # top image

# plot all 3 images on a single figure
# rows read from bottom up: 1 -> 3
# Display row 1 -> BOTTOM
im1 = plt.figimage(imB, xo=0, yo=0)
# Display row 2 -> MIDDLE
im2 = plt.figimage(imM, xo=0, yo=600)
# Display row 3 -> TOP
im3 = plt.figimage(imT, xo=0, yo=1200)
#plt.savefig('Current_AGN_fits/silicate/continuum_stacked/test.png')
plt.show()