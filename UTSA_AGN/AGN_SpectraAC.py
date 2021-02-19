# Plotting the CORRECTED wavelength and flux for spectral cassis,
# herschel and hawc data for several various AGN in our sample set.
# To run this program with available data, just un-comment which
# section you wish to use. Finally, once you are comfortable with the
# plot, uncomment the 'save.fig' function at the end of the code. Once
# the fig is saved, moved the completed plot to the correct folder for
# filing and reviewing.

import matplotlib.pyplot as plt
import numpy as np
import numpy as np


plt.rcParams['figure.figsize'] = (10, 8)

# Cassis data: awaveO = wavelength observed, flux = flux, err = error
# BREAK 1:
waveO1, flux1, err1 = np.loadtxt('CASSIS\CASSIS_BC\Break1\cassis_3c273_break1.txt', unpack = True)
# BREAK 2:
waveO2, flux2, err2 = np.loadtxt('CASSIS\CASSIS_BC\Break2\cassis_3c273_break2.txt', unpack = True)
# BREAK 3: USE AS NEEDED -> COMMENT OUT WHEN UN-NEEDED
#waveO3, flux3, err3 = np.loadtxt('CASSIS\CASSIS_BC\Break3\cassis_ngc3783_break3.txt', unpack = True)

# input hawc data
#d, e, f, g = np.loadtxt('HAWC\hawcngc4388.txt', unpack = True)

# input herschel data
#h, i = np.loadtxt('HERSHEL\herschelngc4388.txt', unpack = True)

# WAVELENGTH CORRECTION: correct 5-15um and 15-38um
# wavelength(emitted)=wavelength(observed)/1+z
# z = redshift (value found in spreadsheet)
# creates two separate arrays of corrected wavelengths/fluxes for ranges
# waveE1 = 5-15um and waveE2 = 15-38um
z = 0.00476
for i in waveO1:
    waveE1 = waveO1 / (1 + z)

for i in waveO2:
    waveE2 = waveO2 / (1 + z)

#break 3 correction: UNCOMMENT AS NEEDED
#for i in waveO3:
#    waveE3 = waveO3 / (1 + z)

# Flux Correection: SHORT->LONG
# corrects the break in flux from changing filters, calculated by hand
# SCALING UP FACTOR (SF): SF = f_2 / f_1 (f = flux and break wavelength)
for i in flux1:
    SF = 1.167838416
    fluxf1 = flux1 * SF

# SCALING DOWN FACTOR: UNCOMMENT AS NEEDED
#for i in flux2:
#    SF = 1.10898574
#    fluxf2 = flux2 * SF

# break3 corrections: UNCOMMENT AS NEEDED
#for i in flux2:
#    SF2 = 1.105475748
#    fluxf2 = flux2 * SF2

# Recombine wavelengths and fluxes for plotting:
# 1) 5-15um combination: BREAK1 COMBO

with open ('CASSIS\CASSIS_ACbreak1\cassis_ngc7314AC1.txt', 'w') as fh:
    for a, b, c in zip(waveE1, fluxf1, err1):
        print('%.7f  %.7f  %.7f ' % (a, b, c), file = fh)

# 2) 15-38um combination: BREAK2 COMBO
with open ('CASSIS\CASSIS_ACbreak2\cassis_ngc7314AC2.txt', 'w') as fh:
    for a, b, c in zip(waveE2, flux2, err2):
        print('%.7f  %.7f  %.7f ' % (a, b, c), file = fh)

# 3) 15-38um combination: BREAK3 COMBO -> UNCOMMENT AS NEEDED
#with open ('CASSIS\CASSIS_ACbreak3\cassis_ngc3783AC3.txt', 'w') as fh:
#    for a, b, c in zip(waveE3, flux3, err3):
#        print('%.7f  %.7f  %.7f ' % (a, b, c), file = fh)

# combine both created .txt files into one readable .txt file:
data = data2 = ''
with open('CASSIS\CASSIS_ACbreak1\cassis_ngc7314AC1.txt') as fp:
    data = fp.read()

with open('CASSIS\CASSIS_ACbreak2\cassis_ngc7314AC2.txt') as fp:
    data2 = fp.read()

# break3 combo -> UNCOMMENT AS NEEDED
#data3 = ''
#with open('CASSIS\CASSIS_ACbreak3\cassis_ngc3783AC3.txt') as fp:
#    data3 = fp.read()

# merge files;
data += ''
data += data2
# break 3 addition -> UNCOMMENT AS NEEDED
#data += data3

with open('CASSIS\CASSIS_FINAL\cassis_ngc7314.txt', 'w') as fp:
    fp.write(data)

# Re-read corrected wavelength/flux for cassis and plot:
waveC, fluxC, errC = np.loadtxt('CASSIS\CASSIS_FINAL\cassis_ngc7314.txt', unpack = True)

# plot cassis
#plt.subplot(221)
plt.plot(waveC, fluxC, color = 'C2', linewidth = 2)

# plot hawc
#plt.scatter(d, e, color = 'r', marker = 'x', linestyle = 'None')

# plot herschel
#plt.subplot(222)
#plt.scatter(h, i, color = 'm', marker = 'o', linestyle = 'None')

plt.xlabel('Wavelength[um]')
plt.ylabel('Flux[Jy]')
plt.title('NGC 7314 Sy 1.9/2')

#plt.loglog(wave, flux)

#plt.loglog(j,k)

#plt.savefig('ngc7314.png')
plt.show()