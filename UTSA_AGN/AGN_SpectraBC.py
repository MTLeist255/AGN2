# Plotting the UN-CORRECTED wavelength and flux for spectral cassis,
# herschel and hawc data for several various AGN in our sample set.
# To run this program with available data, just un-comment which
# section you wish to use. Finally, once you are comfortable with the
# plot, uncomment the 'save.fig' function at the end of the code. Once
# the fig is saved, moved the completed plot to the correct folder for
# filing and reviewing.

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker


plt.rcParams['figure.figsize'] = (10, 8)

# Cassis data: awaveO = wavelength observed, flux = flux, err = error
waveC, fluxC, errC = np.loadtxt('CASSIS/CASSIS_BC/BC_combo/cassis_circinusHiRES_combo.txt', unpack = True)
print(waveC)

# input hawc data
#waveH, fluxH, errH, g = np.loadtxt('HAWC/hawcngc4388.txt', unpack = True)

# input herschel data
#waveHer, fluxHer = np.loadtxt('HERSHEL/herschelngc4388.txt', unpack = True)

# WAVELENGTH CORRECTION -> UNCOMMENT AS NEEDED
z = 0.0010332
waveE = waveC / (1 + z)

#Process to read-in, correct wavelength, and create final corrected data -> UNCOMMENT AS NEEDED
data = ''
with open ('CASSIS/CASSIS_BC/BC_combo/cassis_circinusHiRES_BC.txt', 'w') as fh:
   for a, b, c in zip(waveE, fluxC, errC):
       print('%.7f  %.7f  %.7f ' % (a, b, c), file = fh)

with open('CASSIS/CASSIS_BC/BC_combo/cassis_circinusHiRES_BC.txt') as fp:
   data = fp.read()

with open('CASSIS/CASSIS_FINAL/cassis_circinusHiRES.txt', 'w') as fp:
   fp.write(data)

waveF, fluxF, errF = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_circinusHiRES.txt', unpack = True)

# plot cassis
#plt.subplot(221)

tick_spacing = 1
fig, ax = plt.subplots(1, 1)
# Corrected wavelength -> UNCOMMENT AS NEEDED
ax.plot(waveF, fluxF)
#ax.plot(waveC, fluxC)

# plot hawc
#plt.scatter(d, e, color = 'r', marker = 'x', linestyle = 'None')

# plot herschel
#plt.subplot(222)
#plt.scatter(h, i, color = 'm', marker = 'o', linestyle = 'None')

ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
plt.xlabel('Wavelength[um]')
plt.ylabel('Flux[Jy]')
plt.title('Circinus Sy 2')
#plt.loglog(wave, flux)

#plt.loglog(j,k)

#plt.savefig('ngc7314_before.png')
plt.show()