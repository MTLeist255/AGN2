# Read in 9/18 um spectra, add a straight line fit to data, then plot

import numpy as np
import matplotlib.pyplot as plt

from astropy.modeling import models
from astropy import units as u

from specutils.spectra import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum

#9 um data -> redshift corrected spectra broken down between 6-12 um
wave9, flux9, errC = np.loadtxt('CASSIS/CASSIS_9um/cassis_3c273_9um.txt', unpack = True)
# 18 um data -> redshift corrected spectra broken down between 15-22 um
wave18, flux18, errC = np.loadtxt('CASSIS/CASSIS_18um/cassis_3c273_18um.txt', unpack = True)

# Add straight line continuum: 9/18 um
# 9um: bounds 7-11 um
x_9 = [7.0263930, 10.9744104] # x = xstart -> xfinish
y_9 = [ 0.2879851, 0.4079096] # y = ystart -> yfinish
# 18um: bounds 16-21 um
x_18 = [16.0109949, 21.9018499]
y_18 = [0.5584640, 0.7031290]

plt.figure(figsize=(12,6))
# 9um plot
plt.subplot(1,2,1)
plt.plot(wave9, flux9)
plt.plot(x_9, y_9)
plt.title('test 2')
plt.xlabel('wavelength (um)')
plt.ylabel('Flux (mJy)')
# 18um plot
plt.subplot(1,2,2)
plt.plot(wave18, flux18)
plt.plot(x_18, y_18)
plt.title('test 2')
plt.xlabel('wavelength (um)')
plt.ylabel('Flux (mJy)')

#plt.savefig('Current_AGN_fits/silicate/9_18_continuum/test4.png')
plt.show()