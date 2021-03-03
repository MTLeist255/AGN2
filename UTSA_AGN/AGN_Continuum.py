# Fit a generic continuum to the redhsift corrected spectra

import numpy as np
import matplotlib.pyplot as plt

from astropy.modeling import models
from astropy import units as u

from specutils.spectra import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum

waveC, fluxC, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_mrk231.txt', unpack = True)

spectrum = Spectrum1D(flux=fluxC*u.Jy, spectral_axis=waveC*u.um) # creates bounds to fit continuum

g1_fit = fit_generic_continuum(spectrum) # adds fit

y_continuum_fitted = g1_fit(waveC*u.um) # adds fit to flux data

plt.plot(waveC, fluxC) # plot spectrum
plt.plot(waveC, y_continuum_fitted) # plot continuum
plt.title('Continuum Fitting')
plt.grid(True)
#plt.savefig('Current_AGN_fits/silicate/continuum/test.png')
plt.show()