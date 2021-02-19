import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import astropy.io.ascii as ascii

#read in spitzer spectra
infile_spitzer = ascii.read('out_signatures.txt')
a = infile_spitzer['File']
b = infile_spitzer['Object']
type= infile_spitzer['Type']

#read in hi-res spectra
infile_hires = ascii.read('out_hires.txt')
c = infile_hires['HiRes']

#create subplots
fig,axes = plt.subplots(nrows=4,ncols=1,sharex=True,figsize=(4,8),gridspec_kw={'hspace':0})
plt.ylabel=('Flux')
plt.xlabel=('Wavelength [$\mu$m]')
for ii in range(len(infile_spitzer['File'])):

    wave,flux,err = np.loadtxt(a[ii],unpack=True)
    axes[ii].plot(wave,flux,color='black')
    axes[ii].set_ylabel('Flux [Jy]')
    axes[ii].text(5,0.75*max(flux),b[ii])
    axes[ii].text(5,0.60*max(flux),type[ii])
    axes[ii].xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    axes[ii].yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))

    wave2,flux2,err2=np.loadtxt(c[ii],unpack=True)
    axes[ii].plot(wave2,flux2,color='red')

fig.suptitle('Outlier AGN')
fig.savefig('hires_out_agn.pdf')
plt.show()
