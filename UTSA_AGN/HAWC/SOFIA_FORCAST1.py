# GRISM SPECTRA: BASIC INSPECTION AND ASSESSMENT

# INTRODUCTION:
# This recipe provides an overview and sample python code for plotting
# and assessing FORCAST grism data. The FORCAST observing modes are
# described in the SOFIA Observer’s Handbooks, available from the
# Proposing and Observing page on the SOFIA website; and the FORCAST
# data products are described in the FORCAST GO Data Handbook. Raw
# FORCAST data suffers from several instrumental artifacts. Nearly all
# of the artifacts are removed or corrected by the FORCAST pipeline,
# including: bad pixels; the Droop effect; non-linear pixel response;
# and the “jailbar” effect. In addition, the grism pipeline extracts the
# one-dimensional spectra, applies both wavelength and flux calibration,
# and corrects for telluric absorption. For point sources, an "optimal
# extraction" algorithm is used, while for extended (non-point-like)
# sources, a standard summation over a fixed aperture is used (see
# recipe on custom extractions for additional information). See the
# FORCAST GO Data Handbook for details regarding the artifacts, pipeline
# algorithms, and the flux-calibration process. Level 3 FORCAST grism
# data is written in FITS format as a 2-dimensional data array with a
# single header (for detailed description, see FORCAST GO Data Handbook).
# The data is in a 5 x N array, where N is the number of wavelength
# bins/samples, with the following convention:
#
# 1) The first line of the data array contains the wavelength (μm);
# 2) the second line contains the pipeline-reduced and flux calibrated
#       spectrum (Jy);
# 3) the third line contains the uncertainty (Jy);
# 4) the fourth contains the adopted atmospheric transmission spectrum;
# 5) and the fifth contains the instrumental response used in the flux
#   calibration process (Me-/s/Jy, where "Me-" is 10^6 electrons).
#
# Each FORCAST grism file contains the spectrum for a specific grism,
# which is indicated in the name of the file, e.g. "FORG111" would be a
# file containing data for the 11.1 μm grism. The pipeline does not
# combine specta from different grisms, this is left to the user.
#
# There are two data types for LEVEL 3 FORCAST grism data:
#
# 1) CAL files: flux calibrated results for each raw (LEVEL_1) data file.
# There can be multiple CAL files for a single AOR.
#
# 2) CMB files: if possible, multiple CAL files for a single AOR are
# coadded into a "combined" data file. Usually there is only one CMB file
# per AOR.
#
# Except for the provenance, the two datatypes are identical and hence the
# procedures below are equally valid for both. The "error" spectrum is for
# the statistical uncertainty only and does not reflect systematic
# uncertainties in the absolute flux, for example losses due to slit
# misalignment or poor image quality (see below). If absolute fluxes are
# required, it is highly recommended that additional imaging data are
# obtained at the same time in order to calibrate the grism spectra
# accurately.

# 1) VIEWING THE FITS HEADER
# First we will read a FITS file for one of the grism observations,
# display the full header, and then pick out some important keyword values
# to display in a "summary" table:
from astropy.io import fits

#first open the file...
print('1) Opening a FITS File and Displaying Info:\n')
g063_fits = fits.open('../forcast-sample-data/F0434_FO_GRI_0501381_FORG063_CMB_0228-0229.fits')
#...and get info
g063_fits.info()

# You can see from the FITS file info that the file contains a single
# primary HDU with 5 x 242 element array containing the calibrated data.
# Now display the full header:

#Now display FITS header in full:
print('1b) Display Full FITS Header using keywords:\n')
print(g063_fits[0].header)

# The header is quite large, so we might want to just pull out some
# keywords of interest:
print('1c) FITS Header Useful Bits:\n')
#And now pull out useful summary data from specific keywords:
keywords = ['ALTI_STA', 'ALTI_END', 'AOR_ID ', 'CHPAMP1 ', 'DATE-OBS',
            'MISSN-ID  ', 'NODAMP  ', 'OBJECT ', 'SKY_ANGL ', 'SKYMODE  ',
            'SLIT    ', 'ZA_START  ', 'ZA_END  ', 'TELVPA  ', 'NODSTYLE',
            'TOTINT    ']

#Loop over keywords to print.
for kywd in keywords:
        print(kywd, '=', g063_fits[0].header[kywd], '/',
              g063_fits[0].header.comments[kywd])

# This list of keywords provides a nice summary of the observation.
# There are a few different keywords in the header for integration
# time(s); here we show TOTINT because it reflects the total on-source
# integration time (factoring in chop-nod style). This is the value that
# the GO should use to compare to the FORCAST integration time estimator.
# It is also important to take a look at the HISTORY cards to see a
# summary of the data pipelining and especially the QA notes, which are
# usually at the end of the HISTORY block:
print('1d) Displays the History of Data Pipeline:\n')
#Display HISTORY block:
print(g063_fits[0].header['HISTORY'])

# So in this case, manual optimization of the telluric correction
# failed and the default atran model was used for the reduction. So the
# GO might want to examine the G063 data to see if there are issues with
# the telluric correction in the data.

# 2) SLIT ORIENTATION ON THE SKY
# For some observations, it is important to know the slit orientation
# and geometry on the sky. The slit width (arcsec), height (arcsec),
# and long-axis position angle (degrees, E of N) are logged in the
# keywords SLTW_ARC, SLTH_ARC, and SKY_ANGL, respectively:
print('2) Display Keywords Related to Slit Geometry:\n')
#Display keywords related to slit geometry:
slit_keywords = ['SLTW_ARC', 'SLTH_ARC', 'SKY_ANGL']
#Loop over keywords to print.
for kywd in slit_keywords:
        print(kywd, '=', g063_fits[0].header[kywd], '/',
              g063_fits[0].header.comments[kywd])
#
# SUPER IMPORTANT PROCESS
#
# We can use the slit information and sky position to create a DS9
# region file which can then be imported into DS9 and overplotted on any
# image of the object with a valid WCS. Specification for regions in
# DS9 can be found here: http://ds9.si.edu/doc/ref/region.html. Note
# that either the TELRA/DEC or OBSRA/DEC can be used for the sky
# position: TELRA/DEC is the sky position as reported by the telescope
# at the time of the observation whereas OBSRA/DEC is the sky position
# specified by the observer in the AOR. Due to various issues, the
# TELRA/DEC is not always accurate, especially for chop/nod observations.
# We recommend trying both and contacting your instrument scientist if
# there is a discrepancy:
print('2b) Creating a File to be opened in DS9:\n')
#First store the relevant values in local variables, for convenience.
object = g063_fits[0].header['OBJECT']
#converting from hours to degrees
ra = str(15.0*g063_fits[0].header['OBSRA'])
dec = str(g063_fits[0].header['OBSDEC'])
slitw_arc = str(g063_fits[0].header['SLTW_ARC'])
slith_arc = str(g063_fits[0].header['SLTH_ARC'])
sky_angl = str(g063_fits[0].header['SKY_ANGL'])

#now create a regions file and write the values above in according to
# region file spec:
file = open(object+'_slit.reg', 'w')
file.write('# Region file format: DS9 version 4.1\n')
file.write('global color = green dashlist = 8 3 width = 1 '
           'font = "helvetica 10 normal roman" select = 1 highlite = 1 '
           'dash = 0 fixed = 0 edit = 1 move = 1 delete = 1 include = 1 '
           'source = 1\n')
file.write('fk5\n')
file.write('box(' + ra + 'd,' + dec + 'd,' + slitw_arc + '\",' +
           slith_arc + '\",' + sky_angl + ') # text = {' + object + '}')
file.close()

# To use the region file, first load an image of your target into DS9.
# In this example, we will load a WISE image of AB Aur, using the
# downloaded app 'SOAImageDS9' and the downloaded image. Then click
# "Regions" -> "Load Region...", navigate to the local directory where
# your region file was saved, and click on the region file (in this
# example: "AB Aur (1) slit.reg"). Format should be "ds9"; the slit
# size and orientation should then be overplotted on the image:
#
# FOR MY PURPOSES: use this particular portion to overlay orientation
# slits over AGN axes. Download images of the AGN from the IRSA archive

# 3) READING THE FILE AND LOADING THE DATA INTO A TABLE STRUCTURE
#
# Now we'll read a FITS file for the grism of interest and load the data
# into a convenient set of table structures in Python. We'll start with
# the G227 (22.7 μm) grism data. The file contains a single HDU
# containing a 5 x 236 array. Now place the data portion into a separate
# table structure for convenience:
import numpy as np
from astropy.table import Column
from astropy.table import Table
from astropy import constants as const
from astropy import units as u

print('3) Read a FITS file and Load the data into a Table:\n')
# Define a function for loading fits data into tables with units
def loadFORCASTGrismData(filename):
    # Now open fits file for the sample data...
    data_fits = fits.open(filename)

    # ... read the data portion of the file into a separate array:
    data_tmp = data_fits[0].data

    # ... load into table for convenience:
    data_table = Table([data_tmp[0], data_tmp[1], data_tmp[2],
                        data_tmp[3], data_tmp[4]],
                       names = ('wavelength', 'flux', 'error',
                                'telluric', 'response'), masked = True,
                       meta = {'name': 'Data Table'})

    # ...and assign units:
    data_table['wavelength'].unit = 'micron'
    data_table['flux'].unit = 'Jy'
    data_table['error'].unit = 'Jy'
    # response is (Me-/s)/Jy
    data_table['response'].unit = '1/(Jy*s)'
    # ...and mask NaNs in the flux.
    data_table['flux'].mask = np.isnan(data_table['flux'])
    # return both the FITS structure and the table.
    return data_fits, data_table

# Now read the file for the G227 data:
g227_fits, g227 = loadFORCASTGrismData('../forcast-sample-data/F0434_FO_GRI_0501383_FORG227_CMB_0230-0233.fits')

print(g227.info)

# 4) PLOTTING THE SPECTRA
#
# Start by plotting flux (with errors) and S/N as a function of
# wavelength:
import matplotlib.pyplot as plt

print('4) Plotting Spectra w/ Errors & S/N vs wavelength:\n')
plt.figure(figsize = (15, 10))
plt.errorbar(g227['wavelength'], g227['flux'], yerr = g227['error'],
             fmt = 'o', label = 'Flux')
# This plot will not show if the Telluric Plot is available
plt.title("4A) G227: Wavelenth vs Flux")
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('F$_v$ (Jy)')
plt.legend(loc = 'upper left')
plt.figure(figsize = (15, 10))
plt.step(g227['wavelength'], g227['flux']/g227['error'], label = 'SNR',
         where = 'mid')

plt.title("4B) G227: Wavelength vs S/N")
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('SNR')
plt.legend(loc = 'upper left')
plt.show()

# In order to assess whether any particular features in the spectrum are
# real or not, it is important to compare the flux to the telluric
# transmission. We'll write a short function for comparing flux
# spectrum to telluric transmission from the same data table:
print('4c) Plotting double FORCAST Spectral and Comparing:\n')
#Define a function to plot two spectra from the same FORCAST data table:
def compareTelluric(table):
    fig, ax1 = plt.subplots(figsize = (15, 10))

    ax1.set_xlabel('Wavelength ($\mu$m)')
    ax1.set_ylabel('F$_v$ (Jy)')
    ax1.errorbar(table['wavelength'], table['flux'], yerr = table['error'],
                 fmt = 'o', label = 'Flux')
    plt.legend(loc='lower left')

    color = 'tab:red'
    # instantiate a second axes that shares the same x-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('Transmission', color = color)
    ax2.plot(table['wavelength'], table['telluric'], color = color,
             linestyle = ':', label = 'Telluric Transmission')
    ax2.tick_params(axis = 'y', labelcolor = color)
    ax2.set_ylim(bottom = 0.0)
    #plt.title(table)
    plt.legend(loc = 'upper left')
    plt.title('4C) Multi-Plot -> see Print notes')
    plt.show()

compareTelluric(g227)
print('1st-Comparison:\n', compareTelluric(g227))

# It is clear form this plot that most of the "features" in the spectrum
# are in fact residuals left over from the telluric correction. It is
# sometimes possible to improve the correction by careful modelling of
# the atmospheric lines, but some residuals are almost always present.
# For very low-S/N data, it is useful to look at both the telluric
# transmission and the response. Here we'll take a look at the G329
# grism data as an example:
print('4d) Telluris Transmission and Response for low-S/N:\n')
#Read in data for G329 grism using the function we defined earlier:
g329_fits,g329 = loadFORCASTGrismData('../forcast-sample-data/F0434_FO_GRI_0501384_FORG329_CMB_0239-0246.fits')

#now plot flux and telluric transmission together
print('2nd-Comparison:\n', compareTelluric(g329))
#and now plot flux and response  together
fig, ax1 = plt.subplots(figsize = (15, 10))

ax1.set_xlabel('Wavelength ($\mu$m)')
ax1.set_ylabel('F$_v$ (Jy)')
ax1.errorbar(g329['wavelength'], g329['flux'], yerr = g329['error'],
             fmt = 'o', label = 'Flux')
plt.legend(loc = 'lower left')
color = 'tab:red'
# instantiate a second axes that shares the same x-axis
ax2 = ax1.twinx()
ax2.set_ylabel('Response', color = color)
ax2.plot(g329['wavelength'], g329['response'], color = color,
         linestyle = ':', label = 'Response')
ax2.tick_params(axis = 'y', labelcolor = color)
plt.legend(loc = 'upper left')
plt.title('4D) Multi Plot -> see Print notes')
plt.show()

# Due to the strong telluric absorption coupled with relatively low
# response, the flux values for λ > 35 μm should be viewed with some
# caution.

# 5) MASKING REGIONS WITH STRONG TELLURIC FEATURES
#
# Regions with particularly strong telluric absorption generally suffer
# from very poor correction and should be masked out. For example, the
# G111 grism generally suffers from poor correction at 9.6 μm due to the
# very deep telluric O3 band there:
print('5) Un-masked Telluric O3 band plot:\n')
#Read in data for G111 grism using the function we defined earlier:
g111_fits,g111 = loadFORCASTGrismData('../forcast-sample-data/F0434_FO_GRI_0501382_FORG111_CMB_0234-0238.fits')

#and use our plotting function from earlier to compare flux to telluric
# transmission:
compareTelluric(g111)
print(compareTelluric(g111))

# Now mask out all flux values for which the transmission is less than
# some threshold (0.7 in this example) and re-plot. Note that matplotlib
# automatically applies the mask to the data when plotting
print('5b) Masking out all flux values for trasmission < 0.7 plot:\n')
#Pick masking threshold for telluric transmission:
thresh = 0.7

#and generate masks for flux and error columns.
g111['flux'].mask  = g111['telluric'] < thresh
g111['error'].mask = g111['telluric'] < thresh

#Re-plot; masks are applied automatically.
compareTelluric(g111)
print(compareTelluric(g111))

# 6) PLOTTING THE COMPLETE SED
#
# Now we can load the last grism data file (G063) and plot the full SED
# from 5 -- 40 mic.:
print('6) Plotting the compete SED of all 4 sources:\n')
#And finally, read in the G063 data:

g063_fits,g063 = loadFORCASTGrismData('../forcast-sample-data/F0434_FO_GRI_0501381_FORG063_CMB_0228-0229.fits')

#Now plot all grism data together as a line plot, semi-log
plt.figure(figsize = (15, 10))
plt.semilogx(g063['wavelength'], g063['flux'], label = 'G063')
plt.semilogx(g111['wavelength'], g111['flux'], label = 'G111')
plt.semilogx(g227['wavelength'], g227['flux'], label = 'G227')
plt.semilogx(g329['wavelength'], g329['flux'], label = 'G329')

plt.title("6) AB Aur: All 4 Sources")
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('F$_v$ (Jy)')
plt.legend()
plt.show()

# (Notice that Python automatically uses the masked data for G111.) The
# first thing we notice is that there appears to be a flux break or
# discrepancy between the G227 and G329 grisms, which is most likely
# due to additional slit losses for the G329 data.

# 7) CALIBRATION USING INDEPENDANT PHOTOMETRY
#
# The absolute flux calibration accuracy of the grism data can be
# uncertain by >10% due to systematic uncertainties at the time of
# observation (e.g. poor atmospheric conditions, bad seeing, or slit
# mis-alignment). In this section, we show how to use independant
# FORCAST photometric observations to "re-calibrate" the grism data
# (assuming the imaging observations are considered accurate). The most
# common reason for systematic error in the grism flux data is
# mis-alignment of the slit during the observation. Since the width of
# the FORCAST slits are smaller than the SOFIA PSF in general, the slit
# only "samples" a fraction of the flux from the source. The FORCAST
# pipeline accounts for this "loss" automatically for each grism using
# a default PSF model. However, if the target is not centered on the
# slit during the observation, or if the image quality is poor, there
# can be an additional loss that is not accounted for in the pipeline
# processing. These additional "slit losses" can be assessed and "fixed"
# using trusted absolute photometry from the same spectral region. For
# the sample data, we also have photometric fluxes for each grism band
# pass which we use to compare to the grism spectra below:
print('7) All 4 Sources w/ Photometric Values Overlaid:\n')
#Create a table of values for the photometry
filters = ['F064', 'F111', 'F197', 'F315']
#microns
waves = [6.34774, 11.0888, 19.6993, 31.3615]
#filter bandpass (FWHM, microns)
delwaves = [0.14, 0.95, 5.5, 5.7]
#Jy
fluxes = [12.854, 26.778, 47.129, 77.264]
#absolute uncertainty (ratio)
relerrs = [0.063705, 0.073197, 0.063439, 0.089142]

phot = Table([filters, waves, fluxes, relerrs, delwaves],
                  names = ('filter', 'wave', 'flux', 'relerr', 'delwave'),
                  meta = {'name': 'Photometry Data Table'})

phot['wave'].unit = 'micron'
phot['delwave'].unit = 'micron'
phot['flux'].unit = 'Jy'

#Now plot all grism data together as a line plot, using semi-log (x-axis)
fp = 'tab:blue'
plt.figure(figsize = (15, 10))
plt.semilogx(g063['wavelength'], g063['flux'], label = 'G063')
plt.semilogx(g111['wavelength'], g111['flux'], label = 'G111')
plt.semilogx(g227['wavelength'], g227['flux'], label = 'G227')
plt.semilogx(g329['wavelength'], g329['flux'], label = 'G329')

plt.errorbar(phot['wave'], phot['flux'], xerr = 0.5 * phot['delwave'],
             yerr = phot['flux'] * phot['relerr'], fmt = 'o',
             color = fp, label = 'Photometry')

plt.title("7) AB Aur: All 4 Sources w/ Photometry Overlaid")
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('F$_v$ (Jy)')
plt.legend()
plt.show()

# The grism spectra match up very well with the photometry except for
# G329 which appears too low, most likely due to slit misalignment
# during the observation. Assuming the losses are independant of
# wavelength, the grism data can be shifted up to match the photometric
# value by simply averaging the grism data (weighted by the filter
# transmission), calculating the offset between the average grism flux
# and the photometric flux, and applying that offset to the grism flux
# data. Filter transmission curves can be found in the FORCAST
# Observers Handbook, which is available from the SOFIA Proposal
# Resources page:
print('7b) Matching Photometric and GRISM Data:\n')
import numpy as np
import numpy.ma as ma

#read in filter transmission file from SOFIA website.
f315 = Table.read('https://www.sofia.usra.edu/sites/default/files/31pt5mu.txt',
                  format = 'ascii', names = ('wave', 'trans'))


#re-bin transmission onto our grism data wavelength grid and save in
# new array.
f315_wts = np.interp(g329['wavelength'], f315['wave'], f315['trans'])

#and now plot flux and response  together
fig, ax1 = plt.subplots(figsize=(15, 10))

ax1.set_xlabel('Wavelength ($\mu$m)')
ax1.set_ylabel('F$_v$ (Jy)')
ax1.errorbar(g329['wavelength'], g329['flux'], yerr = g329['error'],
             fmt = 'o', label = 'Flux')
plt.legend(loc = 'lower left')

color = 'tab:red'
# instantiate a second axes that shares the same x-axis
ax2 = ax1.twinx()
ax2.set_ylabel('F315 Filter Transmission', color = color)
ax2.plot(g329['wavelength'], f315_wts, color = color, linestyle = ':',
         label = 'Filter Transmission')
ax2.tick_params(axis = 'y', labelcolor = color)
plt.legend(loc = 'upper left')

#now compute weighted average using the filter transmission for the
# weights.
g329_av = ma.average(g329['flux'], weights = f315_wts)

#and shift the spectrum up by the difference
g329_f_adj = g329['flux'] + (phot['flux'][3] - g329_av)

#and add to the data table
g329['flux_adj'] = g329_f_adj
g329['flux_adj'].unit = 'Jy'

# Here we've added an extra column for the new "adjusted" flux spectrum.
# Now re-plot with the adjusted values for G329:
fp = 'tab:blue'
plt.figure(figsize = (15, 10))

# plot all grism settings
plt.semilogx(g063['wavelength'], g063['flux'], label = 'G063')
plt.semilogx(g111['wavelength'], g111['flux'], label = 'G111')
plt.semilogx(g227['wavelength'], g227['flux'], label = 'G227')

# plot adjusted spectrum for g329
# Adjusted spectrum:
plt.semilogx(g329['wavelength'], g329['flux_adj'], label = 'G329 Adj')
plt.errorbar(phot['wave'], phot['flux'], xerr = 0.5 * phot['delwave'],
             yerr = phot['flux'] * phot['relerr'], fmt = 'o',
             color = fp, label = 'FORCAST Phot (2017)')

plt.title("7b) AB Aur: Photometric w/ GRISM Correction")
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('F$_v$ (Jy)')
plt.legend()
plt.show()

# The adjusted G329 spectrum now matches the G227 spectrum very nicely.

# 8) CONVERTING TO λF_λ
#
# We often want to work with λF_λ in order to assess where the object is
# emitting the most power (per logade or decade). AstroPy quantities
# and constants makes this easy. Remember that (λF_λ = vFv = c*Fv/λ):
print('8) Working w/ λF_λ to access where the object is emitting the '
      'most power (per logade or decade):\n')
#Calculate lam f_lam (= nu F_nu = c * F_nu / lambda), convert to W/m^2,
# and add to table with mask if it exists.
# calculate:
lamflam=const.c * g063['flux'] / g063['wavelength']
# and convert to W/m^2:
lamflam=lamflam.to(u.W / (u.m * u.m))
# ...put in table...
g063['lamflam'] = lamflam
# calculate:
lamflam=const.c * g111['flux'] / g111['wavelength']
# and convert to W/m^2:
lamflam=lamflam.to(u.W / (u.m * u.m))
# ...put in table...
g111['lamflam'] = lamflam
# ...and dont forget to mask lamflam as well:
g111['lamflam'].mask = g111['flux'].mask
# calculate:
lamflam=const.c * g227['flux'] / g227['wavelength']
# and convert to W/m^2:
lamflam=lamflam.to(u.W / (u.m * u.m))
# ...put in table...
g227['lamflam'] = lamflam
# calculate using adjusted flux:
lamflam=const.c * g329['flux_adj'] / g329['wavelength']
# and convert to W/m^2:
lamflam=lamflam.to(u.W / (u.m * u.m))
# ...put in table...
g329['lamflam'] = lamflam

g111.info
print(g111.info)

# We've retained MKS units, but it is easy to convert to CGS from here
# if desired. And now re-plot with λF_λ:
fp = 'tab:blue'
plt.figure(figsize = (15, 10))

# plot all grism settings
plt.semilogx(g063['wavelength'], g063['lamflam'], label = 'G063')
plt.semilogx(g111['wavelength'], g111['lamflam'], label = 'G111')
plt.semilogx(g227['wavelength'], g227['lamflam'], label = 'G227')

# plot adjusted spectrum for g329
# Adjusted spectrum:
plt.semilogx(g329['wavelength'], g329['lamflam'], label = 'G329 Adj')

plt.title("8) AB Aur: Working w/ λF_λ to access where the object is "
          "emitting the most power (per logade or decade)")
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('$\lambda$F$_{\lambda}$ (W/m^2)')
plt.legend()
plt.show()