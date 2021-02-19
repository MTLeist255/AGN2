# HAWC+ DATA RECIPE / 30 DORADUS DATA

# In this jupyter recipe, we will explore the HAWC+ data cube and describe some
# of the basic analysis techniques involving imaging polarimetry data.

# INTRODUCTION

# The Strategic Director's Discretionary Time (S-DDT) for SOFIA is aimed at
# providing the astronomical community with data sets of high scientific
# interest over a broadrange of potential research topics without any
# proprietary period. These observingsessions allow the general user community
# access to high-level data products that aremeant not only for general
# understanding of SOFIA data and its packaging but also for inclusion in
# published scientific work. The S-DDT target have been selected on a
# non-interference basis with existing programs and in terms of SOFIA
# flight planning. The 76_0001 program, "Community Science: HAWC+ Polarimetry
# of 30 Dor," was designed and scheduled to provide the community with SOFIA
# polarimetry data of an important and relatively bright source. The observing
# strategy also provided significantly increased scheduling efficiency for the
# OC6I (HAWC+) flights in July 2018. The west-bound observing legs for 30
# Doradus allowed a larger fraction of the highest ranked Cycle 6 targets,
# predominantly in the inner Galaxy, to be scheduled and flown. To enhance the
# scientific exploitation of these data products, we present here an overview
# of the observations, visualizations of the data, and preliminary analysis of
# their quality. The purpose of this document is to provide a template for a
# jupyter notebook session, or even just some code snippets, that will aid in
# analysis of SOFIA/HAWC+ data. In this notebook, we focus on imaging
# polarimetry data; however, some of the techniques mentioned here would also
# apply to HAWC+ total intensity imaging as well. For more information on HAWC+
# observing modes, visit our observer's handbook at
# https://www.sofia.usra.edu/science/instruments/hawc.

# 1) DATA STRUCTURE
# For this analysis, we require the standard numpy/scipy/matplotlib stack as
# well the astropy and aplpy modules. With just a few lines of code, we can
# explore the HAWC+ fits data cubes and plot the images:

from astropy.io import fits

print('1) HAWC+ Data Cubes Info:\n')
efile = 'sofia_data/F0485_HA_POL_7600019_HAWEHWPE_PMP_055-075.fits'
dfile = 'sofia_data/F0481_HA_POL_7600012_HAWDHWPD_PMP_050-083.fits'
cfile = 'sofia_data/F0484_HA_POL_7600018_HAWCHWPC_PMP_022-114.fits'


afile = 'sofia_data/F0485_HA_POL_76000110_HAWAHWPA_PMP_043-052.fits'
hawc = fits.open(afile)
hawc.info()

# 2) STOKES I
# Stokes I---the zeroth extension in the fits file---represents the total
# intensity of the image, let us go ahead and plot this extension:
import matplotlib.pyplot as plt
from aplpy import FITSFigure

print('2) Plot the Zeroth Extension FITS File (or Total Intensity of the '
      'Image):\n')
# set colormap for all plots
cmap = 'rainbow'
# or hawc[0]. Note the extension is from the hawc.info() table above
stokes_i = hawc['STOKES I']
# load HDU into aplpy figure
axs = FITSFigure(stokes_i)
# display the data with WCS projection and chosen colormap
axs.show_colorscale(cmap = cmap)

# FORMATTING
axs.tick_labels.set_font(size = 'small')

# Add colorbar
axs.add_colorbar()
axs.colorbar.set_axis_label_text('Flux (Jy/pix)')
plt.title('2) Stokes I Plot: Total Intensity of the Image')
plt.show()

# 3) STOKES Q AND U
# Similarly, we can plot the Stokes Q and Stokes U images: Stokes Q and U
# represents the unpolarized intensity of the measurement:
print('3) Plot Stokes Q and U: Unpolarized Intensities')
stokes_q = hawc['STOKES Q']
stokes_u = hawc['STOKES U']

# generate FITSFigure as subplot to have two axes together
axq = FITSFigure(stokes_q, subplot = (1, 2, 1))
# show Q
axq.show_colorscale(cmap = cmap)
axu = FITSFigure(stokes_u, subplot = (1, 2, 2), figure = plt.gcf())
# show U
axu.show_colorscale(cmap = cmap)

# FORMATTING
axq.set_title('3) Stokes Q')
axu.set_title('3) Stokes U')
axu.axis_labels.set_yposition('right')
axu.tick_labels.set_yposition('right')
axq.tick_labels.set_font(size = 'small')
axq.axis_labels.set_font(size = 'small')
axu.tick_labels.set_font(size = 'small')
axu.axis_labels.set_font(size = 'small')
plt.show()

# We can additionally plot the associated error maps for each extension:
print('3b) Associated Error Maps for Stokes Q:\n')
stokes_q = hawc['STOKES Q']
error_q = hawc['ERROR Q']
# generate FITSFigure as subplot to have two axes together
axq = FITSFigure(stokes_q, subplot = (1, 2, 1))
# show Q
axq.show_colorscale(cmap = cmap)

axe = FITSFigure(error_q, subplot = (1, 2, 2), figure = plt.gcf())
# show error
axe.show_colorscale(cmap = cmap)

# FORMATTING
axq.set_title('3b) Stokes Q')
axe.set_title('3b) Error Q')
# hide axis/tick labels
axq.axis_labels.hide()
axe.axis_labels.hide()
axq.tick_labels.hide()
axe.tick_labels.hide()
plt.show()

print('3c) Associated Error Maps for Stokes U:\n')
stokes_u = hawc['STOKES U']
error_u = hawc['ERROR U']
# generate FITSFigure as subplot to have two axes together
axq = FITSFigure(stokes_u, subplot = (1, 2, 1))
# show Q
axq.show_colorscale(cmap = cmap)

axe = FITSFigure(error_u, subplot = (1, 2, 2), figure = plt.gcf())
# show error
axe.show_colorscale(cmap = cmap)

# FORMATTING
axq.set_title('3c) Stokes U')
axe.set_title('3c) Error U')
# hide axis/tick labels
axq.axis_labels.hide()
axe.axis_labels.hide()
axq.tick_labels.hide()
axe.tick_labels.hide()
plt.show()

# 4) POLARIZED INTENSITY I_p
# Level 4 HAWC+ additionally provides extensions with the linear polarization
# percentage (p), angle (θ), and their associated errors (σ). Percent
# polarization (p) and error (σ_p) are calculated as:
#
# p = 100 * ((Q / I)^2 + (U / I)^2)^1/2
#
# σ_p = (100/I)*((1 /(Q^2+U^2)((Q*σ_Q)^2+(U*σ_U)^2 + 2*q*u*σ_QU)+[(Q/I)^2+(U/I)^2]*σ^2_I-2(Q/I)*σ_QI-2(U/I)*σ_UI)^1/2
#
# Note that p here represents the percent polarization as opposed to the more
# typical convention for p as the fractional polarization.Maps of these data
# are found in extensions 7 (PERCENT POL) and 9 (ERROR PERCENT POL). Polarized
# intensity, Ip, can then be calculated as
#
# I_p = (I × p) / 100
#
# which is included in extension 13 (POL FLUX). Also included is the debiased
# polarization percentage (p′) calculated as:
#
# p′ = (p^2−σ^2_p)^1/2
#
# found in extension 8 (DEBIASED PERCENT POL). We similarly define the
# debiased polarized intensity as
#
# I_p′ = (I × p′) / 100
#
# which is included in extension 15 (DEBIASED POL FLUX). These values are
# also included in table form in extension 17 (POL DATA):
print('4) Intensity (I) vs Polaraized Intensity (I_p):\n')
stokes_ip = hawc['DEBIASED POL FLUX']

axi = FITSFigure(stokes_i, subplot = (1, 2, 1))
# show I
axi.show_colorscale(cmap = cmap)

axp = FITSFigure(stokes_ip, subplot = (1, 2, 2), figure = plt.gcf())
# show Ip
axp.show_colorscale(cmap = cmap)

# FORMATTING
axi.set_title('4) Intensity' r'$I$')
axp.set_title('Polarized Intensity  'r'$I_{p^\prime}$')
axp.axis_labels.set_yposition('right')
axp.tick_labels.set_yposition('right')
axi.tick_labels.set_font(size = 'small')
axi.axis_labels.set_font(size = 'small')
axp.tick_labels.set_font(size = 'small')
axp.axis_labels.set_font(size = 'small')
plt.show()

# 5) PLOTTING POLARIZATION VECTORS
# From the Q and U maps, the polarization angle θ is calculated in the
# standard way:
#
# θ = (90/π)tan^-1(U/Q)
#
# σ_θ=(90/π(Q^2+U^2))((Qσ_Q)^2+(Uσ_U)^2-2QUσ_QU)^1/2
#
# The angle map is stored in extension 10 (POL ANGLE) in degrees, with its
# error in extension 12 (ERROR POL ANGLE). As part of the HAWC+ reduction
# pipeline, θ is corrected for the vertical position angle of the instrument
# on the sky, the angle of the HWP plate, as well as an offset angle that is
# calibrated to each filter configuration. θ=0∘ corresponds to the North-South
# direction, θ=90∘ corresponds to the East-West direction, and positive values
# follow counterclockwise rotation. We also provide the PA of polarization
# rotated by 90∘, θ_90, in extension 11 (ROTATED POL ANGLE). This PA of
# polarization needs to be used with caution. If the measured polarization is
# dominated by magnetically-aligned dust grains, then the PA of polarization,
# θ, can be rotated by 90∘ to visualize the magnetic field morphology. For
# more details, see Hildebrand et al. 2000; Andersson et al. 2015. We can now
# use the p′ and θ_90 maps to plot the polarization vectors. First, however,
# let us make a quality cut. Rather than defining a σ cut on the polarization
# vectors themselves, it is more useful to define a signal-to-noise cut on
# total intensity, I, the measured quantity. Starting with the propagated
# error on the polarization fraction:
#
# σ_p = (100/I)*((1 /(Q^2+U^2)((Q*σ_Q)^2+(U*σ_U)^2 + 2*q*u*σ_QU)+[(Q/I)^2+(U/I)^2]*σ^2)^1/2
#
# Let's assume the errors in Q, U, and I are comparable so that there are no
# covariant (cross) terms in the error expansion. Therefore:
#
# σ_Q = σ_U = σ_Q,U
# σ_QI = σ_QU = σ_UI = 0
#
# σ_p = (1/I)*(σ^2_Q,U + σ^2_I*p^2)^1/2
#
# If we assume that p is relatively small (e.g. the source is not highly
# polarized), and that the errors in I are small, then the second term
# (σ^2_I*p^2) is negligible:
#
# σ_p = (σ_Q,U/I)
#
# By design, the HAWC+ optics split the incident radiation into two orthogonal
# linear polarization states that are measured with two independent detector
# arrays. The total intensity, Stokes I, is recovered by linearly adding both
# polarization states. If the data is taken at four equally-spaced HWP angles,
# and assuming 100% efficiency of the instrument, then the uncertainty in I is
# related to the uncertanties in Q and U:
#
# σ_Q ~ σ_U ~ (2)^1/2 * σ_I
#
# This simplifies our error on p to:
#
# σ_p ~ [(2)^1/2]/(S/N)_I
#
# if we desire an accuracy of σp ∼ 0.5%, we require a S/N in total intensity
# I of ∼283. This S/N cut in I is very conservative. In the Level 4 HAWC+
# data, the last extension, FINAL POL DATA, contains a table of values similar
# to POL DATA, with somewhat less restrictive quality cuts applied. This
# extension includes vectors satisfying the following three criteria:
#
# 1) (S/NI) > 200
# 2) (S/Np) > 3
# 3) p<50%
#
# Since we include maps of all measurable polarization information with the
# full data set, you are free to decide on any quality cuts that satisfy your
# scientific needs. In this next panel, we include a single quality cut on
# S/N > 100, by performing the following steps:
#
# 1) use the Stokes I image as the background for the vector plot
# 2) perform a quality cut on Stokes I/σI>100 to make a mask
# 3) mask out low S/N vectors
# 4) plot remaining polarization vectors
# 5) add contours to better visualize changes in flux across the map
from astropy.io import fits
import numpy as np
from aplpy import FITSFigure

print('5) Single quality cut on S/N > 100:\n')
p = hawc['DEBIASED PERCENT POL']
def make_polmap(filename, title = None, figure = None, subplot = (1, 1, 1)):
      hawc = fits.open(filename)
      # %
      p = hawc['DEBIASED PERCENT POL']
      # deg
      theta = hawc['ROTATED POL ANGLE']
      # I
      stokes_i = hawc['STOKES I']
      # error I
      error_i = hawc['ERROR I']

      # 1. plot Stokes I
      # convert from Jy/pix to Jy/sq. arcsec
      # map scale in arcsec/pix.CDELT2 always in deg
      pxscale = stokes_i.header['CDELT2'] * 3600
      stokes_i.data /= pxscale ** 2
      error_i.data /= pxscale ** 2

      fig = FITSFigure(stokes_i, figure = figure, subplot = subplot)

      # 2. perform S/N cuts on I/sigma_I, and p/sigma_p
      i_err_lim = 100
      mask = np.where(stokes_i.data / error_i.data < i_err_lim)

      # 3. mask out low S/N vectors by setting masked indices to NaN
      p.data[mask] = np.nan

      # 4. plot vectors
      # 1pix = scalevec * 1% pol -> scale vectors to make it easier to see
      scalevec = 0.4
      # step size = display every 'step' vectors
      fig.show_vectors(p, theta, scale=scalevec, step=2)
      # step size of 2 is effectively Nyquist sampling -> close to the beam
      # size

      # 5. plot contours
      ncontours = 30
      fig.show_contour(stokes_i, cmap = cmap, levels = ncontours, filled = True,
                       smooth = 1, kernel = 'box')
      fig.show_contour(stokes_i, colors = 'gray', levels = ncontours,
                       smooth = 1, kernel = 'box', linewidths = 0.3)

      # Show image
      fig.show_colorscale(cmap = cmap)

      # If title, set it
      if title:
            fig.set_title(title)

      # Add colorbar
      fig.add_colorbar()
      fig.colorbar.set_axis_label_text('Flux (Jy/arcsec$^2$)')

      # Add beam indicator
      fig.add_beam(facecolor = 'red', edgecolor = 'black', linewidth = 2,
                   pad = 1, corner = 'bottom left')
      fig.add_label(0.02, 0.02, 'Beam FWHM', horizontalalignment = 'left',
                    weight = 'bold', relative = True, size = 'small')

      # Add vector scale: polarization vectors are displayed such that
      # 'scalevec' * 1% pol is 1 pix long must translate pixel size to
      # angular degrees since the 'add_scalebar' function assumes a physical
      # scale
      vectscale = scalevec * pxscale / 3600
      fig.add_scalebar(5 * vectscale, "p = 5%", corner = 'top right',
                       frame = True)

      # FORMATTING
      fig.tick_labels.set_font(size = 'small')
      fig.axis_labels.set_font(size = 'small')
      plt.title('5) Single quality cut on S/N > 100')
      plt.show()

      return stokes_i, p, mask, fig

make_polmap(afile)

# 6) PLOTTING POLARIZATION FRACTION
# We can also plot the polarization fraction p to better visualize the
# structure of 30 Doradus. We plot the same contours from total intensity I
# in the background:
print('6) Polarization Fraction of 30 Doradus w/ Intensity in the '
      'background:\n')
# could not get to match the given code (?)
fig = FITSFigure(p)

# Show image
fig.show_colorscale(cmap = cmap)

# Plot contours
ncontours = 30
fig.show_contour(stokes_i, colors = 'gray', levels = ncontours, smooth = 1,
                 kernel = 'box', linewidths = 0.3)

# Add colorbar
fig.add_colorbar()
fig.colorbar.set_axis_label_text('$p^\prime$ (%)')
plt.title('6) Polarization Fraction of 30 Doradus w/ Intensity in the background')
plt.show()

#7) HAWC+ POLARIZATION MAPS
print('7) HWAC+ Polarization Maps:\n')
files = [afile, cfile, dfile, efile]
titles = ['A','C','D','E']

for file, title in zip(files, titles):
    make_polmap(file, title)
