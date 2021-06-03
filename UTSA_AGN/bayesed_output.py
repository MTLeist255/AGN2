# Simple program to plot the output SEDs of BAYESED
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams['figure.figsize'] = (8, 6)
waveA, fluxA, errorA = np.loadtxt('bayesed_output/ngc4151B_photom.txt', unpack = True)
wave1, waveB, fluxB = np.loadtxt('bayesed_output/model_sum.txt', unpack = True)
wave2, waveC, fluxC, trans1, trans2, trans3 = np.loadtxt('bayesed_output/ssh.txt', unpack = True)
wave3, waveD, fluxD, trans4, trans5, trans6 = np.loadtxt('bayesed_output/clumpy_torus.txt', unpack = True)

# uncomment for gb/bb additions
# bb
wave3, waveE, fluxE, trans7, trans8, trans9 = np.loadtxt('bayesed_output/bb.txt', unpack = True)
#
# gb
wave4, waveF, fluxF, trans10, trans11, trans12 = np.loadtxt('bayesed_output/gb.txt', unpack = True)

# seperate by input
waveG, fluxG, errorG = np.loadtxt('bayesed_output/ngc4151BHST_photom.txt', unpack = True)
waveH, fluxH, errorH = np.loadtxt('bayesed_output/ngc4151Bukirt_photom.txt', unpack = True)
waveI, fluxI, errorI = np.loadtxt('bayesed_output/ngc4151BSOFIA_photom.txt', unpack = True)
waveJ, fluxJ, errorJ = np.loadtxt('bayesed_output/ngc4151Bherschel_photom.txt', unpack = True)
# # converting mJy -> Jy: model sum
# for i in fluxB:
#     fluxB1 = fluxB/100000
#
# # converting mJy -> Jy: ssh
# for i in fluxC:
#     fluxC1 = fluxC/100000
#
# # converting mJy -> Jy: clumpy
# for i in fluxD:
#     fluxD1 = fluxD/100000

# # Converting Jy -> mJy
# for i in fluxA:
#    fluxA1 = fluxA*100000

#
# print('converted ->', waveA)
# print('flux ->', fluxC)
#
# # converting mJy -> Jy: add/subtract inputs as needed
# with open ('cenA/total_corrected.txt', 'w') as fh:
#     for a, b in zip(waveB, fluxB1):
#         print(' %.5f %.15f  ' % (a, b), file = fh)
# #
# wMS, fMS = np.loadtxt('cenA/total_corrected.txt', unpack = True)

# # Converting Jy -> mJy
# with open ('cenA/total_correctedPHOTOM.txt', 'w') as fh:
#     for a, b in zip(waveA, fluxA1):
#         print(' %.5f %.15f ' % (a, b), file = fh)
#
# wPhot, fPhot = np.loadtxt('cenA/total_correctedPHOTOM.txt', unpack = True)

# Calculating characteristic temperature of blackbody curve
# b & lambda_max are in meters(m)
lambda_max = 0.000111081
b = 0.0029
T = b / lambda_max

print('Characteristic BB Temp =', T)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

# PHOTOMETRY
# uncorrected photometry
#ax.scatter(waveA, fluxA, label = 'full photom')
# corrected photometry
#ax.scatter(wPhot, fPhot)

# individual photom inputs, uncomment as needed
ax.scatter(waveG, fluxG, marker = 'D', color = 'purple', label = 'HST')
ax.scatter(waveH, fluxH, marker = 'p', color = 'red', label = 'ukirt')
ax.scatter(waveI, fluxI, color = 'blue', label = 'SOFIA')
ax.scatter(waveJ, fluxJ, marker = 'X', color = 'green', label = 'herschel')

# MODEL SUM
# uncorrected model sum
ax.loglog(waveB, fluxB, color = 'red', label = 'sum')
# corrected model sum
#ax.loglog(wMS, fMS, color = 'red', label = 'sum')

# SSP -> un-comment whenever SSP is added
# uncorrected SSP
#ax.loglog(waveC, fluxC, color = 'blue', label = 'ssp')
# corrected SSP
# ax.loglog(wSSH, fSSH, color = 'blue', label = 'ssp')

# CLUMPY
# uncorrected clumpy
ax.loglog(waveD, fluxD, color = 'black', label = 'clumpy')
# corrected clumpy
# ax.loglog(wClum, fClum, color = 'black', label = 'clumpy')

# BB curve -> UNCOMMENT AS NEEDED
#ax.loglog(waveE, fluxE, color = 'green', label = 'bb')

# GB curve -> UNCOMMENT AS NEEDED
ax.loglog(waveF, fluxF, color = 'orange', label = 'graybody')

plt.xlabel('Wavelength[$\mu$m]')
plt.ylabel('Flux[Jy]')

#plt.title('Mrk 573')
# w/ BB temp
plt.title('NGC 4151')
plt.legend()

# stellar + clumpy sed
#plt.savefig('bayesed_output/sed_noNIR/ngc3081/ngc3081_total_SSPC.png')
# stellar + clumpy + bb sed
#plt.savefig('bayesed_output/sed_noNIR/ngc3081/ngc3081_total_SSPCBB.png')
# clumpy + bb sed
#plt.savefig('bayesed_output/sed_noNIR/ngc3081/ngc3081_total_CBB.png')
# stellar + clumpy + gb sed
#plt.savefig('bayesed_output/sed_noNIR/ngc3081/ngc3081_total_SSPCGB.png')
# clumpy + gb sed
#plt.savefig('bayesed_output/sed/ngc4151/ngc4151B_TRIM_CGB.png')

# Repeat process for knn outputs

# Repeat process for ann outputs

plt.show()
plt.close()