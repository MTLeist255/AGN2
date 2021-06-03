# read in micron then convert to angstrom
import matplotlib.pyplot as plt
import numpy as np
import csv

waveC, fluxC = np.loadtxt('sofia/herschel/spire/PSW.pb', unpack = True)
#print('original ->', waveC)
# print('flux ->', fluxC)

# convert micron to angstrom

for i in waveC:
    waveA = waveC * 0.0001

#
#print('converted ->', waveA)
# print('flux ->', fluxC)
#
# # take wavelength (angstrom) and flux error and create new txt file
with open ('sofia/herschel/spire_correct/PSW.txt', 'w') as fh:
    for a, b in zip(waveA, fluxC):
        print('%.4f %.17f' % (a, b), file = fh)

# test
# waveC2, fluxC2 = np.loadtxt('../FOR_F054N.txt', unpack = True)
# print('new output txt ->', waveC2)