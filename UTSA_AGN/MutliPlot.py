# Recreate multi plot from Fuller AAS 237
import matplotlib.pyplot as plt
import numpy as np


# load redshift correct data
wave1, flux1, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_mrk231_trim.txt', unpack = True)
wave2, flux2, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_circinus_hires_trim.txt', unpack = True)
wave3, flux3, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_ngc3281_trim.txt', unpack = True)
wave4, flux4, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_ngc5506_trim.txt', unpack = True)
wave5, flux5, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_cenA_trim.txt', unpack = True)
# wave6, flux6, errC = np.loadtxt('CASSIS/CASSIS_FINAL/cassis_eso323g77_trim.txt', unpack = True)

# flux scaling
for i in flux1:
    f1 = (1.0 + flux1) / flux1[242]

for i in flux2:
    f2 = (-2.0 + flux2) / flux2[875]

for i in flux3:
    f3 = (0.0 + flux3) / flux3[235]

for i in flux4:
    f4 = (-0.9 + flux4) / flux4[234]

for i in flux5:
    f5 = (-1.5 + flux5) / flux5[233]

plt.plot(wave1, f1, label = 'Mrk 231')
plt.plot(wave2, f2, label = 'Circinus')
plt.plot(wave3, f3, label = 'NGC 3281')
plt.plot(wave4, f4, label = 'NGC 5506')
plt.plot(wave5, f5, label = 'Centaurus A')
#plt.plot(wave6, flux6, label = '')
plt.legend()
plt.title('10 and 18 μm Absorption Feature')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Scaled Flux[Jy]')
plt.savefig('CASSIS/1018um.png')
plt.show()

