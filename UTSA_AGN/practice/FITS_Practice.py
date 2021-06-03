# Create/open a FITS file using astropy.io.fits -> organized for CIGALE format
from astropy.io import fits
import numpy as np
import os
from pathlib import Path
import pandas as pd
import csv

# create a simple sequence of floats from 0.0 -> 99.9
n = np.arange(100.0)

# create primaryHDU object to encapsulate all the data
hdu = fits.PrimaryHDU(n)

# read in input file and seperate into different arrays to be read into FITS file
csv_file = 'mcg52316c.csv'

test = []
id = []
redshift = []
filter = []
filter_err = []
filter1 = []
filter_err1 = []
filter2 = []
filter_err2 = []
filter3 = []
filter_err3 = []
filter4 = []
filter_err4 = []
filter5 = []
filter_err5 = []
filter6 = []
filter_err6 = []
filter7 = []
filter_err7 = []
filter8 = []
filter_err8 = []
filter9 = []
filter_err9 = []
filter10 = []
filter_err10 = []
filter11 = []
filter_err11 = []
filter12 = []
filter_err12 = []
distance = []

with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        test.append(row.get('test'))
        id.append(row.get('id'))
        redshift.append(row.get('redshift'))
        filter.append(row.get('filter'))
        filter_err.append(row.get('filter_err'))
        filter1.append(row.get('filter1'))
        filter_err1.append(row.get('filter_err1'))
        filter2.append(row.get('filter2'))
        filter_err2.append(row.get('filter_err2'))
        filter3.append(row.get('filter3'))
        filter_err3.append(row.get('filter_err3'))
        filter4.append(row.get('filter4'))
        filter_err4.append(row.get('filter_err4'))
        filter5.append(row.get('filter5'))
        filter_err5.append(row.get('filter_err5'))
        filter6.append(row.get('filter6'))
        filter_err6.append(row.get('filter_err6'))
        filter7.append(row.get('filter7'))
        filter_err7.append(row.get('filter_err7'))
        filter8.append(row.get('filter8'))
        filter_err8.append(row.get('filter_err8'))
        filter9.append(row.get('filter9'))
        filter_err9.append(row.get('filter_err9'))
        filter10.append(row.get('filter10'))
        filter_err10.append(row.get('filter_err10'))
        filter11.append(row.get('filter11'))
        filter_err11.append(row.get('filter_err11'))
        filter12.append(row.get('filter12'))
        filter_err12.append(row.get('filter_err12'))
        distance.append(row.get('distance'))

print(test)
print(id)
print(type(id[0]))
print(type(redshift[0]))
print(type(filter6[0]))
#print(type(id[0]))
print(redshift)
print(distance)

#a = np.array(id, dtype = np.float64)
#print(a)
#print(type(a))




# f = open('mcg_5_23_16.txt', 'r')
# for line in f:
#     line = line.strip()
#     columns = line.split()
#     id = float(columns[0])
#     print('id:', id)


# dir = "/Users/masonleist/PycharmProjects/AGN2/UTSA_AGN/practice/mcg_5_23_16.txt"
# data = pd.read_csv(dir, sep = '\t', header = None)
# df = pd.DataFrame(data)
#
#
# print(df[0])
# print(df[1])
#
# # individual input arrays -> df.iloc[row, column]
# id = df.iloc[:,0].tolist()
# redshift = df.iloc[:,1].tolist()
#
# # here, add/subtract filters/filterErr depending on input file/column number
# filter = df.iloc[:,2].tolist()
# filterErr = df.iloc[:,3].tolist()
# filter1 = df.iloc[:,4].tolist()
# filterErr1 = df.iloc[:,5].tolist()
# filter2 = df.iloc[:,6].tolist()
# filterErr2 = df.iloc[:,7].tolist()
# filter3 = df.iloc[:,8].tolist()
# filterErr3 = df.iloc[:,9].tolist()
# filter4 = df.iloc[:,10].tolist()
# filterErr4 = df.iloc[:,11].tolist()
# filter5 = df.iloc[:,12].tolist()
# filterErr5 = df.iloc[:,13].tolist()
# filter6 = df.iloc[:,14].tolist()
# filterErr6 = df.iloc[:,15].tolist()
# filter7 = df.iloc[:,16].tolist()
# filterErr7 = df.iloc[:,17].tolist()
# filter8 = df.iloc[:,18].tolist()
# filterErr8 = df.iloc[:,19].tolist()
# filter9 = df.iloc[:,20].tolist()
# filterErr9 = df.iloc[:,21].tolist()
# filter10 = df.iloc[:,22].tolist()
# filterErr10= df.iloc[:,23].tolist()
# filter11 = df.iloc[:,24].tolist()
# filterErr11 = df.iloc[:,25].tolist()
# filter12 = df.iloc[:,26].tolist()
# filterErr12 = df.iloc[:,27].tolist()
#
# distance = df.iloc[:,28].tolist()


# # Test to make sure file is read correctly
# print(id)
# print(redshift)
# print(filter)
# print(filterErr)
# print(distance)

# write our bin table more consisely w/o creating intermediate variables for the individual columns and
# w/0 manually creating a ColDefs object:
# here, add/subtract filters/filterErr depending on input file
hdu = fits.BinTableHDU.from_columns([fits.Column(name = 'id', format = '15A', array = id),
                                     fits.Column(name = 'redshift', format = 'D', array = redshift),
                                     fits.Column(name = 'WISE1', format = 'D', array = filter),
                                     fits.Column(name = 'WISE1_err', format = 'D', array = filter_err),
                                     fits.Column(name = 'WISE2', format = 'D', array = filter1),
                                     fits.Column(name = 'WISE2_err', format = 'D', array = filter_err1),
                                     fits.Column(name = 'WISE3', format = 'D', array = filter2),
                                     fits.Column(name = 'WISE3_err', format = 'D', array = filter_err2),
                                     fits.Column(name = 'WISE4', format = 'D', array = filter3),
                                     fits.Column(name = 'WISE4_err', format = 'D', array = filter_err3),
                                     fits.Column(name='sofia.forcast.F315', format='D', array=filter4),
                                     fits.Column(name='sofia.forcast.F315_err', format='D', array=filter_err4),
                                     fits.Column(name='sofia.hawc.A53', format='D', array=filter5),
                                     fits.Column(name='sofia.hawc.A53_err', format='D', array=filter_err5),
                                     fits.Column(name = 'herschel.pacs.70', format = 'D', array = filter6),
                                     fits.Column(name = 'herschel.pacs.70_err', format = 'D', array = filter_err6),
                                     fits.Column(name='sofia.hawc.C89', format='D', array=filter7),
                                     fits.Column(name='sofia.hawc.C89_err', format='D', array=filter_err7),
                                     fits.Column(name='sofia.hawc.D154', format='D', array=filter8),
                                     fits.Column(name='sofia.hawc.D154_err', format='D', array=filter_err8),
                                     fits.Column(name = 'herschel.pacs.160', format = 'D', array = filter9),
                                     fits.Column(name = 'herschel.pacs.160_err', format = 'D', array = filter_err9),
                                     fits.Column(name = 'herschel.spire.PSW', format = 'D', array = filter10),
                                     fits.Column(name = 'herschel.spire.PSW_err', format = 'D', array = filter_err10),
                                     fits.Column(name = 'herschel.spire.PMW', format = 'D', array = filter11),
                                     fits.Column(name = 'herschel.spire.PMW_err', format = 'D', array = filter_err11),
                                     fits.Column(name='herschel.spire.PLW', format='D', array=filter12),
                                     fits.Column(name='herschel.spire.PLW_err', format='D', array=filter_err12),
                                     fits.Column(name = 'distance', format = 'D', array = distance)])
#
# #writing this table to the HDU directly
#hdu.writeto('mcg14.fits')
#
#
#
#
#
