# Reading/printing ASCII data
from astropy.io import ascii

# Open file
f = open('sources.txt', 'r')

# Read and ignore header lines
header1 = f.readline()
header2 = f.readline()
header3 = f.readline()

# Loop over lines and extract variables of interest
for line in f:
    line = line.strip()
    columns = line.split()
    name = columns[2]
    j = float(columns[3])
    print(name, j)

f.close()