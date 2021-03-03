import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image

import matplotlib.pyplot as plt


fig = plt.figure(figsize=(12,18))
ima = Image.open('Current_AGN_fits/silicate/9_18_continuum/test.png')
imb = Image.open('Current_AGN_fits/silicate/9_18_continuum/test3.png')
imc = Image.open('Current_AGN_fits/silicate/9_18_continuum/test4.png')
# Display row 1
im1 = plt.figimage(ima, xo=0, yo=0)
# Display row 2
im2 = plt.figimage(imb, xo=0, yo=600)
# Display row 3
im3 = plt.figimage(imc, xo=0, yo=1200)

plt.show()