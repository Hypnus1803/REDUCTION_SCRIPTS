 #!/usr/bin/env python
# -*- coding: utf8 -*-


from astropy.io import fits
import time
from scipy.ndimage.interpolation import rotate
import matplotlib.pyplot as plt
import astropy.units as u
import os, glob
import numpy as np


# Making cubes...

print "Creating a cube data with correct dimensions and a pseudo tracking"

folderMean = "/Users/joseivan/Desktop/datos_pruebas/"
folderSub = folderMean+"submaps_Data/"
folderCubes = folderMean+'20110410/CubesbyHour/'

#Reading submaps
Submaps = sorted(glob.glob(folderSub+"*.npy"))

t3 = time.time()

subBy_h = Submaps[0:80]

cube = np.zeros([80,300,400],dtype=float)

#~ print subBy_h[0]

for k in range(len(subBy_h)):
	smap = np.load(subBy_h[k])
	cube[k,:,:] = smap
	

hdu = fits.PrimaryHDU(cube)
hdulist = fits.HDUList([hdu])
hdulist.writeto(folderCubes+'cube_filered'+subBy_h[0][-19:-4]+'.fits')
print "Time writing cube -->", time.time()-t3
