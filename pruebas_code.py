#!/usr/bin/env python
# -*- coding: utf8 -*-

from astropy.io import fits
import time
from scipy.ndimage.interpolation import rotate
import matplotlib.pyplot as plt
import astropy.units as u
import os, glob
import numpy as np




print "It's working"
folderMean = "/Users/joseivan/Desktop/datos_pruebas/"
folderRot = folderMean+"Processed_Data/"
folderSub = folderMean+"submaps_Data/"
folderCubes = folderMean+'20110410/CubesbyHour/'



t0=time.time()


#~ hdu = fits.open("/Users/joseivan/Desktop/datos_pruebas/20110410/hmi.Ic_45s.20110410_230000_TAI.2.continuum.fits",checksum=True)
#~ data = hdu[1].data
#~ hdr = hdu[1].header
#~ rot_angle = hdr["CROTA2"]
#~ print rot_angle

#Create a list with 'Raw Data' 

RawData = sorted(glob.glob(folderMean+'20110410/*.fits'))

dateObs = (RawData[0])[-36:-21]
t1 = time.time()

# Rotating HMI continumm images
RawByH = RawData[0:80]
dcn1 = "ImaRot_"
print "Reading and rotating images..."
for i in RawByH:
	print "Reading and rotating image -->",i[-36:-21]
	data = (fits.open(i,checksum=True))[1].data
	hdr = (fits.open(i,checksum=True))[1].header
	rot_angle = hdr["crota2"]
	im = rotate(data,rot_angle,mode = "nearest",prefilter=False)
	np.save(folderRot+dcn1+i[-36:-21],im)
print 50*"*"
print "End rotation section.. \n"
print "Time elapsed rotating images -->", time.time() - t1
print 50*"*" 
del(data)
del(hdr)
del(im)

# Charging rotated images to create submaps

RotData = sorted(glob.glob(folderRot+"*.npy"))

t2 = time.time()

print "Reading rotated images to create submaps.."

dcn2 = "Submaps_"

dx0 = 0.104/0.504246

by_h=RotData[0:80]
x0=0
#Creating submaps
for j in by_h:
	print "Cutting image --> ",j[-19:-4] 
	
	im = np.load(j)
	dx = x0*dx0
	xi,xf = 500,900
	yi,yf = 2500,2800
	smap = im[yi:yf,xi+dx:xf+dx]
	np.save(folderSub+dcn2+j[-19:-4],smap)
	
	x0 = x0 + 1

print 50*"*"
print "End submapping section.. \n"
print "Time elapsed cutting images -->", time.time() - t2
print 50*"*" 

del(smap)
del(im)




print 50*"-"
print "End of code..."
print "Total time elapsed -->", time.time() - t0
print "_"



	
	
	










	
	








#~ data2 = smap.Map("/Users/joseivan/Desktop/datos_pruebas/20110410/hmi.Ic_45s.20110410_230000_TAI.2.continuum.fits")

#~ t1 = time.time()
#~ im1 = data2.rotate(angle=rot_angle*u.degree,use_scipy=True)
#~ tf1 = time.time() - t1
#~ print "Tiempo en imrotate -->", tf1

#~ t2 = time.time()
#~ im2 = rotate(data,rot_angle,mode = "nearest",prefilter=False)
#~ tf2 = time.time()-t2
#~ print "Tiempo en rotate -->", tf2


#~ plt.figure()
#~ plt.imshow(data,origin="lower",cmap="gray")
#~ plt.figure()
#~ plt.imshow(im1,origin="lower",cmap="gray")
#~ plt.figure()
#~ plt.imshow(im2,origin="lower",cmap="gray")
#~ plt.show()
