#!/usr/bin/env python
# -*- coding: utf8 -*-

from astropy.io import fits
import numpy as np
import os
import time
from bignfft import *
import shutil

def bigsonic(first,last,bxdim,bydim):
	if not os.path.exists("work"):
		os.makedirs("work")
	if not os.path.exists("filter"):
		os.makedirs("filter")
		
	str0 = "work/"
	strf = "filter/"
	str2 = "/home/hypnus1803/Desktop/Filtered_Images/"
	cubename = raw_input("Please, put the full path of the cube:")
	structure =fits.open(cubename)
	cube = structure[0].data
	dimx = cube.shape[2]
	dimy = cube.shape[1]
	xdim=dimx
	ydim=dimy
	x_anf=0
	y_anf=0
	scale = 0.504251
	
#-----------------------------------------------------------------

	if xdim%2 != 0: xdim = xdim - 1	
	if ydim%2 != 0: ydim = ydim - 1	

	if (last - first + 1)%2 != 0: last = last - 1	

	tdim = last - first + 1

	print "-------"
	#Mean time separation between images [s]
	t_step = 45.
	
	# Maximum phase velocity [km/s]
	v_ph = 4.
	
	ap = 0
	
	cut = 0
	
	t1 = time.time()
	
	#-----------------------------------------
	
	print tdim, " images to be filtered"
	
	# subsonic construction
	
	#Definin' unit in the Fourier domain
	
	print "Spatial resolution ->", scale, " arcsec/pixel"
	
	print "---"
	
	kx_step = 1./(scale*725.*xdim)
	ky_step = 1./(scale*725.*ydim)
	w_step = 1./(t_step*tdim)
	
	
	# Prepare de filter
	
	nx = xdim/2 + 1
	ny = ydim/2 + 1
	nt = tdim/2 + 1
	filter_mask = np.zeros([ny,nx,nt],dtype=float)
	
	print "Now calculatin' filter..."
	for j in range(ny):
		for i in range(nx):
			k_by_v = np.sqrt((i*kx_step)**2+(j*ky_step)**2)*v_ph
			for k in range(1,nt):
				if k*w_step <= k_by_v:
					filter_mask[j,i,k]=1.
	# Create  slices (kx,ky) of the filter for different w-values
	
	for i in range(nt):
		filter_slice = np.zeros([ydim,xdim])
		if i == 0:
			filter_slice = filter_slice + 1.
		else:
			filter_slice[0:ny,0:nx]=filter_mask[:,:,i]
			filter_slice[ny:ydim,:]=filter_slice[1:ny-1,:][::-1,:]
			filter_slice[:,nx:xdim]=filter_slice[:,1:nx-1][:,::-1]
		dcn = str(i).zfill(3)
		np.save(strf+dcn,filter_slice)
		
		if i != 0:
			dcn = str(tdim-i+first).zfill(3)
			np.save(strf+dcn,filter_slice)
	print "The filter was written to "+strf
	
	del(filter_mask)
	del(filter_slice)
	
	# Loop of reading, optional apodization and writing images
	for n in range(first,last+1):
		dcn = str(n).zfill(3)
		print "Reding images -->", "cube["+str(n)+"]"
		ima = cube[n,:,:]
		ima = ima.astype(float)
		np.save(str0+"apo"+dcn,ima)
		print "apodized images were written to "+str0
	
	del(ima)
	
	# Direct FFT
	
	print 50*"="
	print "Calling bignfft"
	print 50*"-"
	
	bignfft(-1, first,last,dimx,dimy,bxdim,bydim)
	
	# Multiplying transformed images by filter images
	
	ima = np.zeros([ydim,xdim],dtype="complex64")
	filter = np.zeros([ydim,xdim])
	
	for n in range(first,last+1):
		dcn = str(n).zfill(3)
		ima = np.load(str0+"fft"+dcn+".npy")
		filter = np.load(strf+dcn+".npy")
		
		print "multiplying by the filter ",str0+"fft"+dcn
		
		ima = ima*filter
		np.save(str0+"fft"+dcn,ima)
	del(ima)
	del(filter)
	
	#Inverse FFT
	print 50*"="
	print "Calling bignfft"
	print 50*"-"
	
	bignfft(1, first,last,dimx,dimy,bxdim,bydim)
	
	# Saving results
	
	
	ima = np.zeros([ydim,xdim],dtype="complex64")
	
	cube_new = np.zeros([tdim,dimy,dimx],dtype=int)
	
	for n in range(first,last+1):
		dcn = str(n).zfill(3)
		ima = np.load(str0+"F"+dcn+".npy")
		im = (ima.real).astype(int)
		cube_new[n,:,:] = im
	
	hdu = fits.PrimaryHDU(cube_new)
	hdulist = fits.HDUList([hdu])
	hdulist.writeto(str2+'cube_filered.fits')
	
	print "---"
	print "Total elapsed time from begining = ", np.round(time.time()-t1,2)
	print " "
	print "Erasing directories work and filter"
	
	shutil.rmtree("work/")
	shutil.rmtree("filter")

bigsonic(0,80,260/2,450/2)
		
	
	
	
	
	
		
	
		
	
	
	
	
	
	
