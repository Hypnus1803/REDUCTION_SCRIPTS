from astropy.io import fits
import time
from scipy.misc import imrotate
from scipy.ndimage.interpolation import rotate
import matplotlib.pyplot as plt
import sunpy.map as smap

hdu = fits.open("/home/hypnus1803/Desktop/hmi.Ic_45s.20110413_075830_TAI.2.continuum.fits",checksum=True)
data = hdu[1].data
hdr = hdu[1].header
rot_angle = hdr["CROTA2"]
print rot_angle
#~ data2 = smap.Map("/home/hypnus1803/Desktop/hmi.Ic_45s.20110413_075830_TAI.2.continuum.fits")
#~ 
#~ t1 = time.time()
#~ im1 = data2.rotate()#angle=rot_angle)#,use_scipy=True)
#~ tf1 = time.time() - t1
#~ print "Tiempo en imrotate -->", tf1

t2 = time.time()
im2 = rotate(data,rot_angle,mode = "nearest",prefilter=False)
tf2 = time.time()-t2
print "Tiempo en rotate -->", tf2


#~ plt.figure()
#~ plt.imshow(data,origin="lower",cmap="gray")
#~ plt.figure()
#~ plt.imshow(im1,origin="lower",cmap="gray")
#~ plt.figure()
#~ plt.imshow(im2,origin="lower",cmap="gray")
#~ plt.show()


