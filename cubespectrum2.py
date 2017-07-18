#! /usr/bin/env python
#
#  Load a FITS cube , extract the spectrum at a (or reference) pixel
#  and operate and plot some and then more....
#
#
#  22-jun-2017  PJT             summer project - cloned off cubespectrum.py
#  july-2017    Thomas/Peter    various improvements
#
#  @todo
#     - have optional RESTFRQ or RESTFREQ as 3rd argument [done]
#     - output the spectrum in a table, much like testCubeSpectrum.tab [done]
#     - resample the gauss finer (not 5 points but may be 10x more?)


import os, sys, math
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.units import Quantity
c = 299792.458     # [km/s] there should be a way to get 'c' from astropy.units ?



#if len(sys.argv) == 1:
#    print("Usage: %s fitsfile [xpos ypos] [restfreq [vmin vmax]]" % sys.argv[0])
#    sys.exit(1)
#elif len(sys.argv) == 3:
#    print("need ypos")
#    sys.exit(1)
#elif len(sys.argv) > 3:
#    fitsfile = sys.argv[1]
#    pos = [int(sys.argv[2]),int(sys.argv[3])]
#elif len(sys.argv) == 2:
#    fitsfile = sys.argv[1]
#    pos = None




na = len(sys.argv)

if na == 7:
    fitsfile = sys.argv[1] 
    pos = [int(sys.argv[2]),int(sys.argv[3])]
    restfreq = int(sys.argv[4])
    print(restfreq)
    vmin = int(sys.argv[5])
    vmax = int(sys.argv[6])
    use_vel = True

elif na == 5: 
    fitsfile = sys.argv[1]
    pos = [int(sys.argv[2]),int(sys.argv[3])]
    restfreq = sys.argv[4]
    print(restfreq)
    use_vel = True

elif na == 4:
    fitsfile = sys.argv[1]
    pos = [int(sys.argv[2]),int(sys.argv[3])]
    restfreq = 112e9
    vmin = None
    use_vel = True

elif na == 4:
    fitsfile = sys.argv[1]
    pos = [int(sys.argv[2]),int(sys.argv[3])]
    use_vel = False

elif na == 2:
    fitsfile = sys.argv[1]
    pos = None
    use_vel = False

else:
    sys.exit(1)

use_vel = False



# open the fits file
hdu = fits.open(fitsfile)
print(len(hdu))

# get a reference to the header and data.  Data should be 3dim numpy array now
h = hdu[0].header
d = hdu[0].data.squeeze()
print(d.shape)


#  grab the restfreq, there are at least two ways how this is done
#if 'RESTFRQ' in h:
#    restfreq=h['RESTFRQ']
#elif 'RESTFREQ' in h:
#    restfreq=h['RESTFREQ']
#else:
#    restfreq=None
#print("RESTFREQ",restfreq)





if pos == None:
    # the FITS reference pixel is always a good backup
    xpos = int(h['CRPIX1'])
    ypos = int(h['CRPIX2'])
    print("No position given, using reference pixel %g %g" % (xpos,ypos))
else:
    xpos = pos[0]
    ypos = pos[1]





flux     = d[:,ypos,xpos]

nchan    = d.shape[0]       
channeln = np.arange(nchan)
zero     = np.zeros(nchan)

cdelt3 = h['CDELT3']
crval3 = h['CRVAL3']
crpix3 = h['CRPIX3']
# to convert the channel to frequency
channelf = (channeln-crpix3+1)*cdelt3 + crval3
# to convert the Frequency to velocity
#channelv = (1.0-channelf/restfreq) * c
#print (channelf)
#print (channelv)

# what we plot
#channel = channelv 
#channel = channelf
#channel = channeln

if use_vel:
# to convert the Frequency to velocity
    channelv = (1.0-channelf/restfreq) * c
    channel = channelv
else:
    channel = channelf


ipeak = flux.argmax()
xpeak = channel[ipeak]
ypeak = flux[ipeak]

# moments around the peak
m = 50
x = channel[ipeak-m:ipeak+m]
y = flux[ipeak-m:ipeak+m]
xmean = (x*y).sum() / y.sum()
xdisp = (x*x*y).sum() / y.sum() - xmean*xmean
if xdisp > 0:
    xdisp = math.sqrt(xdisp)
fwhm = 2.355 * xdisp    
print("MEAN/DISP/FWHM:",xmean,xdisp,fwhm)
ymodel = ypeak * np.exp(-0.5*(x-xmean)**2/(xdisp*xdisp))







if use_vel:
   plt.figure()  
   if vmin != None:
       channelv = ma.masked_outside(channelv,vmin,vmax)
       plt.xlim([vmin,vmax])
   plt.plot(channelv,flux,'o-',markersize=2,label='data')
   plt.plot(channelv,zero)
#  plt.plot(x,ymodel,label='gauss')
   plt.xlabel("Velocity (km/s)")
   plt.ylabel("Flux")
   plt.title(fitsfile + restfreq + "Spectrum at position %g %g" % (xpos,ypos))
   plt.legend()
   plt.show()

else:  
   plt.figure()
   plt.plot(channelf,flux,'o-',markersize=2,label='data')
   plt.plot(channelf,zero)
   plt.xlabel("Frequency (GHz)")
   plt.ylabel("Flux")
   plt.title(fitsfile + "Spectrum at position %g %g" % (xpos,ypos))
   plt.legend()
   plt.show()


#to create a table of the frequency and flux
xtab = channelf /1e9 #to set the freqency to GHz
ytab = flux 
np.savetxt('Frequency_Flux.tab',np.c_[xtab,ytab], delimiter='  ',header=("Frequency""       " "Flux"),comments='#',fmt='%.8f')













