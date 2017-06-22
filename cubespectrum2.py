#! /usr/bin/env python
#
#  Load a FITS cube , extract the spectrum at a (or reference) pixel
#  and operate and plot some and then more...
#
#
#  22-jun-2017  PJT             summer project - cloned off cubespectrum.py
#
#
#  @todo
#     - have optional RSTFREQ as 3rd argument
#     - resample the gauss finer (not 5 points but may be 10x more?)


import os, sys
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.units import Quantity
c = 299792.458     # there should be a way to get 'c' from astropy.units ?

if len(sys.argv) == 1:
    print("Usage: %s fitsfile [xpos ypos]" % sys.argv[0])
    sys.exit(1)

if len(sys.argv) > 3:
    fitsfile = sys.argv[1]
    pos = [int(sys.argv[2]),int(sys.argv[3])]
elif len(sys.argv) == 2:
    fitsfile = sys.argv[1]
    pos = None

# open the fits file
hdu = fits.open(fitsfile)
print(len(hdu))

# get a reference to the header and data.  Data should be 3dim numpy array now
h = hdu[0].header
d = hdu[0].data.squeeze()
print(d.shape)

#  grab the restfreq, there are at least two way how this is done
if 'RESTFRQ' in h:
    restfreq=h['RESTFRQ']
elif 'RESTFREQ' in h:
    restfreq=h['RESTFREQ']
else:
    restfreq=None
print("RESTFREQ",restfreq)

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
channelf = (channeln-crpix3+1)*cdelt3 + crval3
channelv = (1.0-channelf/restfreq) * c

# what we plot
channel = channelv

ipeak = flux.argmax()
xpeak = channel[ipeak]
ypeak = flux[ipeak]

# moments around the peak
m = 5
x = channel[ipeak-m:ipeak+m]
y = flux[ipeak-m:ipeak+m]
xmean = (x*y).sum() / y.sum()
xdisp = (x*x*y).sum() / y.sum() - xmean*xmean
print("MOMENTS:",xmean,xdisp)
ymodel = ypeak * np.exp(-0.5*(x-xmean)**2/xdisp)



plt.figure()
plt.plot(channel,flux,'o-',markersize=2,label='data')
plt.plot(channel,zero)
plt.plot(x,ymodel,label='gauss')
plt.xlabel("Velocity (km/s)")
plt.ylabel("Flux")
plt.title("Spectrum at position %g %g" % (xpos,ypos))
plt.legend()
plt.show()


