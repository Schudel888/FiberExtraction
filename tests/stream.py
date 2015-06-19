
from astropy.io import fits
import numpy as np

x = np.random.rand(43,55).astype(np.float64)

fn = 'streamed.fits'
header = fits.Header()
header['BITPIX']=-64
header['NAXIS']=3
header['NAXIS1']=5
header['NAXIS2']=43
header['NAXIS3']=55

i=0
with fits.StreamingHDU(fn, header) as hdu:
    while not hdu.write(x):
        i +=1
        print str(i)

hdulist=fits.open(fn)

print hdulist[0].data[0]