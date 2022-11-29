# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.io import fits
from matplotlib.colors import LogNorm


PATH = '/work/chuck/chen/obs'
try:
  img = f'{PATH}/20190822/reduced/light_ab_and-0018V5.fit'
  image_data, header = fits.getdata(img,header=True)
except:
  print('no real image, use fake one')
  from astropy.utils.data import get_pkg_data_filename
  img = get_pk('tutorials/FITS-images/HorseHead.fits')
  image_data, header = fits.getdata(img,header=True)
#hdu = fits.open(img)
#image_data = hdu[0].data
#image_data = fits.getdata(image_file)

#hdu.info()
#hdu[0].header 
#print('filter',hdu[0].header['FILTER'])
#print('midpoint time',header['JD-HELIO'])

print(image_data.shape)
plt.figure()
plt.imshow(image_data, cmap='gray',norm=LogNorm())
#plt.imshow(image_data, cmap='gray')
plt.colorbar()
plt.show()
