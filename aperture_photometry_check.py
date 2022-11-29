# Photometry
# REFERENCE : 
# - PHOTOMETRY https://photutils.readthedocs.io/en/stable/getting_started.html
# - WRITING TABLES https://docs.astropy.org/en/stable/io/ascii/write.html

import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
import photutils

from photutils import DAOStarFinder
from astropy.stats import mad_std
from photutils import aperture_photometry, CircularAperture
from photutils.utils import make_random_cmap
from astropy.io import ascii

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rcParams['figure.figsize'] = [15, 15]


PATH = '/work/chuck/chen/obs'
img = f'{PATH}/20190827/reduced/light_ab_and-0022V4.fit'
img = f'{PATH}/20190822/reduced/light_ab_and-0022V4.fit'

fwhm = 11.
sigma = 5.
# fitshape is an odd number for the rectangular shape used for the fitting. It should be as large as the biggest star in the image
fitshape = (51,51)
de_blending = True
de_blending = False


print('checking aperture photometry settings by single plot')
#PATH = '/work/chuck/jortuno/Observation_Laboratory_Course/20190823/Red_images'

hdul = fits.open(img)
#hdul.info()
#hdul[0].header 
#print('filter',hdul[0].header['FILTER'])
#print('midpoint time',hdul[0].header['JD-HELIO'])
image, header = fits.getdata(img,header=True)
print(repr(hdul[0].header))


bkg_sigma = mad_std(image)  
#daofind = DAOStarFinder(fwhm=10., threshold=500)#5.*bkg_sigma)  
daofind = DAOStarFinder(fwhm=fwhm, threshold=sigma*bkg_sigma)  
sources = daofind(image)  
sources['JD-HELIO'] = header['JD-HELIO']
for col in sources.colnames:  
    sources[col].info.format = '%.8g'  # for consistent table output
print(sources)  


print('  do aperture photometry, de-blending =', de_blending)
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))  
if(de_blending):
  # https://photutils.readthedocs.io/en/stable/api/photutils.psf.BasicPSFPhotometry.html#photutils.psf.BasicPSFPhotometry
  # groupmaker : photutils.psf.DAOGroup. It will set some star groups if they are too close
  daogroup = photutils.psf.DAOGroup(fwhm*2.0)
  # bkg_estimator : photutils.background.MMMBackground or photutils.background.MedianBackground
  median_bkg = photutils.background.MedianBackground()
  mmm_bkg = photutils.background.MMMBackground()
  # psf_model : photutils.psf.sandbox.DiscretePRF or photutils.psf.IntegratedGaussianPRF 
  #      photutils.psf.IntegratedGaussianPRF needs gaussian sigma, so use: astropy.stats.gaussian_fwhm_to_sigma
  #      photutils.psf.sandbox.DiscretePRF needs 
  gaussian_psf_model = photutils.psf.IntegratedGaussianPRF(sigma = fwhm*(astropy.stats.gaussian_fwhm_to_sigma) )
  #gaussian_psf_model = photutils.psf.IntegratedGaussianPRF(sigma = 15.0) 
  # finder : photutils.detection.DAOStarFinder or photutils.detection.IRAFStarFinder (this last one does NOT work fine)
  # fitter : Least square method for fitting the model to the data. For instance: 
  #   astropy.modeling.fitting.LevMarLSQFitter
  #   astropy.modeling.fitting.LinearLSQFitter 
  LevMar_fitter = astropy.modeling.fitting.LevMarLSQFitter()
  linear_fitter = astropy.modeling.fitting.LinearLSQFitter()
  # aperture_radius, i.e., radius (in units of pixels) used to compute initial estimates for the fluxes of sources
  aperture_radius = fwhm
  apertures = photutils.psf.BasicPSFPhotometry(group_maker=daogroup, 
                  bkg_estimator=mmm_bkg,
                  psf_model=gaussian_psf_model,
                  fitshape=fitshape,
                  finder=daofind,
                  fitter=LevMar_fitter,
                  aperture_radius=aperture_radius)
  phot_table = apertures(image)
  print('residue image')
  plt.imshow(apertures.get_residual_image(), cmap='viridis',
  aspect=1, interpolation='nearest', origin='lower')  

  #star_groups = daogroup(phot_table)
  #plt.imshow(apertures.get_residual_image(), origin='lower', interpolation='nearest',cmap='Greys_r')
  #cmap = make_random_cmap()
  #for i, group in enumerate(star_groups.groups):
  #  xypos = np.transpose([group['x_0'], group['y_0']])
  #  ap = CircularAperture(xypos, r=fwhm)
  #  ap.plot(color=cmap.colors[i])

else:
  apertures = CircularAperture(positions, r=fwhm)
  phot_table = aperture_photometry(image, apertures) 

  #plt.imshow(image, cmap='gray', norm=LogNorm())
  plt.imshow(image, cmap='gray_r', origin='lower')
  apertures.plot(color='blue',lw=20, alpha=0.5)

for col in phot_table.colnames:  
    phot_table[col].info.format = '%.8g'  # for consistent table output
ascii.write(phot_table, f'ap_check.dat', overwrite=True)
print(phot_table)  
plt.show()

