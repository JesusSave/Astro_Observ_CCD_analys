#===================================================================
#
#         APERTURE PHOTOMETRY
#
# This script is used to apply basic aperture photometry to the reduced .fits files.
# The stars are ruled out by FWHM and sigma value.
#
# "roundness" is a measure of the difference of the second moments in each direction (zero means "round", +1 or -1 means "elongated),
# "sharpness" is a measure of the central pixel value to its neighbors.
#
# The "roundness" is a measure of the symmetry of the stellar image. It is defined as
#                      xwidth - ywidth
#   roundness =  2 * (-----------------)
#                      xwidth + ywidth
# The "sharpness" parameter is similar to, but not exactly the same as, the sharpness parameter in DAOPHOT. It is here
#
#  (image value at star centroid) - (mean of image values around centroid)
#  -----------------------------------------------------------------------
#                    (image value at star centroid)
#
# WARNING:
# If the targets stars overlap significantly, should use de-blending techniques
# like PSF fitting photometry.
# Before running this, should use aperture_photometry_check.py to identify fwhm value, inorder get proper aperture 
#
# VERSION: 24 Sep 2020
# AUTHOER: QIANG CHEN chen@camk.edu.pl 
#
# PYTHON ENVIRONMENT REQUIREMENT:
# - pip install astropy
# - pip install --no-deps ccdproc
# - pip install photutils
# alternatively can use (CAMK):
#   source /home/chen/python-chen/chen3.6/bin/activate
#
# REFERENCE:
# - APERTURE PHOTOMETRY https://photutils.readthedocs.io/en/stable/getting_started.html
# - WRITING TABLES https://docs.astropy.org/en/stable/io/ascii/write.html
# - EDITING TABLES https://docs.astropy.org/en/stable/table/modify_table.html
#===================================================================

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import astropy
import photutils

from astropy.io import fits
from astropy.stats import mad_std
import photutils
from photutils import DAOStarFinder
from astropy.io import ascii
from photutils import aperture_photometry, CircularAperture


PATH = '/work/chuck/chen/obs'

folder,fwhm,sigma, de_blending = '20190822/reduced', 10., 7, False 
folder,fwhm,sigma, de_blending = '20190823/reduced', 11., 5, False 
folder,fwhm,sigma, de_blending = '20190827/reduced', 11., 5, False 

# fitshape is an odd number for the rectangular shape used for the fitting (de_blening mode). 
#    It should be as large as the biggest star in the image
fitshape = (51,51) 
# Master dark too big, science frame background negative! cannot de-blending!

print('read reduced fits files')
imgs = glob.glob(f'{PATH}/{folder}/light*.fit')

for image in imgs:
  print('  find stars in figure:',imgs.index(image)+1,'/',len(imgs))
  path = os.path.dirname(os.path.realpath(image))
  infile = os.path.basename(image)
  infilename = os.path.splitext(infile)[0]

  data, header = fits.getdata(image,header=True)
  #print('  midpoint time',header['JD-HELIO'])
  bkg_sigma = mad_std(data)  
  daofind = DAOStarFinder(fwhm=fwhm, threshold=sigma*bkg_sigma)  
  sources = daofind(data)
  sources['JD-HELIO'] = header['JD-HELIO']
  for col in sources.colnames:  
    sources[col].info.format = '%f'  # for consistent table output
  #print(sources)
  ascii.write(sources, f'{path}/stars_{infilename}.dat', overwrite=True)


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
    #        photutils.psf.IntegratedGaussianPRF needs gaussian sigma, so use: astropy.stats.gaussian_fwhm_to_sigma
    #        photutils.psf.sandbox.DiscretePRF needs 
    gaussian_psf_model = photutils.psf.IntegratedGaussianPRF(sigma = fwhm*(astropy.stats.gaussian_fwhm_to_sigma) )
    #gaussian_psf_model = photutils.psf.IntegratedGaussianPRF(sigma = 15.0) 
    # finder : photutils.detection.DAOStarFinder or photutils.detection.IRAFStarFinder (this last one does NOT work fine)
    # fitter : Least square method for fitting the model to the data. For instance: 
    #     astropy.modeling.fitting.LevMarLSQFitter
    #     astropy.modeling.fitting.LinearLSQFitter 
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
    phot_table = apertures(data)
  else:
    apertures = CircularAperture(positions, r=fwhm)
    phot_table = aperture_photometry(data, apertures) 

  phot_table['JD-HELIO'] = header['JD-HELIO'] 
  for col in phot_table.colnames:  
    phot_table[col].info.format = '%f'  # for consistent table output
  #print(phot_table) 
  ascii.write(phot_table, f'{path}/ap_{infilename}.dat', overwrite=True)
