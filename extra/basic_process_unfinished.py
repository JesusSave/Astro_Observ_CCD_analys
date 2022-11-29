#===================================================================
#
#           OBSERVATIONAL FITS FILES REDUCTION
#
# This script is used to overscan, correct, trim, gain correct, 
# substract master bias, master dark, and flat-field the FIT img.
# The script rearanges automatically.
#
# WARNING: MAINTENANCE SUSPEND!
#          The fit file stacking and substracing may be wrong, 
#          Please use ccdproc combine and subtraction instead !!
#
# VERSION: 27 Sep 2020
# AUTHOER: QIANG CHEN chen@camk.edu.pl 
#
# PYTHON ENVIRONMENT REQUIREMENT:
# - pip install astropy
# - pip install --no-deps ccdproc
# alternatively can use (CAMK):
#   source /home/chen/python-chen/chen3.6/bin/activate
#
# REFERENCE:
# - (ATUO) Reduction toolbox (astropy) https://ccdproc.readthedocs.io/en/latest/reduction_toolbox.html
# - (RECOMMENDED) REDUCTION EXAMPLE https://nbviewer.jupyter.org/gist/mwcraig/06060d789cc298bbb08e
# - CCD DATA REDUCTION GUIDE https://mwcraig.github.io/ccd-as-book/00-00-Preface.html
# - SAVE TO FITS FORMAT (astropy) https://docs.lightkurve.org/tutorials/03-making-fits-files.html
# - EDIT HEADER http://learn.astropy.org/FITS-header.html
#===================================================================

import glob, os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.nddata import CCDData
import ccdproc
from ccdproc import combine
from matplotlib.colors import LogNorm

PATH = '/work/chuck/chen/obs'

folder = '20190822'
data_shape = (2199, 2749)
filters = ['V','B','U','I']

num_fig = 2 # should >=5, low memory @_@

biass = []
darks = []
flats = []
lights = []
excludes = []

imgs = glob.glob(f'{PATH}/{folder}/*.fit')
print('TOTAL IMAGES',len(imgs))

print('Removing strange and re-aranging images')
for img in imgs:
  data,header = fits.getdata(img,header=True)
  #print('  img',img, imgs.index(img),'/',len(imgs))
  if (data.shape != data_shape):
    #print('  Exclude:',img,', data shape',data.shape,', TYP',header['IMAGETYP'],', EXPTIME :', header['EXPTIME'])
    excludes.append(img)
  else:
    if (header['IMAGETYP'] == 'Dark Frame'):
      darks.append(img)
    elif (header['IMAGETYP'] == 'Light Frame'):
      lights.append(img)
    elif (header['IMAGETYP'] == 'Bias Frame'):
      biass.append(img)
    elif (header['IMAGETYP'] == 'Flat Field'):
      flats.append(img)
      #print('  img',img, data.shape, header['IMAGETYP'],'filter',header['FILTER'],', EXPTIME :', header['EXPTIME'])
    else:
      print('  strange:',img,', data shape',data.shape,'TYP',header['IMAGETYP'],', EXPTIME :', header['EXPTIME'])

print('  Excludes',len(excludes))
print('  Total residue',len(imgs))
print('  Dark img', len(darks))
print('  Light img',len(lights))
print('  Bias img',len(biass))
print('  Flats img',len(flats))
#imgs =  [x for x in imgs if x not in excludes]


print('Checking bias and creating master bias')
bias_list = []
for img in biass[0:num_fig]:
  data,header = fits.getdata(img,header=True)
  print('  img :',img,', data shape',data.shape,', EXPTIME :', header['EXPTIME'])
  #hdul.info()
  #print(repr(hdul[0].header))
  ccd = CCDData.read(img, unit=u.adu)
  bias_list.append(ccd)

master_bias = combine(bias_list, method='median')
#biases = ccdproc.Combiner(bias_list)
#master_bias = biases.average_combine()



print('Checking DARKs and creating master dark')
darks_list = []
expts = []
for img in darks[0:num_fig]:
  data,header = fits.getdata(img,header=True)
  dark_exptime = header['EXPTIME']
  expts.append(dark_exptime)
  print('  img :',img,data.shape,', exposure time :', dark_exptime)
  ccd = CCDData.read(img, unit=u.adu)
  darks_list.append(ccd)

if (len(set(expts)) != 1):
  print('ERROR: darks have different exposure times!')
  exit()

img_concat = [fits.getdata(img) for img in darks]
master_dark = combine(darks_list, method='median')
#master_dark = CCDData(master_dark, unit=u.electron)
master_dark.header['exposure'] = dark_exptime



excludes_flats = []
excludes_lights = []
for filt in filters:
  print('PROCESS BY FILTER TYPE:', filt)
  print('  making master flat')
  flats_sub = []
  flats = [x for x in flats if x not in excludes_flats]
  for img in flats[0:num_fig]:
    data,header = fits.getdata(img,header=True)
    if (header['FILTER'] == filt):
      print('  img:',img, data.shape,filt,', exposure time :',header['EXPTIME'])
      flats_sub.append(img, unit=u.adu)
      excludes_flats.append(img)
    img_concat = [fits.getdata(img)-master_bias-master_dark for img in flats_sub]
    master_flat = CCDData(flats_sub, unit=u.electron)
    #master_flat = CCDData(master_flat, unit=u.electron)


  print('  Basic Processing Loop')
  for img in lights:
    data,header = fits.getdata(img,header=True)
    if (header['FILTER'] == filt):
      print('    i',lights.index(img),'/',len(lights))
      print('      img :',img,fits.getdata(img).shape,', filter :',header['FILTER'],', exposour time :', header['EXPTIME'])
      #excludes_lights.append(img)
      path = os.path.dirname(os.path.realpath(img))
      infile = os.path.basename(img)
      outpath = f'{path}/output'
      if not os.path.exists(outpath):
        os.makedirs(outpath)
      outfile = f'{outpath}/reduced_'+infile

      ccd = CCDData(data, unit=u.adu)
      #print(ccd.shape)
      ccd.header['exposure'] = header['EXPTIME']  # for dark subtraction
      nccd = ccdproc.ccd_process(ccd,# oscan='[201:232,1:100]',
                           #trim='[1:200, 1:100]',
                           error=True,
                           gain=1.0*u.electron/u.adu, # 1
                           readnoise=0*u.electron, # 0
                           bias_frame=master_bias,
                           dark_frame=master_dark,
                           exposure_key='exposure',
                           exposure_unit=u.second,
                           dark_scale=True,
                           master_flat=master_flat)

      print('      Updating header and Saving figure to FIT file')
      header['SWOWNER'] = "Q. Chen"
      fits.writeto(outfile, data, header, clobber=True)

