# REFERENCE: https://mwcraig.github.io/ccd-as-book/02-04-Combine-bias-images-to-make-master.html

from pathlib import Path
from astropy.nddata import CCDData
from astropy.io import fits

import ccdproc as ccdp
from astropy.stats import mad_std
import numpy as np
from astropy import units as u

IN_PATH = '/work/chuck/chen/obs/20190822'
IN_PATH = '/work/chuck/chen/obs/20190827'
IN_PATH = '/work/chuck/chen/obs/20190823'

outfolder = 'reduced'
reduced_path = Path(IN_PATH, outfolder)
reduced_path.mkdir(exist_ok=True)

flat_image_type = 'Flat Field'
naxis1 = 2749 # the size we needed
naxis2 = 2199 


print('read raw data')
raw_data_directory = Path(IN_PATH)
raw_data = ccdp.ImageFileCollection(raw_data_directory)


print('check bias frames')
for a_flat, fname in raw_data.hdus(imagetyp='Bias Frame',naxis1=naxis1,naxis2=naxis2, return_fname=True):
    print(f'  {fname}',a_flat.shape,'exposure', a_flat.header['EXPOSURE'], 'standard deviation ', a_flat.data.std())
calibrated_biases = raw_data.files_filtered(imagetyp='Bias Frame',naxis1=naxis1,naxis2=naxis2,include_path=True)

#print('copy out chosen frames to a new folder')
#import shutil
#for bias in calibrated_biases:
#    shutil.copy(bias,outfolder)


print('create master bias, if memomry overflow, please shut all the other applications!')
combined_bias = ccdp.combine(calibrated_biases,
                             method='average',
                             unit=u.adu, # CCDData requires a unit for the image if it is not in the header
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                             mem_limit=350e16
                            )

combined_bias.meta['combined'] = True
combined_bias.write(reduced_path / 'combined_bias.fit')


