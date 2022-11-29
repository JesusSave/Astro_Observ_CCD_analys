# REFERENCE : https://mwcraig.github.io/ccd-as-book/03-06-Combine-darks-for-use-in-later-calibration-steps.html


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


print('check dark frames')
for a_flat, fname in raw_data.hdus(imagetyp='Dark Frame', return_fname=True):
    print(f'  {fname}',a_flat.shape,'exposure', a_flat.header['EXPOSURE'], 'standard deviation ', a_flat.data.std())
calibrated_darks = raw_data.files_filtered(imagetyp='Dark Frame',naxis1=naxis1,naxis2=naxis2,include_path=True)

#print('copy out chosen frames to a new folder')
#import shutil
#for bias in calibrated_biases:
#    shutil.copy(bias,outfolder)


darks = (raw_data.summary['imagetyp'] == 'Dark Frame') & (raw_data.summary['naxis1'] == naxis1)
dark_times = set(raw_data.summary['exptime'][darks])
print('  check dark exposure times\n', dark_times)


print('create master dark, if memomry overflow, please shut all the other applications!')
for exp_time in sorted(dark_times):
    combined_dark = ccdp.combine(calibrated_darks,
                                 method='average',
                                 unit=u.adu,
                                 sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                 sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std,
                                 mem_limit=350e16
                                )

    combined_dark.meta['combined'] = True
    dark_file_name = 'combined_dark_{:6.3f}.fit'.format(exp_time)
    combined_dark.write(reduced_path / dark_file_name)

