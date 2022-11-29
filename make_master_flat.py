# REFERENCE: https://mwcraig.github.io/ccd-as-book/05-03-Calibrating-the-flats.html
# When dark exposure is different with flat frames, need to reduce the bias.


def find_nearest_dark_exposure(image, dark_exposure_times, tolerance=0.5):
    """
    Find the nearest exposure time of a dark frame to the exposure time of the image,
    raising an error if the difference in exposure time is more than tolerance.
    
    Parameters
    ----------
    image : astropy.nddata.CCDData
        Image for which a matching dark is needed.
    dark_exposure_times : list
        Exposure times for which there are darks.
    tolerance : float or ``None``, optional
        Maximum difference, in seconds, between the image and the closest dark. Set
        to ``None`` to skip the tolerance test.
    Returns
    -------
    float
        Closest dark exposure time to the image.
    """

    dark_exposures = np.array(list(dark_exposure_times))
    idx = np.argmin(np.abs(dark_exposures - image.header['exptime']))
    closest_dark_exposure = dark_exposures[idx]

    if (tolerance is not None and 
        np.abs(image.header['exptime'] - closest_dark_exposure) > tolerance):
        raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
                           'time {}.'.format(closest_dark_exposure, a_flat.header['exptime']))
    return closest_dark_exposure


from pathlib import Path
from astropy import units as u
from astropy.nddata import CCDData
import ccdproc as ccdp
from matplotlib import pyplot as plt
import numpy as np


IN_PATH = '/work/chuck/chen/obs/20190822'
IN_PATH = '/work/chuck/chen/obs/20190827'
IN_PATH = '/work/chuck/chen/obs/20190823'

outfolder = 'reduced'
reduced_path = Path(IN_PATH, outfolder)
reduced_path.mkdir(exist_ok=True)

flat_image_type = 'Flat Field'
naxis1 = 2749 # the size we needed
naxis2 = 2199 
exptime_tolerance = 300 # If have better choice, reduce it!


print('read dark data and making dictionary according to exptime')
dark_path = Path(IN_PATH, outfolder)
dark_reduced = ccdp.ImageFileCollection(dark_path)
combined_darks = {ccd.header['exptime']: ccd for ccd in dark_reduced.ccds(imagetyp='Dark Frame', combined=True)}
combined_dark_files = dark_reduced.files_filtered(imagetyp='Dark Frame', combined=True)

print('read bias data')
bias_path = Path(IN_PATH, outfolder)
bias_reduced = ccdp.ImageFileCollection(bias_path)
combined_bias = list(bias_reduced.ccds(combined=True, imagetyp='Bias Frame'))[0]

print('read raw data')
raw_data_directory = Path(IN_PATH)
raw_data = ccdp.ImageFileCollection(raw_data_directory)

print('check master dark exposure time')
n_combined_dark = len(combined_dark_files)
n_dark_expected = 1
expected_exposure_times = set([300])

if n_combined_dark < n_dark_expected:
    raise RuntimeError('One or more combined dark is missing. Please re-run the dark notebook.')
elif n_combined_dark > n_dark_expected:
    raise RuntimeError('There are more combined dark frames than expected.')
    
actual_exposure_times = set(h['exptime'] for h in dark_reduced.headers(imagetyp='Dark Frame', combined=True))

if (expected_exposure_times - actual_exposure_times):
    raise RuntimeError('Encountered unexpected exposure time in combined darks. '
                       'The unexpected times are {}'.format(actual_exposure_times - expected_exposure_times))

print('check flat frames')
for a_flat, fname in raw_data.hdus(imagetyp=flat_image_type,naxis1=naxis1, return_fname=True):
    print(f'  {fname}',a_flat.shape,'fileter',a_flat.header['FILTER'], 'exposure', a_flat.header['EXPOSURE'], 'standard deviation ', a_flat.data.std())
calibrated_darks = raw_data.files_filtered(imagetyp='Dark Frame',naxis1=naxis1,naxis2=naxis2,include_path=True)

flats = (raw_data.summary['imagetyp'] == flat_image_type) & (raw_data.summary['naxis1'] == naxis1)
print('checking flats exposure times:\n', set(raw_data.summary['exptime'][flats]))


print('calibrite all the flats in the folder')
for ccd, file_name in raw_data.ccds(imagetyp=flat_image_type, naxis1=naxis1, naxis2=naxis2,  # Get the flats frames of chosen size
                                   return_fname=True,                                        # Provide the file name too.
                                   ccd_kwargs={'unit': 'adu'}, # CCDData requires a unit for the image if it is not in the header
                                  ): 
    print('  calibrite',file_name)
    #ccd = ccdp.trim_image(ccd[:, :4096])                                                    # Trim the overscan
    closest_dark = find_nearest_dark_exposure(ccd, actual_exposure_times, tolerance=exptime_tolerance) # Find the correct dark exposure
    if (closest_dark-ccd.header['exptime']>100):
      print('    need substract bias due to master dark exptime differs a lot with target')
      combined_bias = list(bias_reduced.ccds(combined=True, imagetyp='Bias Frame'))[0]
      ccd = ccdp.subtract_bias(ccd, combined_bias)
    ccd = ccdp.subtract_dark(ccd, combined_darks[closest_dark],
                             exposure_time='exptime', exposure_unit=u.second, scale=True)    # Subtract the dark current 
    ccd.meta['combined'] = True
    ccd_name = f'flat_{file_name}'
    ccd.write(reduced_path / ccd_name)                                                      # Save the result
