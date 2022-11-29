# REFERNCE https://mwcraig.github.io/ccd-as-book/06-00-Reducing-science-images.html



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
from matplotlib import pyplot as plt
import numpy as np
from astropy import units as u
import ccdproc as ccdp
#from convenience_functions import show_image


IN_PATH = '/work/chuck/chen/obs/20190822'
IN_PATH = '/work/chuck/chen/obs/20190827'
IN_PATH = '/work/chuck/chen/obs/20190823'

outfolder = 'reduced'
reduced_path = Path(IN_PATH, outfolder)
reduced_path.mkdir(exist_ok=True)

science_imagetyp = 'Light Frame'
flat_imagetyp = 'Flat Field'
exposure = 'exposure'
naxis1 = 2749 # the size we needed
naxis2 = 2199 
exptime_tolerance = 300 # If have better choice, reduce it!

print('read science images')
raw_data = ccdp.ImageFileCollection(IN_PATH)
lights = raw_data.summary[(raw_data.summary['imagetyp'] == science_imagetyp) & (raw_data.summary['naxis1'] == naxis1)]
#print('  Science Images:\n', lights['date-obs', 'file', 'object', 'filter', exposure])

print('read calibrated images, make dictionary according to Filter')
reduced_data = ccdp.ImageFileCollection(reduced_path)
combo_calibs = reduced_data.summary[reduced_data.summary['combined'].filled(False).astype('bool')]
print('  Calibrates Images:\n', combo_calibs['date-obs', 'file', 'imagetyp', 'filter', exposure])

combined_darks = {ccd.header[exposure]: ccd for ccd in reduced_data.ccds(imagetyp='Dark Frame', combined=True)}
combined_flats = {ccd.header['FILTER']: ccd for ccd in reduced_data.ccds(imagetyp=flat_imagetyp, combined=True)}
combined_bias = [ccd for ccd in reduced_data.ccds(imagetyp='Bias Frame', combined=True)][0]

print('Reducing science images by filter and exptime in Loop')
i = 1
for ccd, file_name in raw_data.ccds(imagetyp=science_imagetyp,naxis1=naxis1,ccd_kwargs={'unit': 'adu'},return_fname=True):
    print('  reducing', file_name, i,'/',len(lights))
    i+=1
    #ccd = ccdp.trim_image(ccd[:, :4096])
    '''    
    Note that the first argument in the remainder of the ccdproc calls is
    the *ccd* image, so that the calibration steps are cumulative.
    '''
    closest_dark = find_nearest_dark_exposure(ccd, combined_darks.keys(), tolerance=exptime_tolerance)
    if (closest_dark-ccd.header['exptime']>100):
      #print('    need substract bias due to master dark exptime differs a lot with target')
      ccd = ccdp.subtract_bias(ccd, combined_bias)
    ccd = ccdp.subtract_dark(ccd, combined_darks[closest_dark], 
                                 exposure_time='exptime', exposure_unit=u.second)
    master_flat = combined_flats[ccd.header['FILTER']]
    ccd = ccdp.flat_correct(ccd, master_flat)
    ccd.meta['combined'] = True
    ccd_name = f'light_{file_name}'
    ccd.write(reduced_path / ccd_name)

