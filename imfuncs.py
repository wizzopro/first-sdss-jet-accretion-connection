#!/usr/bin/env python
import numpy as np
from astropy.io import fits


# ###################
# HELPER FUNCTIONS for image cutouts
# Version 1.0.dev1     P v Oers 10-3-2017
# ###################


def makefitsimfromimagehdu(filename, imagehdu):
    """
    Write an imageHDU to a file using astropy fits package
    Overwrites old files
    """
    hdu = fits.HDUList()
    hdu.append(imagehdu)
    hdu.writeto(filename, clobber=True)
    # print 'wrote',filename
    hdu.close()


def imcheck(imdata, boxdiam, repeats=25):
    """
    Image quality checks for cutouts. Cutouts report an error if:
    - There is a NAN in the cutout
    - There are more than nrepeats pixels with EXACTLY the same flux value
    - The dimensions of the cutout are not square / do not match the required dimensions (probably cutout close to edge)
    """
    if (np.isnan(np.sum(imdata)) is True or np.unique(imdata, return_counts=True)[1].max() > repeats or imdata.shape != (boxdiam, boxdiam)):
        error = 1
    else:
        error = 0
    return error
