#!/usr/bin/env python

import os
import string
import itertools
import numpy as np

import parmap
import folders
import settings
import imfuncs

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord

from astropy.utils.exceptions import AstropyWarning
import warnings
warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', category=AstropyWarning)


# #### To start, download all available data from the VLA (Very Large Array) FIRST (Faint Images of the Radio Sky at Twenty-Centimeters) survey (http://sundog.stsci.edu/).
# #### These images are available through anonymous ftp from ftp://archive.stsci.edu/pub/vla_first/data and comprise a total size of about 340 GB.
# #### Set the folder where the data is downloaded in the folders.firstdata_folder variable
# #### This module consequently selects the last version of every available FIRST image, gathers and constructs a table with important necessary of these images.
# #### Most importantly the ra/dec covered by all useable first images as well as the RMS as estimated by the survey itself.
# #### Also do a new estimate for RMS from multiple areas in a FIRST image (needed later for statistics).
# #### This module can be run independently to create the table, but is also called from main module (first.py)
# #### (v1.0.dev1) PvO 23-01-2017


def makenewimtab(imchunksize=100):  # Main routine to start table making
    """
    Initiate creation if large table with necessary image information (to produce cutouts at later stage).
    Uses parallel. per core to gather rows.
    """
    IMSUSE = getlastims()  # Get last version of each FIRST image
    print 'Getting info on', len(IMSUSE), 'FIRST images'
    ROWS = []
    # Split image list into chunks to save memory /  avoid waiting cores at end of list
    IMCHUNKS = [IMSUSE[x: x + imchunksize] for x in range(0, len(IMSUSE), imchunksize)]
    for i in range(len(IMCHUNKS)):
        imchunk = IMCHUNKS[i]
        print 'doing chunk', i, '/', len(IMCHUNKS) - 1, 'length', len(imchunk)
        NEWROWS = parmap.starmap(getiminfo, zip(imchunk))
        ROWS = ROWS + NEWROWS
        i += 1
    return ROWS


def getiminfo(im, n_rmscutouts=10):  # Routines to get RA/DEC coverage from a FIRST image and a good estimate of its the RMS
    """
    Go through a FIRST image and determine its RA/DEC coverage from FITS headers
    Return coordinates between [0,360] degrees
    Also estimates RMS using getfirstrms()
    Eeturns rows
    """
    imfile = folders.firstdata_folder + im[0:5] + '/' + im
    hdulist = fits.open(imfile, memmap=None)
    # fitsio.read_header takes 2x longer for 100 ims (even tho memmap = False can speed it up slightly)
    imhdr = hdulist[0].header
    rac, decc = imhdr['CRVAL1'], imhdr['CRVAL2']
    w = WCS(imhdr).celestial
    # rac, decc = w.wcs.crval[0], w.wcs.crval[1]  # Also possible but need verification
    ras, decs = w.all_pix2world([1, imhdr.get('NAXIS1')], [1, imhdr.get('NAXIS2')], 1)

    rac = rac % 360.
    ras = ras % 360.
    ramin, ramax = cutra(ras)
    decmin, decmax = min(decs), max(decs)

    # See if RMs estimate is available from header
    try:
        rmshdr = float([x.split()[4] for x in imhdr['HISTORY'][0:20] if 'RMSNOISE' in x][0])
    except:
        rmshdr = -1.

    # imdata = fitsio.read(imfile, memmap=False)  # data read DOES seem faster in fitsio...due to using getdata()..highly ineff acc to manual
    imdata = hdulist[0].data
    rms = getfirstrms(imdata, rac, decc, w, boxdiam=settings.setboxdiam, n_rmscutouts=n_rmscutouts)
    hdulist.close()
    # del imdata, imhdr

    return rac, decc, ramin, ramax, decmin, decmax, rms, rmshdr, im


# easiest/best to call getfirstrms from getiminfo (to avoid having to supply rac,decc,w but rather just one coordinate)
def getfirstrms(imdata, rac, decc, w, boxdiam=settings.setboxdiam, n_rmscutouts=10):
    """
    Asses RMS at chosen coordinates by sampling an median-averaging local RMS at coordinates steps of 0.05 degree away from those coordinates.
    The average is of n_rmscutouts images with sides of boxdiam length.)
    Supply both RA/DEC as well as WCS coordinates
    """
    tstra = [rac, rac + 0.05, rac - 0.05, rac + 0.1, rac - 0.1, rac + 0.15, rac - 0.15, rac + 0.2, rac - 0.2]
    tstdec = [decc, decc + 0.05, decc - 0.05, decc + 0.1, decc - 0.1, decc + 0.15, decc - 0.15]
    tstc = []
    for r in itertools.product(tstra, tstdec):
        # if r != (rac,decc):  # don't have to remove centre ra/dec: in survey as no expected source in centre of fov if any FIRST images
        crd = SkyCoord(r[0], r[1], unit='deg')
        tstc.append(crd)
    if n_rmscutouts > len(tstc):
        print "(***) Warning: increase number of testra, tstdec. Only", len(tstc), "available at the moment."
        print "(***) Setting n_rmscutouts to max"
        n_rmscutouts = len(tstc)
        print n_rmscutouts
    rmsarr = []
    i = 0
    cutoutsize = (boxdiam, boxdiam)
    while len(rmsarr) < n_rmscutouts:
        c = tstc[i]
        cutout = Cutout2D(imdata[0][0], c, size=cutoutsize, wcs=w)  # , mode='strict')  # copy = True)
        if imfuncs.imcheck(cutout.data, boxdiam) == 0:
            rmsarr.append(cutout.data.std())
        else:
            # print 'im not useable'
            pass
        i += 1
    return np.median(rmsarr)


# #### Helper functions


def getlastims():
    """
    Make list to ensure only use the last version of each FIRST image
    Make list of folders from main FIRST db and find last version of each FIRST image
    Uses parallel. per FIRST db folder.
    """
    print "Making list of last version of each FIRST image"
    thedir = folders.firstdata_folder
    FOLDERS = [name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name))]  # ensures subfolders are folders and not e.g. hidden files created by OS
    res = parmap.starmap(lastverperfolder, zip(FOLDERS))
    IMS = list(itertools.chain(*res))
    print "Done making list"
    return IMS


def lastverperfolder(folder):
    """
    Find last version of each FIRST image, within each FIRST db subfolder
    First create dicts to associate last letter in filename root with respective value
    Special case for three images that do not have version letter in their root file name
    """
    lettonumdic = dict(zip(string.letters, [ord(c) % 32 for c in string.letters]))
    lettonumdic['.'] = 27
    numtoletdic = dict(zip([ord(c) % 32 for c in string.letters], string.letters[0:26]))
    numtoletdic[27] = ''

    filetab = Table(rows=[(x[:11], lettonumdic[x[11]]) for x in os.listdir(folders.firstdata_folder + folder)])
    IMCOORDS = list(set(filetab['col0']))
    IMS = []
    for imcoord in IMCOORDS:
        tstval = filetab[filetab['col0'] == imcoord]
        maxval = max(tstval['col1'])
        fname = imcoord + numtoletdic[maxval] + '.fits'
        IMS.append(fname)
    return IMS


def cutra(ras):
    """
    Make sure RA interval lies between [0, 360] and decide how to handle values close to boundaries
    """
    if abs(ras[0] - ras[1]) > 180:  # Correction only necessary if two ra values >180 degrees apart
        ratst = []
        for ra in ras:
            if ra % 360 > 180:
                ratst.append(abs(360 - ra % 360))
            elif ra % 360 <= 180:
                ratst.append(ra % 360)
        replaceidx = np.argmin(ratst)
        if replaceidx == 0:
            keepidx = 1
        else:
            keepidx = 0
        if ras[replaceidx] < ras[keepidx]:
            ras[replaceidx] = 360.
        else:
            ras[replaceidx] = 0.
    ramin, ramax = min(ras), max(ras)
    return ramin, ramax


if __name__ == "__main__":
    ROWS = makenewimtab(64)  # Get wanted info (see namear below) on all last images
    namear = ['ra', 'dec', 'ramin', 'ramax', 'decmin', 'decmax', 'rms', 'rmshdr', 'im']
    imtab = Table(rows=ROWS, names=namear)
    imtab.write('imtab8_tmp.fits', overwrite=True)
