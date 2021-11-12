import glob
import healpy
import astropy.wcs as pywcs
import numpy as np
import fitsio
import multiprocessing as mp
import pickle
import re

nside = 256


def gethpx(f):
    # hdr= pyfits.getheader(f)
    hdr = fitsio.read_header(f)
    # D=pyfits.getdata(f)
    nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
    step = 10
    assert (nx % step == 0)
    assert (ny % step == 0)
    xgrid, ygrid = np.mgrid[0:nx:step, 0:ny:step]
    wc = pywcs.WCS(hdr)
    C = wc.pixel_to_world(xgrid, ygrid)
    hpx = healpy.ang2pix(nside, C.ra.deg, C.dec.deg, nest=True, lonlat=True)
    return np.unique(hpx)


if __name__ == '__main__':
    prefix = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/'
    fs = glob.glob(prefix + '/north/metrics/*/blobs*.fits.gz')
    fs1 = glob.glob(prefix + '/south/metrics/*/blobs*.fits.gz')

    fs = fs + fs1

    pool = mp.Pool(64)
    hpxs = pool.map(gethpx, fs)
    hemis = []
    bricknames = []
    ret = {}
    for f0, hpx in zip(fs, hpxs):
        f = f0.split('/')
        hemis.append(f[-4])
        M = re.match('blobs-(.*).fits.gz', f[-1])
        bricknames.append(M.group(1))
        ret[f0] = dict(brickname=M.group(1), hemi=hemis[-1], hpx=hpx)
    with open('blobdb_%d.pkl' % nside, 'wb') as fp:
        pickle.dump(ret, fp)
