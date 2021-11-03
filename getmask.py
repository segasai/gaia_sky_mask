import numpy as np
import astropy.units as auni
import astropy.coordinates as acoo
import sqlutilpy
import astropy.wcs as pywcs

from cffi import FFI
import _filler

ffi = FFI()
WSDB_HOST = open('WSDB_HOST').read()


def getmask(ra0, dec0, step_asec, npix):
    """
    Return the integer mask and WCS object for a gaia based sky mask
    Arguments
    ra0,dec0 center
    pixel size in arcsec
    number of pixels
    """
    step = step_asec / 3600
    wcdict = dict(CRVAL1=ra0,
                  CRVAL2=dec0,
                  CRPIX1=(npix + 1) / 2,
                  CRPIX2=(npix + 1) / 2,
                  CD1_1=step,
                  CD2_2=step,
                  CD1_2=0,
                  CD2_1=0,
                  CTYPE1='RA---TAN',
                  CTYPE2='DEC--TAN')

    wc = pywcs.WCS(wcdict)
    '''ra, dec = healpy.pix2ang(nsidex, iis, lonlat=True, nest=True)'''
    pad = 0.2  # padding for bright stars
    maxrad = npix * np.sqrt(2) * step + pad
    ra, dec, g = sqlutilpy.get(
        '''select ra,dec,phot_g_mean_mag from gaia_edr3.gaia_source
    where q3c_radial_query(ra,dec,%f,%f,%f)''' % (ra0, dec0, maxrad),
        host=WSDB_HOST)
    # hpx = healpy.ang2pix(nsidex, ra, dec, lonlat=True, nest=True)
    mags = np.minimum(g, 20)

    rads = np.maximum(630. * 1.396**(-mags), 1)
    rads = rads / step_asec  # in pixels now
    nstars = len(ra)

    ys, xs = wc.world_to_pixel(
        acoo.SkyCoord(ra=ra * auni.deg, dec=dec * auni.deg))
    ima = np.zeros((npix, npix), dtype=np.int32)
    ima = np.require(ima, dtype=np.int32, requirements='AOC')

    args = {'xs': xs, 'ys': ys, 'rads': rads}

    Cargs = {}
    Xargs = {}
    for k in args.keys():
        Xargs[k] = np.require(args[k],
                              dtype=np.float32,
                              requirements=["A", "O", "C"])
        Cargs[k] = ffi.cast('float*', _filler.ffi.from_buffer(Xargs[k]))
    Cargs['ima'] = ffi.cast('int*', _filler.ffi.from_buffer(ima))
    _filler.lib.procima(Cargs['ima'], npix, Cargs['xs'], Cargs['ys'],
                        Cargs['rads'], nstars)

    return ima, wc
