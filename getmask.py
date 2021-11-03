import numpy as np
import astropy.units as auni
import astropy.coordinates as acoo
import sqlutilpy
import astropy.wcs as pywcs
import healpy
from cffi import FFI
import _filler
import astropy.io.fits as pyfits
import multiprocessing as mp

ffi = FFI()
WSDB_HOST = open('WSDB_HOST').read()


def getrad(Gmag):
    return 600. * 1.4**(-Gmag)


def getangle(x, y):
    pos = np.argmax(y)
    pos1 = (pos + 3) % 4
    angle = np.arctan2(y[pos] - y[pos1], x[pos] - x[pos1])
    return float(angle)


def getwcs0(ra0, dec0, crpix1, crpix2, angle, step):
    return dict(CRVAL1=ra0,
                CRVAL2=dec0,
                CRPIX1=crpix1,
                CRPIX2=crpix2,
                CD1_1=np.cos(angle) * step,
                CD2_2=np.cos(angle) * step,
                CD1_2=-np.sin(angle) * step,
                CD2_1=np.sin(angle) * step,
                CTYPE1='RA---TAN',
                CTYPE2='DEC--TAN')


def getwcs(nside, hpx, step_asec=1):
    step = step_asec / 3600
    ra0, dec0 = healpy.pix2ang(nside, hpx, lonlat=True, nest=True)
    xyz = healpy.boundaries(nside, hpx, step=1, nest=True)
    edges = np.array([healpy.vec2ang(_, lonlat=True) for _ in xyz.T])
    wc = pywcs.WCS(getwcs0(ra0, dec0, 1, 1, 0, 1))
    x, y = wc.world_to_pixel(
        acoo.SkyCoord(ra=edges[:, 0] * auni.deg, dec=edges[:, 1] * auni.deg))
    angle = getangle(x, y)
    wc = pywcs.WCS(getwcs0(ra0, dec0, 1, 1, angle, step))
    x, y = wc.world_to_pixel(
        acoo.SkyCoord(ra=edges[:, 0] * auni.deg, dec=edges[:, 1] * auni.deg))
    expand = 1.05
    npix = 2 * int(np.ceil(expand * max(np.max(np.abs(x)), np.max(np.abs(y)))))
    # ensure rectangular
    return angle, npix


def getmask_hpx(nside, hpx, ofname):
    ra0, dec0 = healpy.pix2ang(nside, hpx, lonlat=True, nest=True)
    step_asec = 1
    angle, npix = getwcs(nside, hpx, step_asec=step_asec)
    ima, wcdict = getmask(ra0, dec0, step_asec, npix, angle=angle)
    header = pyfits.Header()
    for k, v in wcdict.items():
        header[k] = v
    header['IMAGEW'] = npix
    header['IMAGEH'] = npix
    pyfits.PrimaryHDU(data=ima, header=header).writeto(ofname, overwrite=True)


def doall(nthreads=36):
    nside = 64
    pool = mp.Pool(nthreads)
    res = []
    for hpx in range(12 * nside * nside):
        res.append((hpx,
                    pool.apply_async(
                        getmask_hpx,
                        (nside, hpx, 'data/gaia_skymap-%05d.fits.gz' % hpx))))
    for r in res:
        hpx, r = r
        r.get()
        print(hpx)


def getmask(ra0, dec0, step_asec, npix, angle=0, binning=4, epoch=2022):
    """
    Return the integer mask and WCS object for a gaia based sky mask
    Arguments
    ra0,dec0 center
    pixel size in arcsec
    number of pixels
    """
    minrad_asec, maxrad_asec = 1.5, 1000

    step = step_asec / 3600.
    wcdict = getwcs0(ra0, dec0, (npix + 1) / 2, (npix + 1) / 2, angle, step)
    assert (binning in [1, 2, 4])
    wc = pywcs.WCS(wcdict)
    pad = maxrad_asec / 3600.  # padding in degrees for bright stars
    maxrad = npix * np.sqrt(2) * step + pad
    ra, dec, pmra, pmdec, refep, g = sqlutilpy.get(
        '''select ra,dec,pmra,pmdec,ref_epoch,phot_g_mean_mag from gaia_edr3.gaia_source
    where q3c_radial_query(ra,dec,%f,%f,%f)''' % (ra0, dec0, maxrad),
        host=WSDB_HOST)
    bad = ~np.isfinite(pmra)
    pmra[bad] = 0
    pmdec[bad] = 0
    ra = ra + pmra / 3600e3 / np.cos(np.deg2rad(dec)) * (epoch - refep)
    dec = dec + pmdec / 3600e3 * (epoch - refep)
    rads = np.clip(getrad(g), minrad_asec, maxrad_asec)
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
                        Cargs['rads'], nstars, binning)

    return ima, wcdict
