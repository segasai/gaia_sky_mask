import numpy as np
import os
import astropy.units as auni
import astropy.coordinates as acoo
import sqlutilpy
import astropy.wcs as pywcs
import healpy
from cffi import FFI
import _filler
import astropy.io.fits as pyfits
import multiprocessing as mp
import pickle
from getmask import getwcs, getwcs0
ffi = FFI()
import fitsio

class si:
    DB = None
    nsideDB= None

def getdb(nsidedb):
    D = pickle.load(open('blobdb_%d.pkl' % nsidedb, 'rb'))
    hpx_index = {}
    for fname, v in D.items():
        v['brickname']
        for hp in v['hpx']:
            if hp not in hpx_index:
                hpx_index[hp] = []
            hpx_index[hp].append((fname, v['brickname'], v['hemi']))
    return hpx_index


def getblobs(hpx):
    pass


def fill_me(ima, xs, ys, binning=4):
    ima = np.require(ima, dtype=np.int32, requirements='AOC')
    npix = ima.shape[0]
    args = {'xs': xs, 'ys': ys}
    nstars = len(xs)
    Cargs = {}
    Xargs = {}
    for k in args.keys():
        Xargs[k] = np.require(args[k],
                              dtype=np.float32,
                              requirements=["A", "O", "C"])
        Cargs[k] = ffi.cast('float*', _filler.ffi.from_buffer(Xargs[k]))
    Cargs['ima'] = ffi.cast('int*', _filler.ffi.from_buffer(ima))
    _filler.lib.procima_pix(Cargs['ima'], npix, Cargs['xs'], Cargs['ys'],
                            nstars, binning)


def proc_ima(f, brickname, hemi):
    hdr = fitsio.read_header(f)
    dat = fitsio.read(f)
    wc = pywcs.WCS(hdr)
    xind = dat > -1
    nx = dat.shape[0]
    #ygrid, xgrid = np.mgrid[:nx, :nx]
    y,x=np.nonzero(xind)
    C = wc.pixel_to_world(x,y)#xgrid[xind], ygrid[xind])
    return C.ra.deg, C.dec.deg


def getmask_hpx(nside, hpx, ofname, DB=None, nsideDB=None):
    if DB is None:
        DB, nsideDB = si.DB, si.nsideDB
    ra0, dec0 = healpy.pix2ang(nside, hpx, lonlat=True, nest=True)
    step_asec = 1
    binning = 4
    angle, npix = getwcs(nside, hpx, step_asec=step_asec)
    wcdict = getwcs0(ra0, dec0, (npix + 1) / 2, (npix + 1) / 2, angle,
                     step_asec / 3600.)
    ima = np.zeros((npix, npix), dtype=np.int32)
    xgrid, ygrid = np.mgrid[:npix, :npix]
    wc = pywcs.WCS(wcdict)
    C = wc.pixel_to_world(xgrid, ygrid)
    del xgrid,ygrid
    hpxopt = healpy.ang2pix(nsideDB,
                            C.ra.deg,
                            C.dec.deg,
                            nest=True,
                            lonlat=True)
    del C
    hpxopt = np.unique(hpxopt)
    fs = set()
    for curhp in hpxopt:
        if curhp in DB:
            for curf in DB[curhp]:
                fs.add(curf)

    for curf, curbrick, curhemi in fs:
        curra, curdec = proc_ima(curf, curbrick, curhemi)
        cury, curx = wc.world_to_pixel(
            acoo.SkyCoord(ra=curra * auni.deg, dec=curdec * auni.deg))
        #curx,cury=[np.round(_).astype(int) for _ in [curx,cury]]
        xind = (curx >= -.5) & (cury < npix) & (curx >= -.5) & (cury < npix)
        fill_me(ima, curx[xind], cury[xind], binning)

    header = pyfits.Header()
    for k, v in wcdict.items():
        header[k] = v
    header['IMAGEW'] = npix
    header['IMAGEH'] = npix
    pyfits.PrimaryHDU(data=ima, header=header).writeto(ofname, overwrite=True)


def doall(nthreads=32):
    nside = 64
    res = []
    nsideDB = 256
    DB = getdb(nsideDB)
    si.DB = DB
    si.nsideDB=nsideDB
    pool = mp.Pool(nthreads, maxtasksperchild=4)
    for hpx in range(12 * nside * nside):
        ra,dec = healpy.pix2ang(nside, hpx, lonlat=True, nest=True)
        if dec<-33:
            continue
        ofname = 'data/ls_skymap-%05d.fits.gz' % hpx
        if os.path.exists(ofname):
            # do not overwrite
            continue 
        res.append((hpx,
                    pool.apply_async(
                        getmask_hpx,
                        (nside, hpx, ofname))))
    for r in res:
        hpx, r = r
        r.get()
        print(hpx)
