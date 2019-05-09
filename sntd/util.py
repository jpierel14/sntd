#!/Users/jpierel/anaconda3/envs/astro2/bin python2

import os,sncosmo,glob
from astropy.io import ascii
import numpy as np
from collections import OrderedDict as odict
from scipy.interpolate import splrep,splev
from copy import copy


__current_dir__=os.path.abspath(os.getcwd())
__dir__=os.path.abspath(os.path.dirname(__file__))

NORMAL = 0    # use python zip libraries
PROCESS = 1   # use (zcat, gzip) or (bzcat, bzip2)
PARALLEL = 2  # (pigz -dc, pigz) or (pbzip2 -dc, pbzip2)

__all__=['flux_to_mag','_cast_str','_get_default_prop_name','_isfloat','anyOpen','_props','_findMax','_findMin',
         '_guess_time_delays','_guess_magnifications','__dir__','load_example_data']
_props=odict([
    ('time',{'mjd', 'mjdobs', 'jd', 'time', 'date', 'mjd_obs','mhjd','jds'}),
    ('band',{'filter', 'band', 'flt', 'bandpass'}),
    ('flux',{'flux', 'f','fluxes'}),
    ('fluxerr',{'flux_error', 'fluxerr', 'fluxerror', 'fe', 'flux_err','fluxerrs'}),
    ('zp',{'zero_point','zp', 'zpt', 'zeropoint'}),
    ('zpsys',{'zpsys', 'magsys', 'zpmagsys'}),
    ('mag',{'mag','magnitude','mags'}),
    ('magerr',{'magerr','magerror','magnitudeerror','magnitudeerr','magerrs'})
])


def load_example_data():
    example_files=glob.glob(os.path.join(__dir__,'data','examples','*.dat'))
    return(ascii.read(example_files[0]),ascii.read(example_files[1]))

def _guess_magnifications(curves,referenceImage):
    """Guess t0 and amplitude of the model based on the data.

    Assumes the data has been standardized."""
    ref=None
    mags=dict([])
    all_ims=[referenceImage]
    for k in [x for x in curves.images.keys() if x!=referenceImage]:
        all_ims.append(k)
    for k in all_ims:
        if not curves.images[k].fits:
            bestRatio=-np.inf
            bestBand=None
            for b in np.unique(curves.images[k].table['band']):
                ratio=np.abs(np.max(curves.images[k].table['flux'][curves.images[k].table['band']==b])/np.min(curves.images[k].table['flux'][curves.images[k].table['band']==b]))
                if ratio>bestRatio:
                    bestRatio=ratio
                    bestBand=b
            maxTime,maxValue=_findMax(curves.images[k].table['time'][curves.images[k].table['band']==b],curves.images[k].table['flux'][curves.images[k].table['band']==b])


            if k==referenceImage:
                ref=copy(maxValue)
            mags[k]=maxValue/ref

        else:
            if 'x0' in curves.images[k].fits.res.param_names:
                amplitude='x0'
            else:
                amplitude='amplitude'
            mag=curves.images[k].fits.model.get(amplitude)
            if k==referenceImage:
                ref=copy(mag)
            mags[k]=mag/ref
    return mags


def _guess_time_delays(curves,referenceImage):
    #tds=colorFit(curves)
    #if tds:
    #    return tds
    ref=None
    tds=dict([])
    all_ims=[referenceImage]
    for k in [x for x in curves.images.keys() if x!=referenceImage]:
        all_ims.append(k)
    for k in all_ims:
        if not curves.images[k].fits:
            maxValue,maxFlux=_findMax(curves.images[k].table['time'],curves.images[k].table['flux'])
            if k==referenceImage:
                ref=copy(maxValue)
            tds[k]=maxValue-ref
        else:
            t0=curves.images[k].fits.model.get('t0')
            if k==referenceImage:
                ref=copy(t0)
            tds[k]=t0-ref
    return tds

def _findMax(time,curve):
    #TODO check edge cases
    t0=np.where(curve==np.max(curve))[0][0]
    if t0==0:
        #return(time[0])
        return (None,None)

    elif t0==len(time)-1:

        #return(time[-1])
        return (None,None)

    else:
        fit=splrep(time[t0-1:t0+2],curve[t0-1:t0+2],k=2)

    interptime=np.linspace(time[t0-1],time[t0+1],100)
    flux=splev(interptime,fit)
    return(interptime[flux==np.max(flux)],np.max(flux))


def _findMin(time,curve):
    #TODO check edge cases
    t0=np.where(curve==np.min(curve))[0][0]

    if t0==0:
        #return(time[0])
        return (None,None)

    elif t0==len(time)-1:
        #return(time[-1])
        return (None,None)

    else:
        fit=splrep(time[t0-1:t0+2],curve[t0-1:t0+2],k=2)

    interptime=np.linspace(time[t0-1],time[t0+1],100)
    flux=splev(interptime,fit)
    return(interptime[flux==np.min(flux)],np.max(flux))

def flux_to_mag(table,bandDict,zpsys='AB'):
    """
    Accepts an astropy table of flux data and does the conversion to mags (flux in ergs/s/cm^2/AA)

    :param table: Table containing flux, flux error error, and band columns
    :type table: astropy.Table
    :param bandDict: translates band to sncosmo.Bandpass object (i.e. 'U'-->bessellux)
    :type bandDict: dict
    :param zpsys: magnitude system
    :type zpsys: str,optional
    :returns: astropy.Table object with mag and mag error added
    """
    table=table[table['flux']>0]
    ms=sncosmo.get_magsystem(zpsys)
    table[_get_default_prop_name('mag')] = np.asarray(map(lambda x, y: ms.band_flux_to_mag(x,y), table[_get_default_prop_name('flux')],#/sncosmo.constants.HC_ERG_AA
                                                       np.array([bandDict[z] for z in table[_get_default_prop_name('band')]])))
    table[_get_default_prop_name('magerr')] = np.asarray(map(lambda x, y: 2.5 * np.log10(np.e) * y / x, table[_get_default_prop_name('flux')],
                                                          table[_get_default_prop_name('fluxerr')]))
    table=table[np.abs(table['mag']/table['magerr'])>5]
    return(table)

def _cast_str(s):
    try:
        return int(s)
    except:
        try:
            return float(s)
        except:
            return s.strip()


def _get_default_prop_name(prop):
    for key,value in _props.items():
        if {prop.lower()} & value:
            return key
    return prop


def _isfloat(value):
    try:
        float(value)
        return True
    except:
        return False

def anyOpen(filename, mode='r', buff=1024*1024, external=PARALLEL):
    if 'r' in mode and 'w' in mode:
        return None
    if filename.startswith('!'):
        import subprocess
        if 'r' in mode:
            return subprocess.Popen(filename[1:], shell=True, bufsize=buff,
                                    stdout=subprocess.PIPE).stdout
        elif 'w' in mode:
            return subprocess.Popen(filename[1:], shell=True, bufsize=buff,
                                    stdin=subprocess.PIPE).stdin
    elif filename.endswith('.bz2'):
        if external == NORMAL:
            import bz2
            return bz2.BZ2File(filename, mode, buff)
        elif external == PROCESS:
            if not which('bzip2'):
                return anyOpen(filename, mode, buff, NORMAL)
            if 'r' in mode:
                return anyOpen('!bzip2 -dc ' + filename, mode, buff)
            elif 'w' in mode:
                return anyOpen('!bzip2 >' + filename, mode, buff)
        elif external == PARALLEL:
            if not which('pbzip2'):
                return anyOpen(filename, mode, buff, PROCESS)
            if 'r' in mode:
                return anyOpen('!pbzip2 -dc ' + filename, mode, buff)
            elif 'w' in mode:
                return anyOpen('!pbzip2 >' + filename, mode, buff)
    elif filename.endswith('.gz'):
        if external == NORMAL:
            import gzip
            return gzip.GzipFile(filename, mode, buff)
        elif external == PROCESS:
            if not which('gzip'):
                return anyOpen(filename, mode, buff, NORMAL)
            if 'r' in mode:
                return anyOpen('!gzip -dc ' + filename, mode, buff)
            elif 'w' in mode:
                return anyOpen('!gzip >' + filename, mode, buff)
        elif external == PARALLEL:
            if not which('pigz'):
                return anyOpen(filename, mode, buff, PROCESS)
            if 'r' in mode:
                return anyOpen('!pigz -dc ' + filename, mode, buff)
            elif 'w' in mode:
                return anyOpen('!pigz >' + filename, mode, buff)
    elif filename.endswith('.xz'):
        if which('xz'):
            if 'r' in mode:
                return anyOpen('!xz -dc ' + filename, mode, buff)
            elif 'w' in mode:
                return anyOpen('!xz >' + filename, mode, buff)
    else:
        return open(filename, mode, buff)
    return None

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
