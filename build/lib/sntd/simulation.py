import sncosmo, datetime, sys, os, scipy
import numpy as np
from astropy.time import Time
from astropy.table import Table
from astropy.io import ascii
from copy import deepcopy

from .util import __dir__
from .io import curve,curveDict
from .ml import getDiffCurve

__all__=['createRandMultiplyImagedSN']
def _getAbsoluteDist():
    absolutes=ascii.read(os.path.join(__dir__,'sim','data','absolutes.ref'))
    total=float(np.sum(absolutes['N'][absolutes['type']!='Ia']))
    absDict=dict([])
    for row in absolutes:
        if row['type']=='Ia':
            frac=1
        else:
            frac=float(row['N'])/total
        absDict[row['type']]={'dist':(row['mean'],row['sigma']),'frac':frac}
    return(absDict)

def _getAbsFromDist(dist):
    mu,sigma=dist
    return(np.random.normal(mu,sigma))

def createRandMultiplyImagedSN(
        sourcename, snType, redshift, telescopename='telescope',
        objectName='object', timeDelayRange=(0,30), muRange=(1,20),
        numImages=2, cadence=5, epochs=50, bands=['F140W','F105W','F160W'],
        gain=1000., skynoiseRange=(1,5), mjdRange=None, zpsys='ab', zp=None,
        microlensing=False, cc_av=.9, ia_av=.3, dust='CCM89Dust',
        dust_frames=['rest', 'free'],
        dust_names=['host', 'lens'], dust_ebv=[None, None],
        minsnr=0.0, av_host=.3, z_lens=None,
        av_lens=None, scatter=True):
    """Generate a multiply-imaged SN light curve set, with time delays and
    magnifications drawn randomly from user-defined ranges.
    """
    time_delays = list(np.random.uniform(timeDelayRange[0], timeDelayRange[1],
                                         numImages))
    magnifications = list(np.random.uniform(muRange[0], muRange[1],
                                            numImages))
    return createMultiplyImagedSN(
        sourcename, snType, redshift, telescopename=telescopename,
        objectName='object', time_delays=time_delays,
        magnifications=magnifications,
        numImages=numImages, cadence=cadence, epochs=epochs, bands=bands,
        gain=gain, skynoiseRange=skynoiseRange, mjdRange=mjdRange,
        zpsys=zpsys, zp=zp, microlensing=microlensing, cc_av=cc_av,
        ia_av=ia_av, dust_model=dust, screenz=screenz, screenebv=screenebv,
        dust_frames=dust_frames, dust_names=dust_names, minsnr=minsnr,
        av_host=av_host, z_lens=z_lens, av_lens=av_lens, scatter=scatter)



def createMultiplyImagedSN(
        sourcename, snType, redshift, telescopename='telescope',
        objectName='object', time_delays=[0., 20.], magnifications=[1., 5.],
        numImages=2, cadence=5, epochs=50, bands=['F105W', 'F125W', 'F160W'],
        gain=1000., skynoiseRange=(1, 5), mjdRange=None, zpsys='ab', zp=None,
        microlensing_type=None, microlensing_params=[],
        dust_model='CCM89Dust', av_host=.3, av_lens=None,
        z_lens=None, minsnr=0.0, scatter=True):
    """Generate a multiply-imaged SN light curve set, with user-specified time
    delays and magnifications.

    :param microlensing_type: str  Name any of the valid sncosmo microlensing
    types. Currently ['AchromaticSplineMicrolensing',
    'ChromaticSplineMicrolensing', 'AchromaticMicrolensing']
    :param microlensing_params: list   Parameters for the microlensing effect
    -- see sncosmo models.py for details.  Currently, if using either of the
    spline-based "mock" lensing options, then this params list must give
    three values for [nanchor, sigmadm, nspl].  If the microlensing_type is
    'AchromaticMicrolensing' then we are reading a simulated microlensing
    data file, and this params list must give [filename, mag_type] where
    mag_type is 'multiply' or 'add'.
    """
    if not mjdRange:
        now=np.round(Time(datetime.datetime.now()).mjd,3)
        times=np.linspace(now,now+cadence*epochs,epochs)
    else:
        #times=np.linspace(mjdRange[0],mjdRange[-1],epochs)
        times=np.arange(mjdRange[0],min(mjdRange[0]+cadence*epochs,mjdRange[-1])+cadence,cadence)
    bandList=np.array([np.tile(b,len(times)) for b in bands]).flatten()
    ms=sncosmo.get_magsystem(zpsys)

    zpList=[ms.band_flux_to_mag(1,b) for b in bandList] if not zp else [zp for i in range(len(bandList))]

    obj=curveDict(telescopename=telescopename,object=objectName)
    obj.bands = set(bandList)
    obj.zpDict = dict([(bandList[i], zpList[i]) for i in range(len(bandList))])
    obj.zpsys = zpsys
    obstable = Table({'time':np.tile(times,len(bands)), 'band':bandList,
                      'zpsys':[zpsys.upper() for i in range(len(bandList))],
                      'zp':zpList,
                      'skynoise':np.random.uniform(
                          skynoiseRange[0],skynoiseRange[1],len(bandList)),
                      'gain':[gain for i in range(len(bandList))]})
    #print(ascii.write(obstable['band','time','gain','skynoise','zp','zpsys'],Writer=ascii.Latex,latexdict={'units':{'time':'Days','skynoise':'Counts'}}))
    absolutes=_getAbsoluteDist()

    # Set up the dust_model extinction effects in the host galaxy and lens plane
    # TODO allow additional dust screens, not in the host or lens plane?
    # TODO sample from a prior for host and lens-plane dust A_V?
    # TODO : allow different lens-plane dust_model for each image?
    R_V = 3.1  # TODO: allow user-specified alternate dust R_V
    RV_lens = R_V
    RV_host = R_V
    dust_frames = []
    dust_names = []
    dust_effect_list = []
    if dust_model and (av_lens or av_host):
        dust_effect = {'SFD98Map': sncosmo.SFD98Map,
                       'CCM89Dust': sncosmo.CCM89Dust,
                       'OD94Dust': sncosmo.OD94Dust,
                       'F99Dust': sncosmo.F99Dust}[dust_model]()
        if av_host:
            dust_frames.append('rest')
            dust_names.append('host')
            dust_effect_list.append(dust_effect)
        if av_lens:
            dust_frames.append('free')
            dust_names.append('lens')
            dust_effect_list.append(dust_effect)

    # The following is not needed, but may be resurrected when we allow user
    # to provide additional dust screens.
    #if not isinstance(dust_names, (list, tuple)):
    #    dust_names=[dust_names]
    #if not isinstance(dust_frames, (list, tuple)):
    #    dust_frames=[dust_frames]

    # The sncosmo Model is initially set up with only dust effects, because
    # as currently constructed, dust has the same effect on all images.
    # Microlensing effects are added separately for each SN image below.
    model=sncosmo.Model(source=sourcename, effects=dust_effect_list,
                        effect_names=dust_names, effect_frames=dust_frames)
    model.set(z=redshift)
    if snType in ['IIP','IIL','IIn']:
        absBand='bessellb'
    else:
        absBand='bessellr'
    model.set_source_peakabsmag(_getAbsFromDist(absolutes[snType]['dist']),
                                absBand, zpsys)
    # TODO: allow user to specify parameters like x1, c, t0 if desired.
    t0 = times[0] + .25 * (times[-1] - times[0])
    if snType=='Ia':
        x0=model.get('x0')
        params={'z':redshift, 't0':t0, 'x0':x0,
                'x1':np.random.normal(0.,1.), 'c':np.random.normal(0.,.1)}
    else:
        amp=model.get('amplitude')
        params={'z':redshift, 't0':t0, 'amplitude':amp}
    model.set(**params)
    if av_host:
        ebv_host = av_host/RV_host
        model.set(hostebv=ebv_host, hostr_v=RV_host)
    else:
        ebv_host = 0
    if av_lens:
        if z_lens is None:
            z_lens = redshift / 2.  # TODO : Warn user about missing z_lens
        ebv_lens = av_lens/RV_lens
        model.set(lensz=z_lens, lensebv=ebv_lens, lensr_v=RV_lens)
    else:
        ebv_lens = 0

    # Step through each of the multiple SN images, adding time delays,
    # macro magnifications, and microlensing effects.
    for imnum, td, mu in zip(range(numImages), time_delays, magnifications):
        # Make a separate model_i for each SN image, so that lensing effects
        # can be reflected in the model_i parameters and propagate correctly
        # into realize_lcs for flux uncertainties
        model_i = deepcopy(model)
        params_i = deepcopy(params)
        if snType=='Ia':
            params_i['x0'] *= mu
        else:
            params_i['amplitude'] *= mu
        params_i['t0'] += td

        if microlensing_type is not None:
            # add microlensing effect
            if 'spline' in microlensing_type.lower():
                # Initiate a spline-based mock ml effect (at this point,
                # a set of random splines is generated and the microlensing
                # magnification curve is fixed)
                nanchor, sigmadm, nspl = microlensing_params
                if microlensing_type.lower().startswith('achromatic'):
                    ml_spline_func = sncosmo.AchromaticSplineMicrolensing
                else :
                    ml_spline_func = sncosmo.ChromaticSplineMicrolensing
                ml_effect = ml_spline_func(nanchor=nanchor, sigmadm=sigmadm,
                                           nspl=nspl)
            else:
                # Load a microlensing simulation from a data file
                # TODO : guard against common user entry errors
                # TODO : allow random file selection from a directory
                ml_filename = microlensing_params[0]
                ml_mag_type = microlensing_params[1]
                ml_effect = sncosmo.AchromaticMicrolensing(
                    ml_filename, ml_mag_type)
            model_i.add_effect(ml_effect, 'microlensing', 'free')
            model_i.set(microlensingz=z_lens)
            params_i['microlensingz'] = z_lens # Redundant?
        else:
            ml_effect = None

        # Generate the simulated SN light curve observations, make a `curve`
        # object, and store the simulation metadata
        model_i.set(**params_i)
        #print(np.nanmin(model_i.bandmag('F125W',zpsys,np.arange(model_i.mintime(),model_i.maxtime(),1))))
        table_i = sncosmo.realize_lcs(
            obstable , model_i, [params_i],
            trim_observations=True, scatter=scatter,thresh=minsnr)
        while len(table_i)==0 or len(table_i[0])<numImages:
            table_i = sncosmo.realize_lcs(
                obstable , model_i, [params_i],
                trim_observations=True, scatter=scatter,thresh=minsnr)
        table_i=table_i[0]
        curve_i=curve(zp=zp,zpsys=zpsys)
        curve_i.object=None

        curve_i.table=table_i
        curve_i.bands=list(set(table_i['band']))
        curve_i.simMeta=table_i.meta
        curve_i.simMeta['model']=model_i
        curve_i.simMeta['hostebv']=ebv_host
        curve_i.simMeta['lensebv']=ebv_lens
        curve_i.simMeta['lensz']=z_lens
        curve_i.simMeta['mu']=mu
        curve_i.simMeta['td']=td
        curve_i.simMeta['microlensing'] = ml_effect
        curve_i.simMeta['microlensing_type'] = microlensing_type
        curve_i.simMeta['microlensing_params'] = microlensing_params
        obj.add_curve(curve_i)

        # #do my own scatter, sncosmo's is weird
        # for i in range(len(lc_i)):
        #     temp=np.random.normal(lc_i['flux'][i],lc_i['fluxerr'][i])
        #     # 1-sigma clipping ?!
        #     #while np.abs(temp)>(np.abs(lc_i['flux'][i])+np.abs(lc_i['fluxerr'][i])):
        #     #    temp=np.random.normal(lc_i['flux'][i],lc_i['fluxerr'][i])
        #     lc_i['flux'][i]=temp
        # lc_i=lc_i[(np.abs(lc_i['flux']/lc_i['fluxerr']))>=minsnr]

    # Store the un-lensed model as a component of the lensed SN object.
    model.set(**params)
    obj.model = model
    return(obj)


def _addML(myCurve):
    mlCurves=dict([])
    num=np.random.randint(1,5)
    print(num)
    for b in set(myCurve['band']):
        mlCurve=getDiffCurve(myCurve['time'][myCurve['band']==b],num)
        myCurve['flux'][myCurve['band']==b]/=mlCurve[b]
        myCurve['fluxerr'][myCurve['band']==b]/=mlCurve[b]
        mlCurves[b]=mlCurve[b]

    return(myCurve,mlCurves)