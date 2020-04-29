import sncosmo,os
from copy import deepcopy,copy
from collections import OrderedDict

from astropy.io import ascii
import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
from sncosmo.utils import alias_map

from .util import __filedir__
from .curve_io import curve,curveDict
from .ml import *

__all__=['createMultiplyImagedSN']

OBSERVATIONS_REQUIRED_ALIASES = ('time', 'band', 'zp', 'zpsys', 'gain',
                                 'skynoise')

OBSERVATIONS_ALIASES = OrderedDict([
    ('time', set(['time', 'date', 'jd', 'mjd', 'mjdobs', 'mjd_obs'])),
    ('band', set(['band', 'bandpass', 'filter', 'flt'])),
    ('zp', set(['zp', 'zpt', 'zeropoint', 'zero_point'])),
    ('zpsys', set(['zpsys', 'zpmagsys', 'magsys'])),
    ('gain', set(['gain'])),
    ('skynoise', set(['skynoise']))
])

def _getAbsoluteDist():
    absolutes=ascii.read(os.path.join(__filedir__,'sim','data','absolutes.ref'))
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


def createMultiplyImagedSN(
        sourcename, snType, redshift,z_lens=None, telescopename='telescope',
        objectName='object', time_delays=[10., 50.], magnifications=[2., 1.],sn_params={},av_dists={},
        numImages=2, cadence=5, epochs=30, clip_time=[-30,150],bands=['F105W', 'F160W'],start_time=None,
        gain=200., skynoiseRange=(1, 1.1), timeArr=None,zpsys='ab', zp=None,hostr_v=3.1,lensr_v=3.1,
        microlensing_type=None, microlensing_params=[],ml_loc=[None,None],
        dust_model='CCM89Dust', av_host=.3, av_lens=0,fix_luminosity=False,
        minsnr=0.0, scatter=True,snrFunc=None):
    """
    Generate a multiply-imaged SN light curve set, with user-specified time
    delays and magnifications.

    Parameters
    ----------
    sourcename: :class:`~sncosmo.Source` or str
        The model for the spectral evolution of the source. If a string
        is given, it is used to retrieve a :class:`~sncosmo.Source` from
        the registry.
    snType : str
        The classification of the supernova
    redshift : float
        Redshift of the source
    z_lens : float
        Redshift of the lens
    telescopename : str
        The name of the telescope used for observations
    objectName : str
        The name of the simulated supernova
    time_delays : :class:`~list` of :class:`~float`
        The relative time delays for the multiple images of the supernova. Must
        be same length as numImages
    magnifications : :class:`~list` of :class:`~float`
        The relative magnifications for the multiple images of hte supernova. Must
        be same length as numImages
    sn_params: :class:`~dict`
        A dictionary with SN model parameters as keys, and some sort of pdf or
        other callable function that returns the model parameter.
    av_dists: :class:`~dict`
        A dictionary with 'host' or 'lens' as keys, and some sort of pdf or
        other callable function that returns the relevant av parameter.
    numImages : int
        The number of images to simulate
    cadence : float
        The cadence of the simulated observations (if timeArr is not defined)
    epochs : int
        The number of simulated observations (if timeArr is not defined)
    clip_time: :class:`~list`
        Rest frame phase start and end time, will clip output table to these values.
    bands : :class:`~list` of :class:`~sncosmo.Bandpass` or :class:`~str`
        The bandpass(es) used for simulated observations
    start_time: float
        Start time for the leading image. If None, start will be the first value in 
        the time array, and the peak will be this ± the 
        relative time delay for the leading image. 
    gain : float
        Gain of the telescope "obtaining" the simulated observations (if snrFunc
        not defined)
    skynoiseRange : :class:`~list` of :class:`~float`
        The left and right bounds of sky noise used to define observational noise
        (if snrFunc not defined)
    timeArr : :class:`~list` of :class:`~float`
        A list of times that define the simulated observation epochs
    zpsys : str or :class:`~sncosmo.MagSystem`
        The zero-point system used to define the photometry
    zp : float or :class:`~list` of :class:`~float`
        The zero-point used to define the photometry, list if simulating multiple
        bandpasses. Then this list must be the same length as bands
    hostr_v: float
        The r_v parameter for the host.
    lensr_v: float
        The r_v parameter for the lens.
    microlensing_type : str
        If microlensing is to be included, defines whether it is
        "AchromaticSplineMicrolensing" or "AchromaticMicrolensing"
    microlensing_params: :class:`~numpy.array` or :class:`~list` of :class:`~int`
        If using AchromaticSplineMicrolensing, then this params list must give
        three values for [nanchor, sigmadm, nspl]. If using AchromaticMicrolensing,
        then this must be a microcaustic defined by a 2D numpy array
    ml_loc: :class:`~list`
        List containing tuple locations of SN on microlensing map (random if Nones)
    dust_model : str
        The dust model to be used for simulations, see sncosmo documentation for options
    av_host : float
        The A<sub>V</sub> parameter for the simulated dust effect in the source plane
    av_lens : float
        The A<sub>V</sub> parameter for the simulated dust effect in the lens plane
    fix_luminosity: bool
        Set the luminosity of every SN to be the peak of the distribution
    minsnr : float
        A minimum SNR threshold for observations when defining uncertainty
    scatter : bool
        Boolean that decides whether Gaussian scatter is applied to simulated
        observations
    snrFunc : :class:`~scipy.interpolate.interp1d` or :class:`~dict`
        An interpolation function that defines the signal to noise ratio (SNR)
        as a function of magnitude in the AB system. Used to define the
        observations instead of telescope parameters like gain and skynoise.
        This can be a dictionary, so that it's different for each filter with
        filters as keys and interpolation functions as values.

    Returns
    -------
    MISN: :class:`~sntd.curveDict`
        A curveDict object containing each of the multiply-imaged SN light curves
        and the simulation parameters.
    Examples
    --------
    >>> myMISN = sntd.createMultiplyImagedSN('salt2', 'Ia', 1.33,z_lens=.53, bands=['F110W'],
        zp=[26.8], cadence=5., epochs=35.,skynoiseRange=(.001,.005),gain=70. , time_delays=[10., 78.],
        magnifications=[7,3.5], objectName='My Type Ia SN', telescopename='HST',minsnr=5.0)
    """

    if timeArr is not None:
        times=timeArr
    else:
        times=np.linspace(0,int(cadence*epochs),int(epochs))
        if start_time is not None:
            times+=start_time
    leading_peak=times[0]
    bandList=np.array([np.tile(b,len(times)) for b in bands]).flatten()
    ms=sncosmo.get_magsystem(zpsys)

    if zp is None:
        zpList=[ms.band_flux_to_mag(1,b) for b in bandList]
    elif isinstance(zp,(list,tuple)):
        zpList=np.array([np.tile(z,len(times)) for z in zp]).flatten()
    else:
        zpList=[zp for i in range(len(bandList))]

    #set up object to be filled by simulations
    curve_obj=curveDict(telescopename=telescopename,object=objectName)
    curve_obj.bands = set(bandList)

    #make sncosmo obs table
    obstable = Table({'time':np.tile(times,len(bands)), 'band':bandList,
                      'zpsys':[zpsys.upper() for i in range(len(bandList))],
                      'zp':zpList,
                      'skynoise':np.random.uniform(
                          skynoiseRange[0],skynoiseRange[1],len(bandList)),
                      'gain':[gain for i in range(len(bandList))]})
    absolutes=_getAbsoluteDist()

    # Set up the dust_model extinction effects in the host galaxy and lens plane
    # TODO allow additional dust screens, not in the host or lens plane?
    # TODO sample from a prior for host and lens-plane dust A_V?
    # TODO : allow different lens-plane dust_model for each image?
    RV_lens = lensr_v
    RV_host = hostr_v
    dust_frames = []
    dust_names = []
    dust_effect_list = []
    if dust_model and (av_lens or av_host or len(av_dists)>0):
        dust_effect = {'CCM89Dust': sncosmo.CCM89Dust,
                       'OD94Dust': sncosmo.OD94Dust,
                       'F99Dust': sncosmo.F99Dust}[dust_model]()
        if av_host or 'host' in av_dists.keys():
            dust_frames.append('rest')
            dust_names.append('host')
            dust_effect_list.append(dust_effect)
        if av_lens or 'lens' in av_dists.keys():
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
    #set absolute magnitude in b or r band based on literature
    if snType in ['IIP','IIL','IIn']:
        absBand='bessellb'
    else:
        absBand='bessellr'
    if model.param_names[2] not in sn_params.keys():
        if fix_luminosity:
            model.set_source_peakabsmag(absolutes[snType]['dist'][0],
                                        absBand, zpsys)
        else:
            model.set_source_peakabsmag(_getAbsFromDist(absolutes[snType]['dist']),
                                    absBand, zpsys)
    else:
        model.parameters[2]=sn_params[model.param_names[2]]()
    

    t0=leading_peak
    if snType=='Ia':
        x0=model.get('x0')
        params={'z':redshift, 't0':t0, 'x0':x0}
        if 'x1' not in sn_params.keys():
            params['x1']=np.random.normal(0.,1.)
        else:
            params['x1']=sn_params['x1']() 
        if 'c' not in sn_params.keys():
            params['c']=np.random.normal(0.,.1)
        else:
            params['c']=sn_params['c']()
    else:
        amp=model.get('amplitude')
        params={'z':redshift, 't0':t0, 'amplitude':amp}
    model.set(**params)
    if av_host or 'host' in av_dists.keys():
        if 'host' in av_dists.keys():
            av_host=av_dists['host']()
        ebv_host = av_host/RV_host
        model.set(hostebv=ebv_host, hostr_v=RV_host)
    else:
        ebv_host = 0
    if av_lens or 'lens' in av_dists.keys():
        if z_lens is None:
            print('No z_lens set, assuming half of source z...')
            z_lens = redshift / 2.
        if 'lens' in av_dists.keys():
            av_lens=av_dists['lens']()
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
        model_i._flux=_mlFlux
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
                    ml_spline_func = ChromaticSplineMicrolensing
                ml_effect = ml_spline_func(nanchor=nanchor, sigmadm=sigmadm,
                                           nspl=nspl)
            else:
                #get magnification curve from the defined microcaustic
                mlTime=np.arange(0,times[-1]/(1+redshift)-model_i._source._phase[0]+5,1)

                time,dmag=microcaustic_field_to_curve(microlensing_params,mlTime,z_lens,redshift,loc=ml_loc[imnum])
                dmag/=np.mean(dmag) #to remove overall magnification

                ml_effect = AchromaticMicrolensing(
                    time+model_i._source._phase[0],dmag, magformat='multiply')
                
            model_i.add_effect(ml_effect, 'microlensing', 'rest')
        else:
            ml_effect = None

        # Generate the simulated SN light curve observations, make a `curve`
        # object, and store the simulation metadata
        model_i.set(**params_i)

        table_i = realize_lcs(
            obstable , model_i, [params_i],
            trim_observations=True, scatter=scatter,thresh=minsnr,snrFunc=snrFunc)
        tried=0
        while (len(table_i)==0 or len(table_i[0])<numImages) and tried<50:
            table_i = realize_lcs(
                obstable , model_i, [params_i],
                trim_observations=True, scatter=scatter,thresh=minsnr,snrFunc=snrFunc)
            tried+=1
        if tried==50:
            #this arbitrary catch is here in case your minsnr and observation parameters
            #result in a "non-detection"
            print("Your survey parameters detected no supernovae.")
            return None
        table_i=table_i[0]
        if timeArr is None:
            table_i=table_i[table_i['time']<model_i.get('t0')+clip_time[1]*(1+model_i.get('z'))]
            table_i=table_i[table_i['time']>model_i.get('t0')+clip_time[0]*(1+model_i.get('z'))]
        #create is curve with all parameters and add it to the overall curveDict object from above
        curve_i=curve()
        curve_i.object=None
        curve_i.zpsys=zpsys
        curve_i.table=deepcopy(table_i)
        curve_i.bands=list(set(table_i['band']))
        curve_i.simMeta=deepcopy(table_i.meta)
        curve_i.simMeta['sourcez']=redshift
        curve_i.simMeta['model']=copy(model_i)
        curve_i.simMeta['hostebv']=ebv_host
        curve_i.simMeta['lensebv']=ebv_lens
        curve_i.simMeta['lensz']=z_lens
        curve_i.simMeta['mu']=mu
        curve_i.simMeta['td']=td
        curve_i.simMeta['microlensing'] = ml_effect
        curve_i.simMeta['microlensing_type'] = microlensing_type

        if microlensing_type=='AchromaticSplineMicrolensing':
            curve_i.simMeta['microlensing_params'] = microlensing_params
        elif microlensing_type is not None:
            curve_i.simMeta['microlensing_params'] = interp1d(time+model_i._source._phase[0],dmag)

        curve_obj.add_curve(curve_i)


    # Store the un-lensed model as a component of the lensed SN object.
    model.set(**params)
    curve_obj.model = model

    return(curve_obj)

def realize_lcs(observations, model, params, thresh=None,
                trim_observations=False, scatter=True,snrFunc=None):
    """***A copy of SNCosmo's function, just to add a SNR function
    Realize data for a set of SNe given a set of observations.

    Parameters
    ----------
    observations : :class:`~astropy.table.Table` or :class:`~numpy.ndarray`
        Table of observations. Must contain the following column names:
        ``band``, ``time``, ``zp``, ``zpsys``, ``gain``, ``skynoise``.
    model : `sncosmo.Model`
        The model to use in the simulation.
    params : :class:`~list` (or generator) of dict
        List of parameters to feed to the model for realizing each light curve.
    thresh : float, optional
        If given, light curves are skipped (not returned) if none of the data
        points have signal-to-noise greater than ``thresh``.
    trim_observations : bool, optional
        If True, only observations with times between
        ``model.mintime()`` and ``model.maxtime()`` are included in
        result table for each SN. Default is False.
    scatter : bool, optional
        If True, the ``flux`` value of the realized data is calculated by
        adding  a random number drawn from a Normal Distribution with a
        standard deviation equal to the ``fluxerror`` of the observation to
        the bandflux value of the observation calculated from model. Default
        is True.

    Returns
    -------
    sne : list of :class:`~astropy.table.Table`
        Table of realized data for each item in ``params``.

    Notes
    -----
    ``skynoise`` is the image background contribution to the flux measurement
    error (in units corresponding to the specified zeropoint and zeropoint
    system). To get the error on a given measurement, ``skynoise`` is added
    in quadrature to the photon noise from the source.

    It is left up to the user to calculate ``skynoise`` as they see fit as the
    details depend on how photometry is done and possibly how the PSF is
    is modeled. As a simple example, assuming a Gaussian PSF, and perfect
    PSF photometry, ``skynoise`` would be ``4 * pi * sigma_PSF * sigma_pixel``
    where ``sigma_PSF`` is the standard deviation of the PSF in pixels and
    ``sigma_pixel`` is the background noise in a single pixel in counts.

    """

    RESULT_COLNAMES = ('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys')
    lcs = []

    # Copy model so we don't mess up the user's model.
    model = copy(model)


    # get observations as a Table
    if not isinstance(observations, Table):
        if isinstance(observations, np.ndarray):
            observations = Table(observations)
        else:
            raise ValueError("observations not understood")

    # map column name aliases
    colname = alias_map(observations.colnames, OBSERVATIONS_ALIASES,
                        required=OBSERVATIONS_REQUIRED_ALIASES)

    # result dtype used when there are no observations
    band_dtype = observations[colname['band']].dtype
    zpsys_dtype = observations[colname['zpsys']].dtype
    result_dtype = ('f8', band_dtype, 'f8', 'f8', 'f8', zpsys_dtype)

    for p in params:
        model.set(**p)
        # Select times for output that fall within tmin amd tmax of the model
        if trim_observations:
            mask = ((observations[colname['time']] > model.mintime()) &
                    (observations[colname['time']] < model.maxtime()))
            snobs = observations[mask]
        else:
            snobs = observations

        # explicitly detect no observations and add an empty table
        if len(snobs) == 0:
            if thresh is None:
                lcs.append(Table(names=RESULT_COLNAMES,
                                 dtype=result_dtype, meta=p))
            continue
        flux = model.bandflux(snobs[colname['band']],
                              snobs[colname['time']],
                              zp=snobs[colname['zp']],
                              zpsys=snobs[colname['zpsys']])
        if snrFunc is not None:
            if isinstance(snrFunc,dict):
                fluxerr=np.ones(len(flux))
                for b in snobs[colname['band']]:
                    inds=np.where(snobs[colname['band']]==b)[0]
                    fluxerr[inds]=np.abs(flux[inds]/snrFunc[b](-2.5*np.log10(flux[inds])+snobs[colname['zp']][inds]))
            else:
                fluxerr=np.abs(flux/snrFunc(-2.5*np.log10(flux)+snobs[colname['zp']]))
        else:
            fluxerr = np.sqrt(snobs[colname['skynoise']]**2 +
                              np.abs(flux) / snobs[colname['gain']])

        # Scatter fluxes by the fluxerr
        # np.atleast_1d is necessary here because of an apparent bug in
        # np.random.normal: when the inputs are both length 1 arrays,
        # the output is a Python float!
        if scatter:
            flux = np.atleast_1d(np.random.normal(flux, fluxerr))

        # Check if any of the fluxes are significant
        if thresh is not None:
            if not np.any(flux/fluxerr > thresh):
                continue
            else:
                inds=np.where(flux/fluxerr>thresh)[0]

        data = [snobs[colname['time']][inds], snobs[colname['band']][inds], flux[inds], fluxerr[inds],
                snobs[colname['zp']][inds], snobs[colname['zpsys']][inds]]

        lcs.append(Table(data, names=RESULT_COLNAMES, meta=p))

    return lcs