# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy,sncosmo
import time
from copy import deepcopy
import math
from collections import OrderedDict
import warnings

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d

from sncosmo.photdata import photometric_data, select_data
from sncosmo.utils import Result, Interp1D, ppf,integration_grid
from sncosmo.bandpasses import get_bandpass
from sncosmo.magsystems import get_magsystem

from sncosmo.constants import HC_ERG_AA, MODEL_BANDFLUX_SPACING

from .cython_helpers import *

__all__ = ['nest_lc', 'flatten_result', 'chisq']


class DataQualityError(Exception):
    pass


def generate_chisq(data, model, signature='iminuit', modelcov=False):
    """Define and return a chisq function for use in optimization.

    This function pre-computes and saves the inverse covariance matrix,
    making subsequent evaluations faster. The model covariance (if specified)
    is fixed at the time the chisq function is generated."""

    # precompute inverse covariance matrix
    cov = np.diag(data.fluxerr**2) if data.fluxcov is None else data.fluxcov
    if modelcov:
        _, mcov = model.bandfluxcov(data.band, data.time,
                                    zp=data.zp, zpsys=data.zpsys)
        cov = cov + mcov
    invcov = np.linalg.pinv(cov)

    # iminuit expects each parameter to be a separate argument (including fixed
    # parameters)
    if signature == 'iminuit':
        def chisq(*parameters):
            model.parameters = parameters
            model_flux = model.bandflux(data.band, data.time,
                                        zp=data.zp, zpsys=data.zpsys)
            diff = data.flux - model_flux
            return np.dot(np.dot(diff, invcov), diff)

    else:
        raise ValueError("unknown signature: {!r}".format(signature))

    return chisq


def chisq(data, model, modelcov=False):
    """Calculate chisq statistic for the model, given the data.

    Parameters
    ----------
    model : `~sncosmo.Model`
    data : `~astropy.table.Table` or `~numpy.ndarray` or `dict`
        Table of photometric data. Must include certain columns.
        See the "Photometric Data" section of the documentation for
        required columns.
    modelcov : bool
        Include model covariance? Calls ``model.bandfluxcov`` method
        instead of ``model.bandflux``. The source in the model must therefore
        implement covariance.

    Returns
    -------
    chisq : float
    """
    data = photometric_data(data)
    data.sort_by_time()

    if data.fluxcov is None and not modelcov:
        mflux = model.bandflux(data.band, data.time,
                               zp=data.zp, zpsys=data.zpsys)
        return np.sum(((data.flux - mflux) / data.fluxerr)**2)

    else:
        # need to invert a covariance matrix
        cov = (np.diag(data.fluxerr**2) if data.fluxcov is None
               else data.fluxcov)
        if modelcov:
            mflux, mcov = model.bandfluxcov(data.band, data.time,
                                            zp=data.zp, zpsys=data.zpsys)
            cov = cov + mcov
        else:
            mflux = model.bandflux(data.band, data.time,
                                   zp=data.zp, zpsys=data.zpsys)
        invcov = np.linalg.pinv(cov)
        diff = data.flux - mflux
        return np.dot(np.dot(diff, invcov), diff)


def flatten_result(res):
    """Turn a result from fit_lc into a simple dictionary of key, value pairs.

    Useful when saving results to a text file table, where structures
    like a covariance matrix cannot be easily written to a single
    table row.

    Parameters
    ----------
    res : Result
        Result object from `~sncosmo.fit_lc`.

    Returns
    -------
    flatres : Result
        Flattened result. Keys are all strings, values are one of: float, int,
        string), suitable for saving to a text file.
    """

    flat = Result(success=(1 if res.success else 0),
                  ncall=res.ncall,
                  chisq=res.chisq,
                  ndof=res.ndof)

    # Parameters and uncertainties
    for i, n in enumerate(res.param_names):
        flat[n] = res.parameters[i]
        if res.errors is None:
            flat[n + '_err'] = float('nan')
        else:
            flat[n + '_err'] = res.errors.get(n, 0.)

    # Covariances.
    for n1 in res.param_names:
        for n2 in res.param_names:
            key = n1 + '_' + n2 + '_cov'
            if n1 not in res.cov_names or n2 not in res.cov_names:
                flat[key] = 0.
            elif res.covariance is None:
                flat[key] = float('nan')
            else:
                i = res.cov_names.index(n1)
                j = res.cov_names.index(n2)
                flat[key] = res.covariance[i, j]

    return flat


def _mask_bands(data, model, z_bounds=None):
    if z_bounds is None:
        return model.bandoverlap(data.band)
    else:
        return np.all(model.bandoverlap(data.band, z=z_bounds), axis=1)


def _warn_dropped_bands(data, mask):
    """Warn that we are dropping some bands from the data:"""
    drop_bands = [(b.name if b.name is not None else repr(b))
                  for b in set(data.band[np.invert(mask)])]
    warnings.warn("Dropping following bands from data: " +
                  ", ".join(drop_bands) +
                  "(out of model wavelength range)", RuntimeWarning)


def cut_bands(data, model, z_bounds=None, warn=True):
    mask = _mask_bands(data, model, z_bounds=z_bounds)

    if not np.all(mask):

        # Fail if there are no overlapping bands whatsoever.
        if not np.any(mask):
            raise RuntimeError('No bands in data overlap the model.')

        if warn:
            _warn_dropped_bands(data, mask)

        data = data[mask]

    return data, mask


def t0_bounds(data, model):
    """Determine bounds on t0 parameter of the model.

    The lower bound is such that the latest model time is equal to the
    earliest data time. The upper bound is such that the earliest
    model time is equal to the latest data time."""

    return (model.get('t0') + np.min(data.time) - model.maxtime(),
            model.get('t0') + np.max(data.time) - model.mintime())


def guess_t0_and_amplitude(data, model, minsnr):
    """Guess t0 and amplitude of the model based on the data."""

    # get data above the signal-to-noise ratio cut
    significant_data = data[(data.flux / data.fluxerr) > minsnr]
    if len(significant_data) == 0:
        raise DataQualityError('No data points with S/N > {0}. Initial '
                               'guessing failed.'.format(minsnr))

    # grid on which to evaluate model light curve
    timegrid = np.linspace(model.mintime(), model.maxtime(),
                           int(model.maxtime() - model.mintime() + 1))

    # get data flux on a consistent scale in order to compare to model
    # flux light curve.
    norm_flux = significant_data.normalized_flux(zp=25., zpsys='ab')

    model_lc = {}
    data_flux = {}
    data_time = {}
    for band in set(significant_data.band):
        model_lc[band] = (
            model.bandflux(band, timegrid, zp=25., zpsys='ab') /
            model.parameters[2])
        mask = significant_data.band == band
        data_flux[band] = norm_flux[mask]
        data_time[band] = significant_data.time[mask]

    if len(model_lc) == 0:
        raise DataQualityError('No data points with S/N > {0}. Initial '
                               'guessing failed.'.format(minsnr))

    # find band with biggest ratio of maximum data flux to maximum model flux
    maxratio = float("-inf")
    maxband = None
    for band in model_lc:
        ratio = np.max(data_flux[band]) / np.max(model_lc[band])
        if ratio > maxratio:
            maxratio = ratio
            maxband = band

    # amplitude guess is the largest ratio
    amplitude = abs(maxratio)

    # time guess is time of max in the band with the biggest ratio
    data_tmax = data_time[maxband][np.argmax(data_flux[maxband])]
    model_tmax = timegrid[np.argmax(model_lc[maxband])]
    t0 = model.get('t0') + data_tmax - model_tmax

    return t0, amplitude








def nest_lc(data, model, vparam_names, bounds, guess_amplitude_bound=False,
            minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
            maxiter=None, maxcall=None, modelcov=False, rstate=None,
            verbose=False, warn=True, **kwargs):
    """Run nested sampling algorithm to estimate model parameters and evidence.

    Parameters
    ----------
    data : `~astropy.table.Table` or `~numpy.ndarray` or `dict`
        Table of photometric data. Must include certain columns.
        See the "Photometric Data" section of the documentation for
        required columns.
    model : `~sncosmo.Model`
        The model to fit.
    vparam_names : list
        Model parameters to vary in the fit.
    bounds : `dict`
        Bounded range for each parameter. Bounds must be given for
        each parameter, with the exception of ``t0``: by default, the
        minimum bound is such that the latest phase of the model lines
        up with the earliest data point and the maximum bound is such
        that the earliest phase of the model lines up with the latest
        data point.
    guess_amplitude_bound : bool, optional
        If true, bounds for the model's amplitude parameter are determined
        automatically based on the data and do not need to be included in
        `bounds`. The lower limit is set to zero and the upper limit is 10
        times the amplitude "guess" (which is based on the highest-flux
        data point in any band). Default is False.
    minsnr : float, optional
        Minimum signal-to-noise ratio of data points to use when guessing
        amplitude bound. Default is 5.
    priors : `dict`, optional
        Prior probability distribution function for each parameter. The keys
        should be parameter names and the values should be callables that
        accept a float. If a parameter is not in the dictionary, the prior
        defaults to a flat distribution between the bounds.
    ppfs : `dict`, optional
        Prior percent point function (inverse of the cumulative distribution
        function) for each parameter. If a parameter is in this dictionary,
        the ppf takes precedence over a prior pdf specified in ``priors``.
    npoints : int, optional
        Number of active samples to use. Increasing this value increases
        the accuracy (due to denser sampling) and also the time
        to solution.
    method : {'classic', 'single', 'multi'}, optional
        Method used to select new points. Choices are 'classic',
        single-ellipsoidal ('single'), multi-ellipsoidal ('multi'). Default
        is 'single'.
    maxiter : int, optional
        Maximum number of iterations. Iteration may stop earlier if
        termination condition is reached. Default is no limit.
    maxcall : int, optional
        Maximum number of likelihood evaluations. Iteration may stop earlier
        if termination condition is reached. Default is no limit.
    modelcov : bool, optional
        Include model covariance when calculating chisq. Default is False.
    rstate : `~numpy.random.RandomState`, optional
        RandomState instance. If not given, the global random state of the
        ``numpy.random`` module will be used.
    verbose : bool, optional
        Print running evidence sum on a single line.
    warn : bool, optional
        Issue warning when dropping bands outside the model range. Default is
        True.

        *New in version 1.5.0*

    Returns
    -------
    res : Result
        Attributes are:

        * ``niter``: total number of iterations
        * ``ncall``: total number of likelihood function calls
        * ``time``: time in seconds spent in iteration loop.
        * ``logz``: natural log of the Bayesian evidence Z.
        * ``logzerr``: estimate of uncertainty in logz (due to finite sampling)
        * ``h``: Bayesian information.
        * ``vparam_names``: list of parameter names varied.
        * ``samples``: 2-d `~numpy.ndarray`, shape is (nsamples, nparameters).
          Each row is the parameter values for a single sample. For example,
          ``samples[0, :]`` is the parameter values for the first sample.
        * ``logprior``: 1-d `~numpy.ndarray` (length=nsamples);
          log(prior volume) for each sample.
        * ``logl``: 1-d `~numpy.ndarray` (length=nsamples); log(likelihood)
          for each sample.
        * ``weights``: 1-d `~numpy.ndarray` (length=nsamples);
          Weight corresponding to each sample. The weight is proportional to
          the prior * likelihood for the sample.
        * ``parameters``: 1-d `~numpy.ndarray` of weighted-mean parameter
          values from samples (including fixed parameters). Order corresponds
          to ``model.param_names``.
        * ``covariance``: 2-d `~numpy.ndarray` of parameter covariance;
          indicies correspond to order of ``vparam_names``. Calculated from
          ``samples`` and ``weights``.
        * ``errors``: OrderedDict of varied parameter uncertainties.
          Corresponds to square root of diagonal entries in covariance matrix.
        * ``ndof``: Number of degrees of freedom (len(data) -
          len(vparam_names)).
        * ``bounds``: Dictionary of bounds on varied parameters (including
          any automatically determined bounds).
        * ``data_mask``: Boolean array the same length as data specifying
          whether each observation was used.
          *New in version 1.5.0.*

    estimated_model : `~sncosmo.Model`
        A copy of the model with parameters set to the values in
        ``res.parameters``.
    """

    try:
        import nestle
    except ImportError:
        raise ImportError("nest_lc() requires the nestle package.")

    # experimental parameters
    tied = kwargs.get("tied", None)

    data = photometric_data(data)

    # sort by time
    if not np.all(np.ediff1d(data.time) >= 0.0):
        sortidx = np.argsort(data.time)
        data = data[sortidx]
    else:
        sortidx = None

    model = copy.copy(model)
    bounds = copy.copy(bounds)  # need to copy this b/c we modify it below

    # Order vparam_names the same way it is ordered in the model:
    vparam_names = [s for s in model.param_names if s in vparam_names]

    # Drop data that the model doesn't cover.
    fitdata, data_mask = cut_bands(data, model,
                                   z_bounds=bounds.get('z', None),
                                   warn=warn)

    if guess_amplitude_bound:
        if model.param_names[2] not in vparam_names:
            raise ValueError("Amplitude bounds guessing enabled but "
                             "amplitude parameter {0!r} is not varied"
                             .format(model.param_names[2]))
        if model.param_names[2] in bounds:
            raise ValueError("cannot supply bounds for parameter {0!r}"
                             " when guess_amplitude_bound=True"
                             .format(model.param_names[2]))

        # If redshift is bounded, set model redshift to midpoint of bounds
        # when doing the guess.
        if 'z' in bounds:
            model.set(z=sum(bounds['z']) / 2.)
        _, amplitude = guess_t0_and_amplitude(fitdata, model, minsnr)
        bounds[model.param_names[2]] = (0., 10. * amplitude)

    # Find t0 bounds to use, if not explicitly given
    if 't0' in vparam_names and 't0' not in bounds:
        bounds['t0'] = t0_bounds(fitdata, model)

    if ppfs is None:
        ppfs = {}
    if tied is None:
        tied = {}

    # Convert bounds/priors combinations into ppfs
    if bounds is not None:
        for key, val in bounds.items():
            if key in ppfs:
                continue  # ppfs take priority over bounds/priors
            a, b = val
            if priors is not None and key in priors:
                # solve ppf at discrete points and return interpolating
                # function
                x_samples = np.linspace(0., 1., 101)
                ppf_samples = ppf(priors[key], x_samples, a, b)
                f = Interp1D(0., 1., ppf_samples)
            else:
                f = Interp1D(0., 1., np.array([a, b]))
            ppfs[key] = f

    # NOTE: It is important that iparam_names is in the same order
    # every time, otherwise results will not be reproducible, even
    # with same random seed.  This is because iparam_names[i] is
    # matched to u[i] below and u will be in a reproducible order,
    # so iparam_names must also be.
    iparam_names = [key for key in vparam_names if key in ppfs]
    ppflist = [ppfs[key] for key in iparam_names]
    npdim = len(iparam_names)  # length of u
    ndim = len(vparam_names)  # length of v

    # Check that all param_names either have a direct prior or are tied.
    for name in vparam_names:
        if name in iparam_names:
            continue
        if name in tied:
            continue
        raise ValueError("Must supply ppf or bounds or tied for parameter '{}'"
                         .format(name))

    def prior_transform(u):
        d = {}
        for i in range(npdim):
            d[iparam_names[i]] = ppflist[i](u[i])
        v = np.empty(ndim, dtype=np.float)
        for i in range(ndim):
            key = vparam_names[i]
            if key in d:
                v[i] = d[key]
            else:
                v[i] = tied[key](d)
        return v

    # Indicies of the model parameters in vparam_names
    idx = np.array([model.param_names.index(name) for name in vparam_names])


    fitdata = photometric_data(fitdata)
    fitdata.sort_by_time()



    fluxdata=np.array(fitdata.flux).astype(np.float32)
    fluxerrdata=np.array(fitdata.fluxerr).astype(np.float32)
    cov = (np.diag(fitdata.fluxerr**2) if fitdata.fluxcov is None else fitdata.fluxcov)

    all_wave=[]
    all_dwave=[]
    all_trans=[]
    all_bands=np.zeros((len(fitdata),2),dtype=np.int32)
    all_zps=np.array(fitdata.zp).astype(np.float32)
    all_zpsys=np.array([get_magsystem(fitdata.zpsys[i]).zpbandflux(fitdata.band[i]) for i in range(len(fitdata))])


    for b in set(fitdata.band):
        wave,dwave,trans=_bandflux_single_spacing(b)
        inds=np.where(np.array([x.name for x in fitdata.band])==b.name)[0]
        all_bands[inds,0]=len(all_wave)

        all_bands[inds,1]=len(all_wave)+len(wave)

        all_wave=np.append(all_wave,wave)
        all_dwave=np.append(all_dwave,[dwave]*len(wave))
        all_trans=np.append(all_trans,trans)




    bandDict={'all_wave':all_wave.astype(np.float32),
              'all_dwave':all_dwave.astype(np.float32),
              'all_trans':all_trans.astype(np.float32),
              'all_zps':all_zps,
              'all_zpsys':all_zpsys.astype(np.float32),
              'all_bands':all_bands.astype(np.int32)}

    bandDict2={b:_bandflux_single_spacing(b) for b in set(fitdata.band)}

    def loglike(parameters):
        model.parameters[idx] = parameters
        if fitdata.fluxcov is None and not modelcov:
            mflux = sntd_bandflux(model,fitdata.band, fitdata.time,bandDict,bandDict2,
                                   zp=fitdata.zp, zpsys=fitdata.zpsys).astype(np.float32)
            return -0.5*fast_chisq(fluxdata,fluxerrdata,mflux)

        else:
            if modelcov:
                mflux, mcov = model.bandfluxcov(fitdata.band, fitdata.time,
                                                zp=fitdata.zp, zpsys=fitdata.zpsys)
                cov += mcov
            else:
                mflux = model.bandflux(data.band, data.time,zp=data.zp, zpsys=data.zpsys).astype(np.float32)
            invcov = np.linalg.pinv(cov)
            diff = data.flux - mflux
            return -0.5*np.dot(np.dot(diff, invcov), diff)



    t0 = time.time()
    res = nestle.sample(loglike, prior_transform, ndim, npdim=npdim,
                        npoints=npoints, method=method, maxiter=maxiter,
                        maxcall=maxcall, rstate=rstate,
                        callback=(nestle.print_progress if verbose else None))
    elapsed = time.time() - t0

    # estimate parameters and covariance from samples
    vparameters, cov = nestle.mean_and_cov(res.samples, res.weights)

    # update model parameters to estimated ones.
    model.set(**dict(zip(vparam_names, vparameters)))

    # If we need to, unsort the mask so mask applies to input data
    if sortidx is not None:
        unsort_idx = np.argsort(sortidx)  # indicies that will unsort array
        data_mask = data_mask[unsort_idx]

    # `res` is a nestle.Result object. Collect result into a sncosmo.Result
    # object for consistency, and add more fields.
    res = Result(niter=res.niter,
                 ncall=res.ncall,
                 logz=res.logz,
                 logzerr=res.logzerr,
                 h=res.h,
                 samples=res.samples,
                 weights=res.weights,
                 logvol=res.logvol,
                 logl=res.logl,
                 vparam_names=copy.copy(vparam_names),
                 ndof=len(fitdata) - len(vparam_names),
                 bounds=bounds,
                 time=elapsed,
                 parameters=model.parameters.copy(),
                 covariance=cov,
                 errors=OrderedDict(zip(vparam_names,
                                        np.sqrt(np.diagonal(cov)))),
                 param_dict=OrderedDict(zip(model.param_names,
                                            model.parameters)),
                 data_mask=data_mask)

    return res, model



def sntd_bandflux(self, band, time, bandDict,bandDict2,zp=None, zpsys=None):
    try:
        return _sntd_bandflux(self, band, time, zp, zpsys,bandDict,bandDict2)
    except ValueError as e:
        sncosmo.models._check_for_fitpack_error(e, time, 'time')
        raise e

def change_sncosmo_model(model):
    sntd_model=deepcopy(model)
    if isinstance(model._source,sncosmo.SALT2Source):
        sntd_model._source._flux=salt2_source_flux_fast



    elif isinstance(model._source,sncosmo.TimeSeriesSource):
        pass
        #sntd_model._source._flux=timeseries_source_flux_fast

    elif isinstance(model._source,sncosmo.StretchSource):
        #sntd_model._source._flux=stretch_source_flux_fast
        pass

    elif isinstance(model._source,sncosmo.MLCS2k2Source):
        #sntd_model._source._flux=mlcs_source_flux_fast
        pass

    elif isinstance(model._source,sncosmo.SNEMOSource):
        #sntd_model._source._flux=snemo_source_flux_fast
        pass

    else:
        print("Don't recognize model source, can't speed up fitting.")
        return model


def _bandflux_single_spacing(band):

    # Set up wavelength grid. Spacing (dwave) evenly divides the bandpass,
    # closest to 5 angstroms without going over.
    wave, dwave = integration_grid(band.minwave(), band.maxwave(),
                                   MODEL_BANDFLUX_SPACING)
    trans = band(wave)
    return(wave,dwave,trans)

import time as timeit

def sntd_model_flux(self,time, wave):
    """Replacement for sncosmo Array flux function."""
    a = 1. / (1. + self._parameters[0])
    phase = (time - self._parameters[1]) * a
    minphase = (self.mintime() - self._parameters[1]) * a
    maxphase = (self.maxtime() - self._parameters[1]) * a
    restwave = wave * a
    # Note that below we multiply by the scale factor to conserve
    # bolometric luminosity.
    f = a * self._source._flux(phase, restwave)

    # Pass the flux through the PropagationEffects.
    for effect, frame, zindex in zip(self._effects, self._effect_frames,
                                     self._effect_zindicies):
        if frame == 'obs':
            effect_wave = wave
            effect_phase=phase*(1./a)
        elif frame == 'rest':
            effect_wave = restwave
            effect_phase=phase
        else:  # frame == 'free'
            effect_a = 1. / (1. + self._parameters[zindex])
            effect_wave = wave * effect_a
            effect_phase=phase/a*(1.+self._parameters[zindex])


        f = effect.propagate(effect_phase,effect_wave, f)

    # f = salt2_source_flux_fast(self._source._model['M0'](phase, restwave).astype(np.float32),
    #                            self._source._model['M1'](phase, restwave).astype(np.float32),
    #                            self._source._colorlaw(restwave).astype(np.float32),
    #                            self._source._parameters[0],
    #                            self._source._parameters[1],
    #                            self._source._parameters[2])/(1+self._parameters[0])

    return f



def _sntd_bandflux(model, band, time_or_phase, zp, zpsys,bandDict,bandDict2):
    """Support function for bandflux in Source and Model.
    This is necessary to have outside because ``phase`` is used in Source
    and ``time`` is used in Model, and we want the method signatures to
    have the right variable name.
    """

    # initialize output arrays
    #bandflux = np.zeros(time_or_phase.shape[0], dtype=np.float32)
    ndim = time_or_phase.ndim  # Save input ndim for return val.
    # Loop over unique bands.
    import sys
    # t1=timeit.time()
    # for b in bandDict2.keys():
    #     mask = band == b
    #     wave,dwave,trans=bandDict2[b]
    #     t=timeit.time()
    #     flux=model._flux(time_or_phase,wave)
    #     print('model',timeit.time()-t)
    #     t=timeit.time()
    #     fsum = fast_bandflux_single(flux[mask].astype(np.float32),wave.astype(np.float32),trans.astype(np.float32),dwave,HC_ERG_AA)
    #     print('single',timeit.time()-t)
    #     t=timeit.time()
    #     if zp is not None:
    #         zpnorm = 10.**(0.4 * zp[mask])
    #         bandzpsys = zpsys[mask]
    #         for ms in set(bandzpsys):
    #             mask2 = bandzpsys == ms
    #             ms = get_magsystem(ms)
    #             zpnorm[mask2] = zpnorm[mask2] / ms.zpbandflux(b)
    #         fsum *= zpnorm
    #     print('zp',timeit.time()-t)
    #     bandflux[mask] = fsum
    # print(timeit.time()-t1)

    #c test
    #print('c test')

    t=timeit.time()

    flux=model._flux(time_or_phase.astype(np.double),bandDict['all_wave'].astype(np.double))

    #print('model',timeit.time()-t)
    #t=timeit.time()

    bandflux = fast_bandflux_single(flux.astype(np.float32),
                                bandDict['all_wave'].astype(np.float32),
                                bandDict['all_trans'],
                                bandDict['all_dwave'],bandDict['all_zps'],
                                bandDict['all_zpsys'],bandDict['all_bands'],
                                HC_ERG_AA)
    print('single',timeit.time()-t)
    sys.exit()

    if ndim == 0:
        return bandflux[0]
    return bandflux