import inspect,sncosmo,os,sys,warnings,pyParz,math
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy,copy
from scipy import stats
from astropy.table import Table
from astropy.extern import six
import nestle
from collections import OrderedDict
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
import scipy

from .util import *
from .curve_io import _sntd_deepcopy
from .models import *
from .ml import *

__all__=['fit_data']

__thetaSN__=['z','hostebv','screenebv','screenz','rise','fall','sigma','k','x1','c']
__thetaL__=['t0','amplitude','dt0','A','B','t1','psi','phi','s','x0']


_needs_bounds={'z'}

class newDict(dict):
    """
    This is just a dictionary replacement class that allows the use of a normal dictionary but with the ability
    to access via "dot" notation.
    """
    def __init__(self):
        super(newDict,self).__init__()

    # these three functions allow you to access the curveDict via "dot" notation
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    __getattr__ = dict.__getitem__

    def __getstate__(self):
        """
        A function necessary for pickling
        :return: self
        """
        return self

    def __setstate__(self, d):
        """
        A function necessary for pickling
        :param d: A value
        :return: self.__dict__
        """
        self.__dict__ = d



def fit_data(curves, snType='Ia',bands=None, models=None, params=None, bounds={}, ignore=None, constants=None,
             method='separate',t0_guess=None,refModel=None,effect_names=[],effect_frames=[],fitting_method='nest',
             dust=None,flip=False,guess_amplitude=True,combinedError=None,showPlots=False,microlensing=None,
             kernel='RBF',combinedGrids=None,refImage='image_1',nMicroSamples=100,color_curve=None,**kwargs):

    """
    The main high-level fitting function.

    Parameters
    ----------
    curves: :class:`~sntd.curve_io.curveDict`
        The curveDict object containing the multiple images to fit.
    snType: str
        The supernova classification
    bands: list of :class:`sncosmo.Bandpass` or str, or :class:`sncosmo.Bandpass` or str
        The band(s) to be fit
    models: list of :class:`sncosmo.Model` or str, or :class:`sncosmo.Model` or str
        The model(s) to be used for fitting to the data
    params: list of str
        The parameters to be fit for the models inside of the parameter models
    bounds: :class:`dict`
        A dictionary with parameters in params as keys and a tuple of bounds as values
    ignore: list of str
        List of parameters to ignore
    constants: :class:`dict`
        Dictionary with parameters as keys and the constant value you want to set them to as values
    method: str
        Needs to be 'separate', 'combined', or 'color'
    t0_guess: :class:`dict`
        Dictionary with image names (i.e. 'image_1','image_2') as keys and a guess for time of peak as values
    refModel: :class:`scipy.interpolate.interp1d` or :class:`sncosmo.Model`
        If doing a combined or color fit, a reference model that will be used to fit all the data at once is required.
    effect_names: list of str
        List of effect names if model contains a :class:`sncosmo.PropagationEffect`.
    effect_frames: list of str
        List of the frames (e.g. obs or rest) that correspond to the effects in effect_names
    fitting_method: str
        Must be 'nest', 'mcmc', or 'minuit'. This is only used when you are trying to find the best of a list of models.
    dust: :class:`sncosmo.PropagationEffect`
        An sncosmo dust propagation effect to include in the model
    flip: Boolean
        If True, the time axis of the model is flipped
    guess_amplitude: Boolean
        If True, the amplitude parameter for the model is estimated, as well as its bounds
    combinedError: :class:`scipy.interpolate.interp1d`
        If doing a combined or color fit, this is the uncertainty on the reference model.
    showPlots: Boolean
        If true, :func:`sncosmo.plot_lc` function is called during the fitting
    microlensing: str
        If None microlensing is ignored, otherwise should be str (e.g. achromatic, chromatic)
    kernel: str
        The kernel to use for microlensing GPR
    combinedGrids: :class:`dict`
        A dictionary with 'td' or 'mu' as keys and tuples with additive bounds as values
    refImage: str
        The name of the image you want to be the reference image (i.e. image_1,image_2, etc.)
    nMicroSamples: int
        The number of pulls from the GPR posterior you want to use for microlensing uncertainty estimation
    color_curve: :class:`astropy.Table`
        A color curve to define the relationship between bands for parameterized light curve model.
    Returns
    -------
    fitted_curveDict: :class:`~sntd.curve_io.curveDict`
        The same curveDict that was passed to fit_data, but with new fits and time delay measurements included
    Examples
    --------
    >>> fitCurves=sntd.fit_data(myMISN2,snType='Ia', models='salt2-extended',bands=['F110W','F125W'],
        params=['x0','x1','t0','c'],constants={'z':1.33},bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},
        method='separate',microlensing=None)

    """


    args = locals()

    args['curves']=_sntd_deepcopy(curves)
    args['bands'] = [bands] if bands and not isinstance(bands,(tuple,list)) else bands
    #sets the bands to user's if defined (set, so that they're unique), otherwise to all the bands that exist in curves
    args['bands'] = list(set(bands)) if bands else list(curves.bands)

    models=[models] if models and not isinstance(models,(tuple,list)) else models
    if not models:
        mod,types=np.loadtxt(os.path.join(__dir__,'data','sncosmo','models.ref'),dtype='str',unpack=True)
        modDict={mod[i]:types[i] for i in range(len(mod))}
        if snType!='Ia':
            mods = [x[0] for x in sncosmo.models._SOURCES._loaders.keys() if x[0] in modDict.keys() and modDict[x[0]][:len(snType)]==snType]
        elif snType=='Ia':
            mods = [x[0] for x in sncosmo.models._SOURCES._loaders.keys() if 'salt2' in x[0]]
    else:
        mods=models
    mods=set(mods)

    args['sn_func'] = {'minuit': sncosmo.fit_lc, 'mcmc': sncosmo.mcmc_lc, 'nest': sncosmo.nest_lc}
    #get any properties set in kwargs that exist for the defined fitting function
    args['props'] = {x: kwargs[x] for x in kwargs.keys() if
                     x in [y for y in inspect.getargspec(args['sn_func'][fitting_method])[0]] and x != 'verbose'}

    if method not in ['separate','combined','color']:
        raise RuntimeError('Parameter "method" must be "separate","combined", or "color".')
    if method=='separate':
        curves=_fitSeparate(args['curves'],mods,args,bounds,**kwargs)
    elif method=='combined':

        curves=_fitCombined(args['curves'],mods,args,bounds,combinedGrids,**kwargs)
    else:
        curves=_fitColor(args['curves'],args,bounds,combinedGrids,**kwargs)
    return curves

def _fitColor(curves,args,bounds,grids,guess_amplitude_bound=False,
              minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
              maxiter=None, maxcall=None, modelcov=False, rstate=None,
              verbose=False, warn=True, **kwargs):
    args['bands']=list(args['bands'])
    if len(args['bands'])!=2:
        print("If you want to analyze color curves, you need two bands!")
        sys.exit(1)

    if not curves.combined.table:
        curves.combine_curves(referenceImage=args['refImage'])

    args['curve']=_sntd_deepcopy(curves)


    if not grids:
        print('Need bounds on time delay and magnification (i.e. combinedGrids undefined)')
        sys.exit(1)


    for b in [x for x in np.unique(args['curve'].combined.table['band']) if x not in args['bands']]:
        args['curve'].combined.table=args['curve'].combined.table[args['curve'].combined.table['band']!=b]


    if not args['refModel']:
        if curves.images['image_1'].fits is not None:
            try:
                args['refModel'],bounds=create_composite_model(curves,args['refImage'],weight='logz',bound_vars=bounds.keys())
            except RuntimeError:
                args['refModel'],bounds=create_composite_model(curves,args['refImage'],bound_vars=bounds.keys())

        else:
            raise RuntimeError("Combined fit had no reference model or model name.")
    elif not isinstance(args['refModel'],sncosmo.Model):
        raise RuntimeError("Your reference model needs to be an SNCosmo model object.")


    gridBounds=dict([])
    for k in curves.images.keys():
        for par in grids.keys():
            if par=='td' and k==args['refImage']:
                continue
            if par=='mu' and k==args['refImage']:
                continue

            gridBounds[k+'_'+par]=np.array(grids[par])+args['curve'].combined.meta[par][k]
    args['curve'].color.time_delays,args['curve'].color.magnifications,args['curve'].color.time_delay_errors,args['curve'].color.magnification_errors,bestRes,bestMod= \
        nest_color_lc(args['curve'],
                      vparam_names=gridBounds.keys(),bands=args['bands'],
                      refModel=args['refModel'],bounds=gridBounds,snBounds=bounds,
                      snVparam_names=bounds.keys(),ref=args['refImage'],guess_amplitude_bound=False,maxiter=150,npoints=50)

    #fix this whole color curve thing, the problem is that right now I'm just dividing one column by the other, but
    #I need to separate it out by image first or else it'll get all mixed up...
    temp=_sntd_deepcopy(args['curve'])
    temp.color_table(args['bands'][0],args['bands'][1],time_delays=args['curve'].color.time_delays,magnifications=args['curve'].color.magnifications)


    args['curve'].color.fits=newDict()
    args['curve'].color.fits['model']=bestMod
    args['curve'].color.fits['res']=bestRes

    return args['curve']

def nest_color_lc(curves,vparam_names,bounds,snBounds,snVparam_names,ref,guess_amplitude_bound=False,bands=[],
                  minsnr=5.,refModel=None, priors=None, ppfs=None, npoints=100, method='single',
                  maxiter=None, maxcall=None, modelcov=False, rstate=None,
                  verbose=False, warn=True, **kwargs):

    # experimental parameters
    tied = kwargs.get("tied", None)

    vparam_names=list(vparam_names)
    if ppfs is None:
        ppfs = {}
    if tied is None:
        tied = {}

    # Convert bounds/priors combinations into ppfs
    if bounds is not None:
        for key, val in six.iteritems(bounds):
            if key in ppfs:
                continue  # ppfs take priority over bounds/priors
            a, b = val
            if priors is not None and key in priors:
                # solve ppf at discrete points and return interpolating
                # function
                x_samples = np.linspace(0., 1., 101)
                ppf_samples = sncosmo.utils.ppf(priors[key], x_samples, a, b)
                f = sncosmo.utils.Interp1D(0., 1., ppf_samples)
            else:
                f = sncosmo.utils.Interp1D(0., 1., np.array([a, b]))
            ppfs[key] = f

    # NOTE: It is important that iparam_names is in the same order
    # every time, otherwise results will not be reproducible, even
    # with same random seed.  This is because iparam_names[i] is
    # matched to u[i] below and u will be in a reproducible order,
    # so iparam_names must also be.

    all_delays=dict([])
    all_mus=dict([])
    all_mu_err=dict([])
    all_delay_err=dict([])
    all_mu_err[ref]=0
    all_delay_err[ref]=0
    all_delays[ref]=0
    all_mus[ref]=1

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


    global best_comb_Z
    global best_comb_Mod
    global best_comb_Res

    best_comb_Res = None
    best_comb_Mod = None
    best_comb_Z = -np.inf

    def chisquare(observed_values,expected_values,errors):
        z = (observed_values - expected_values) / errors
        chi2 = np.sum(z ** 2)

        return chi2


    def loglike(parameters):
        tempCurve=_sntd_deepcopy(curves)


        tempTds=dict([])
        tempMus=dict([])
        for i in range(len(parameters)):
            if iparam_names[i][-2:]=='td':
                tempTds[iparam_names[i][0:iparam_names[i].rfind('_')]]=parameters[i]
            elif iparam_names[i][-2:]=='mu':
                tempMus[iparam_names[i][0:iparam_names[i].rfind('_')]]=parameters[i]
        tempTds[ref]=all_delays[ref]
        tempMus[ref]=all_mus[ref]
        tempCurve.color_table(bands[0],bands[1],time_delays=tempTds,magnifications=tempMus)
        tempCurve.combine_curves(time_delays=tempTds,magnifications=tempMus,referenceImage=ref)
        zpDict={b:tempCurve.combined.table['zp'][tempCurve.combined.table['band']==b][0] for b in bands}
        zpsys=tempCurve.combined.table['zpsys'][0]
        tempRes,tempMod=_inner_color_nest_lc(tempCurve.color.table,refModel,bands,vparam_names=snVparam_names,bounds=snBounds,
                                             zp=zpDict,zpsys=zpsys,guess_amplitude_bound=False,maxiter=10,npoints=5)

        global best_comb_Res
        global best_comb_Mod
        global best_comb_Z
        if tempRes.logz>best_comb_Z:
            best_comb_Res=copy(tempRes)
            best_comb_Z=copy(tempRes.logz)
            best_comb_Mod=copy(tempMod)

        #for im in tempCurve.images.keys():
        #    tempCurve.color.table['time'][tempCurve.color.table['image']==im]-=tempTds[im]
        #    tempCurve.color.table[bands[0]+'-'+bands[1]][tempCurve.color.table['image']==im]+=2.5*np.log10(tempMus[im])

        tempColor=tempMod.color(bands[0],bands[1],'ab',tempCurve.color.table['time'])
        tempCurve.color.table=tempCurve.color.table[~np.isnan(tempColor)]
        tempColor=tempColor[~np.isnan(tempColor)]


        chisq=-chisquare(tempCurve.color.table[bands[0]+'-'+bands[1]],
                         tempColor,tempCurve.color.table[bands[0]+'-'+bands[1]+'_err'])/len(tempColor)
        if np.isnan(chisq):
            print('woah')

        return(chisq)




    res = nestle.sample(loglike, prior_transform, ndim, npdim=npdim,
                        npoints=npoints, method=method, maxiter=maxiter,
                        maxcall=maxcall, rstate=rstate,
                        callback=(nestle.print_progress if verbose else None))

    res = sncosmo.utils.Result(niter=res.niter,
                               ncall=res.ncall,
                               logz=res.logz,
                               logzerr=res.logzerr,
                               h=res.h,
                               samples=res.samples,
                               weights=res.weights,
                               logvol=res.logvol,
                               logl=res.logl,
                               vparam_names=copy(vparam_names),
                               bounds=bounds)
    pdf=_get_marginal_pdfs(res,nbins=npoints,verbose=False)
    for im in [x for x in curves.images.keys() if x!=ref]:
        all_delays[im]=pdf[im+'_td'][2]
        all_delay_err[im]=pdf[im+'_td'][3]
        all_mus[im]=pdf[im+'_mu'][2]
        all_mu_err[im]=pdf[im+'_mu'][3]

    return all_delays,all_mus,all_delay_err,all_mu_err,best_comb_Res,best_comb_Mod

def _inner_color_nest_lc(data, model, bands,vparam_names, bounds,zp, zpsys,guess_amplitude_bound=False,
                         minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
                         maxiter=None, maxcall=None, modelcov=False, rstate=None,
                         verbose=False, warn=True, **kwargs):


    try:
        import nestle
    except ImportError:
        raise ImportError("nest_lc() requires the nestle package.")

    # warnings
    if "nobj" in kwargs:
        warnings.warn("The nobj keyword is deprecated and will be removed in "
                      "sncosmo v2.0. Use `npoints` instead.")
        npoints = kwargs.pop("nobj")

    # experimental parameters
    tied = kwargs.get("tied", None)


    sortidx=None
    model = copy(model)
    bounds = copy(bounds)  # need to copy this b/c we modify it below

    # Order vparam_names the same way it is ordered in the model:
    vparam_names = [s for s in model.param_names if s in vparam_names]



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


    if ppfs is None:
        ppfs = {}
    if tied is None:
        tied = {}

    # Convert bounds/priors combinations into ppfs
    if bounds is not None:
        for key, val in six.iteritems(bounds):
            if key in ppfs:
                continue  # ppfs take priority over bounds/priors
            a, b = val
            if priors is not None and key in priors:
                # solve ppf at discrete points and return interpolating
                # function
                x_samples = np.linspace(0., 1., 101)
                ppf_samples = sncosmo.utils.ppf(priors[key], x_samples, a, b)
                f = sncosmo.utils.Interp1D(0., 1., ppf_samples)
            else:
                f = sncosmo.utils.Interp1D(0., 1., np.array([a, b]))
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

    def chisq(data, model,zp,zpsys):
        mCol=model.color(bands[0],bands[1],zpsys,data['time'])
        data=data[~np.isnan(mCol)]
        mCol=mCol[~np.isnan(mCol)]
        diff = data[bands[0]+'-'+bands[1]]- mCol
        chi=diff/data[bands[0]+'-'+bands[1]+'_err']
        chi2 = np.sum(chi ** 2)
        return chi2

    def loglike(parameters):

        model.parameters[idx] = parameters

        chisq_res=chisq(data,model,zp,zpsys)

        return -0.5 * chisq_res

    res = nestle.sample(loglike, prior_transform, ndim, npdim=npdim,
                        npoints=npoints, method=method, maxiter=maxiter,
                        maxcall=50, rstate=rstate,
                        callback=(nestle.print_progress if verbose else None))

    # estimate parameters and covariance from samples
    vparameters, cov = nestle.mean_and_cov(res.samples, res.weights)

    # update model parameters to estimated ones.
    model.set(**dict(zip(vparam_names, vparameters)))


    # `res` is a nestle.Result object. Collect result into a sncosmo.Result
    # object for consistency, and add more fields.
    res = sncosmo.utils.Result(niter=res.niter,
                               ncall=res.ncall,
                               logz=res.logz,
                               logzerr=res.logzerr,
                               h=res.h,
                               samples=res.samples,
                               weights=res.weights,
                               logvol=res.logvol,
                               logl=res.logl,
                               vparam_names=copy(vparam_names),
                               ndof=len(data) - len(vparam_names),
                               bounds=bounds,
                               parameters=model.parameters.copy(),
                               covariance=cov,
                               errors=OrderedDict(zip(vparam_names,
                                                      np.sqrt(np.diagonal(cov)))),
                               param_dict=OrderedDict(zip(model.param_names,
                                                          model.parameters)))

    return res, model

def _fitCombined(curves,mods,args,bounds,grids,guess_amplitude_bound=False,
                 minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
                 maxiter=None, maxcall=None, modelcov=False, rstate=None,
                 verbose=False, warn=True, **kwargs):

    args['bands']=list(args['bands'])
    if not curves.combined.table:
        curves.combine_curves(referenceImage=args['refImage'])
    args['curve']=_sntd_deepcopy(curves)
    if not grids:
        print('Need bounds on time delay and magnification (i.e. combinedGrids undefined)')
        sys.exit(1)


    for b in [x for x in np.unique(args['curve'].combined.table['band']) if x not in args['bands']]:
        args['curve'].combined.table=args['curve'].combined.table[args['curve'].combined.table['band']!=b]


    if not args['refModel']:
        if curves.images['image_1'].fits is not None:
            try:
                args['refModel'],bounds=create_composite_model(curves,args['refImage'],weight='logz',bound_vars=bounds.keys())
            except RuntimeError:
                args['refModel'],bounds=create_composite_model(curves,args['refImage'],bound_vars=bounds.keys())

        else:
            raise RuntimeError("Combined fit had no reference model or model name.")
    elif not isinstance(args['refModel'],sncosmo.Model):
        raise RuntimeError("Your reference model needs to be an SNCosmo model object.")

    gridBounds=dict([])
    for k in curves.images.keys():
        for par in grids.keys():
            if par=='td' and k==args['refImage']:
                continue
            if par=='mu' and k==args['refImage']:
                continue
            gridBounds[k+'_'+par]=np.array(grids[par])+args['curve'].combined.meta[par][k]

    myRes=dict([])
    for b in args['bands']:
        myRes[b]=nest_combined_lc(args['curve'],
                                  vparam_names=gridBounds.keys(),
                                  band=b,refModel=args['refModel'],bounds=gridBounds,snBounds=bounds,
                                  snVparam_names=bounds.keys(),ref=args['refImage'],guess_amplitude_bound=True,maxiter=150,npoints=50)

    weights=np.array([myRes[b][-2].logz for b in args['bands']])
    final_params=dict([])

    for param in myRes[args['bands'][0]][-1].param_names:
        final_params[param]=np.average([myRes[b][-1].get(param) for b in args['bands']],weights=weights)

    args['curve'].combined.time_delays=dict([])
    args['curve'].combined.magnifications=dict([])
    args['curve'].combined.magnification_errors=dict([])
    args['curve'].combined.time_delay_errors=dict([])

    for k in curves.images.keys():
        args['curve'].combined.time_delays[k]=np.average([myRes[b][0][k] for b in args['bands']],weights=weights)
        args['curve'].combined.magnifications[k]=np.average([myRes[b][1][k] for b in args['bands']],weights=weights)
    args['curve'].combined.time_delay_errors[k]=myRes[args['bands'][np.where(weights==np.max(weights))[0][0]]][2]
    args['curve'].combined.magnification_errors[k]=myRes[args['bands'][np.where(weights==np.max(weights))[0][0]]][3]
    bestRes=myRes[args['bands'][np.where(weights==np.max(weights))[0][0]]][4]
    bestMod=myRes[args['bands'][np.where(weights==np.max(weights))[0][0]]][5]

    bestMod.set(**final_params)

    args['curve'].combine_curves(time_delays=args['curve'].combined.time_delays,magnifications=args['curve'].combined.magnifications,referenceImage=args['refImage'])

    for b in [x for x in set(args['curve'].combined.table['band']) if x not in args['bands']]:
        args['curve'].combined.table=args['curve'].combined.table[args['curve'].combined.table['band']!=b]

    args['curve'].combined.fits=newDict()
    args['curve'].combined.fits['model']=bestMod
    args['curve'].combined.fits['res']=bestRes

    return args['curve']


def create_composite_model(curves,ref,weight='chisq',bound_vars=[]):
    weights=np.array([curves.images[im].fits.res[weight] for im in np.sort(list(curves.images.keys()))])
    if weight=='chisq':
        weights=1./weights
    final_vars=dict([])
    bounds=dict([])
    for param in curves.images[ref].fits.model.param_names:
        if param not in bound_vars:
            final_vars[param]=curves.images[ref].fits.model.get(param)
        else:
            final_vars[param]=np.average([curves.images[k].fits.model.get(param) for k in np.sort(list(curves.images.keys()))],
                                         weights=weights)
        if param in bound_vars:
            bounds[param]=np.average([curves.images[k].fits.final_errs[param] for k in np.sort(list(curves.images.keys()))],
                                     weights=weights)/np.sqrt(len(weights))
            bounds[param]=(final_vars[param]-bounds[param],final_vars[param]+bounds[param])

            final_mod=copy(curves.images[list(curves.images.keys())[np.where(weights==np.max(weights))[0][0]]].fits.model)
    final_mod.set(**final_vars)
    return(final_mod,bounds)



def nest_combined_lc(curves,vparam_names,bounds,snBounds,snVparam_names,ref,guess_amplitude_bound=False,
                     minsnr=5.,refModel=False,band=None, priors=None, ppfs=None, npoints=100, method='single',
                     maxiter=None, maxcall=None, modelcov=False, rstate=None,
                     verbose=False, warn=True, **kwargs):

    # experimental parameters
    tied = kwargs.get("tied", None)

    vparam_names=list(vparam_names)
    if ppfs is None:
        ppfs = {}
    if tied is None:
        tied = {}

    # Convert bounds/priors combinations into ppfs
    if bounds is not None:
        for key, val in six.iteritems(bounds):
            if key in ppfs:
                continue  # ppfs take priority over bounds/priors
            a, b = val
            if priors is not None and key in priors:
                # solve ppf at discrete points and return interpolating
                # function
                x_samples = np.linspace(0., 1., 101)
                ppf_samples = sncosmo.utils.ppf(priors[key], x_samples, a, b)
                f = sncosmo.utils.Interp1D(0., 1., ppf_samples)
            else:
                f = sncosmo.utils.Interp1D(0., 1., np.array([a, b]))
            ppfs[key] = f

    # NOTE: It is important that iparam_names is in the same order
    # every time, otherwise results will not be reproducible, even
    # with same random seed.  This is because iparam_names[i] is
    # matched to u[i] below and u will be in a reproducible order,
    # so iparam_names must also be.

    all_delays=dict([])
    all_mus=dict([])
    all_mu_err=dict([])
    all_delay_err=dict([])
    all_mu_err[ref]=0
    all_delay_err[ref]=0
    all_delays[ref]=0
    all_mus[ref]=1

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


    global best_comb_Z
    global best_comb_Mod
    global best_comb_Res

    best_comb_Res = None
    best_comb_Mod = None
    best_comb_Z = -np.inf


    def loglike(parameters):
        tempCurve=_sntd_deepcopy(curves)

        tempTds=dict([])
        tempMus=dict([])
        for i in range(len(parameters)):
            if iparam_names[i][-2:]=='td':
                tempTds[iparam_names[i][0:iparam_names[i].rfind('_')]]=parameters[i]
            elif iparam_names[i][-2:]=='mu':
                tempMus[iparam_names[i][0:iparam_names[i].rfind('_')]]=parameters[i]
        tempTds[ref]=all_delays[ref]
        tempMus[ref]=all_mus[ref]

        tempCurve.combine_curves(time_delays=tempTds,magnifications=tempMus,referenceImage=ref)

        tempCurve.combined.table=tempCurve.combined.table[tempCurve.combined.table['band']==band]
        tempRes,tempMod=sncosmo.nest_lc(tempCurve.combined.table,refModel,
                                        vparam_names=snVparam_names,bounds=snBounds,guess_amplitude_bound=False,maxiter=10,npoints=5)
        global best_comb_Res
        global best_comb_Mod
        global best_comb_Z
        if tempRes.logz>best_comb_Z:
            best_comb_Res=copy(tempRes)
            best_comb_Z=copy(tempRes.logz)
            best_comb_Mod=copy(tempMod)
        return(tempRes.logz)




    res = nestle.sample(loglike, prior_transform, ndim, npdim=npdim,
                        npoints=npoints, method=method, maxiter=maxiter,
                        maxcall=maxcall, rstate=rstate,
                        callback=(nestle.print_progress if verbose else None))

    res = sncosmo.utils.Result(niter=res.niter,
                               ncall=res.ncall,
                               logz=res.logz,
                               logzerr=res.logzerr,
                               h=res.h,
                               samples=res.samples,
                               weights=res.weights,
                               logvol=res.logvol,
                               logl=res.logl,
                               vparam_names=copy(vparam_names),
                               bounds=bounds)
    pdf=_get_marginal_pdfs(res,nbins=npoints,verbose=False)
    for im in [x for x in curves.images.keys() if x!=ref]:
        all_delays[im]=pdf[im+'_td'][2]
        all_delay_err[im]=pdf[im+'_td'][3]
        all_mus[im]=pdf[im+'_mu'][2]
        all_mu_err[im]=pdf[im+'_mu'][3]

    return all_delays,all_mus,all_delay_err,all_mu_err,best_comb_Res,best_comb_Mod




def _fitSeparate(curves,mods,args,bounds,npoints=100,maxiter=None,**kwargs):
    resList=dict([])
    fitDict=dict([])
    if 't0' in args['bounds']:
        t0Bounds=copy(args['bounds']['t0'])
    if 'amplitude' in args['bounds']:
        ampBounds=copy(args['bounds']['amplitude'])
    for d in curves.images.keys():
        #print(curves.images[d].simMeta)
        args['curve']=copy(curves.images[d])
        for b in [x for x in np.unique(args['curve'].table['band']) if x not in args['bands']]:
            args['curve'].table=args['curve'].table[args['curve'].table['band']!=b]

        if 't0' in args['bounds']:
            if args['t0_guess'] is not None:
                args['bounds']['t0']=(t0Bounds[0]+args['t0_guess'][d],t0Bounds[1]+args['t0_guess'][d])
            else:
                maxFlux=np.max(args['curve'].table['flux'])
                maxTime=args['curve'].table['time'][args['curve'].table['flux']==maxFlux]
                args['bounds']['t0']=(t0Bounds[0]+maxTime,t0Bounds[1]+maxTime)

        if 'amplitude' in args['bounds'] and args['guess_amplitude']:
            args['bounds']['amplitude']=(ampBounds[0]*np.max(args['curve'].table['flux']),ampBounds[1]*np.max(args['curve'].table['flux']))



        curves.images[d].fits=newDict()
        if True or len(args['curve'].table)>63 or len(mods)==1 or args['snType']=='Ia':
            fits=[]
            for mod in mods:
                if mod =='BazinSource' or isinstance(mod,BazinSource):
                    fits.append(param_fit(args,mod))


                else:
                    if len(mods)==1:
                        doFit=False
                    else:
                        doFit=True
                    args['doFit']=doFit
                    fits.append(_fit_data_wrap((mod,args)))
        else:
            args['doFit']=True
            fits=pyParz.foreach(mods,_fit_data,args)

        if len(fits)>1:
            bestChisq=np.inf
            for f in fits:
                if f:
                    res=f['res']
                    mod=f['model']
                    if res.chisq <bestChisq:
                        bestChisq=res.chisq
                        bestFit=mod
                        bestRes=res
        else:
            bestFit=fits[0]['model']
            bestRes=fits[0]['res']


        fitDict[d]=[fits,bestFit,bestRes]

    if not all([fitDict[d][1]._source.name==fitDict[list(fitDict.keys())[0]][1]._source.name for d in fitDict.keys()]):
        print('All models did not match, finding best...')
        bestChisq=np.inf
        bestMod=None
        for d in fitDict.keys():
            chisq=fitDict[d][2].chisq
            if chisq<bestChisq:
                bestChisq=chisq
                bestMod=fitDict[d][1]._source.name

        for d in fitDict.keys():
            for f in fitDict[d][0]:

                if f and f['model']._source.name==bestMod:
                    fitDict[d][1]=f['model']
                    fitDict[d][2]=f['res']
                    break
    dofs=dict([])

    for d in fitDict.keys():
        _,bestFit,bestMod=fitDict[d]
        tempTable=deepcopy(curves.images[d].table)
        for b in [x for x in np.unique(tempTable['band']) if x not in args['bands']]:
            tempTable=tempTable[tempTable['band']!=b]
        if args['flip']:
            tempTable['flux']=np.flip(tempTable['flux'],axis=0)

        if 't0' in args['bounds']:
            if args['t0_guess'] is not None:
                args['bounds']['t0']=(t0Bounds[0]+args['t0_guess'][d],t0Bounds[1]+args['t0_guess'][d])
            else:
                maxFlux=np.max(tempTable['flux'])
                maxTime=tempTable['time'][tempTable['flux']==maxFlux]
                args['bounds']['t0']=(t0Bounds[0]+maxTime,t0Bounds[1]+maxTime)
        if args['snType']=='Ia' and 'x1' not in args['bounds']:
            args['bounds']['x1']=(-1,1)
        if 'amplitude' in args['bounds'] and args['guess_amplitude']:
            args['bounds']['amplitude']=(ampBounds[0]*np.max(tempTable['flux']),ampBounds[1]*np.max(tempTable['flux']))
        if 'amplitude' not in args['bounds']:
            guess_amp_bounds=True
        else:
            guess_amp_bounds=False

        nest_res,nest_fit=_nested_wrapper(curves,tempTable,bestFit,vparams=bestRes.vparam_names,bounds=args['bounds'],
                                          guess_amplitude_bound=guess_amp_bounds,microlensing=args['microlensing'],
                                          zpsys=curves.images[d].zpsys,kernel=args['kernel'],maxiter=maxiter,npoints=npoints,nsamples=args['nMicroSamples'])




        if nest_res.ndof != len(tempTable)- len(nest_res.vparam_names):
            dofs[d]=len(tempTable)- len(nest_res.vparam_names)
        else:
            dofs[d]=nest_res.ndof


        resList[d]=nest_res
        curves.images[d].fits=newDict()
        curves.images[d].fits['model']=nest_fit
        curves.images[d].fits['res']=nest_res


    joint=_joint_likelihood(resList,verbose=False)

    for d in np.sort(list(curves.images.keys())):
        errs=dict([])
        if 'micro' in curves.images[d].fits.res.errors.keys():
            errs['micro']=curves.images[d].fits.res.errors['micro']
        else:
            errs['micro']=0.
        tempTable=copy(curves.images[d].table)

        for b in [x for x in np.unique(tempTable['band']) if x not in args['bands']]:
            tempTable=tempTable[tempTable['band']!=b]
        if not np.all([x in __thetaL__ for x in joint.keys()]):# or args['microlensing'] is not None:
            bds=dict([])

            final_vparams=[]
            for p in joint.keys():
                if p in curves.images[d].fits.res.param_names:
                    if isinstance(joint[p],dict):
                        if p=='t0':
                            bds[p]=(-3+joint[p][d][0],3+joint[p][d][0])
                        elif p!='amplitude' and p!='x0':

                            if joint[p][d][1]==0 or np.round(joint[p][d][1]/joint[p][d][0],6)==0:
                                bds[p]=(joint[p][d][0]-.05*joint[p][d][0],joint[p][d][0]+.05*joint[p][d][0])
                            else:
                                bds[p]=(joint[p][d][0]-joint[p][d][1],joint[p][d][0]+joint[p][d][1])
                        errs[p]=joint[p][d][1]
                        final_vparams.append(p)



                    else:
                        curves.images[d].fits.model.set(**{p:joint[p][0]})
                        errs[p]=joint[p][1]

            finalRes,finalFit=sncosmo.nest_lc(tempTable,curves.images[d].fits.model,final_vparams,bounds=bds,guess_amplitude_bound=True,maxiter=None)

            finalRes.ndof=dofs[d]
            curves.images[d].fits=newDict()
            curves.images[d].fits['model']=finalFit
            curves.images[d].fits['res']=finalRes

        else:
            for p in joint.keys():
                if p in curves.images[d].fits.res.param_names:
                    if isinstance(joint[p],dict):
                        curves.images[d].fits.model.set(**{p:joint[p][d][0]})
                        errs[p]=joint[p][d][1]

                    else:
                        curves.images[d].fits.model.set(**{p:joint[p][0]})
                        errs[p]=joint[p][1]
        curves.images[d].fits['final_errs']=errs

    tds,td_errs,mags,mag_errs,times,fluxes,time_errors,flux_errors=timeDelaysAndMagnifications(curves)
    curves.time_delays=tds
    curves.time_delay_errors=td_errs
    curves.magnifications=mags
    curves.magnification_errors=mag_errs
    curves.measurements={'t0':times,'A':fluxes,'t0_err':time_errors,'A_err':flux_errors}

    if args['showPlots']:
        for d in curves.images.keys():
            tempTable=copy(curves.images[d].table)
            for b in [x for x in np.unique(tempTable['band']) if x not in args['bands']]:
                tempTable=tempTable[tempTable['band']!=b]
            tempMod=copy(curves.images[d].fits.model)


            sncosmo.plot_lc(tempTable,model=tempMod)

            #plt.savefig(nest_fit._source.name+'_'+tempTable['band'][0]+'_refs_'+d+'.pdf',format='pdf',overwrik4ite=True)
            plt.savefig('example_plot_dust_image_'+str(d[-1])+'.png',format='png',overwrite=True)
            plt.show()
            plt.clf()
            plt.close()

    return curves

def _micro_uncertainty(args):
    sample,other=args
    nest_fit,data,colnames,x_pred,vparam_names,bounds=other
    data=Table(data,names=colnames)

    temp_nest_mod=deepcopy(nest_fit)
    tempMicro=AchromaticMicrolensing(x_pred/(1+nest_fit.get('z')),sample,magformat='multiply')
    temp_nest_mod.add_effect(tempMicro,'microlensing','rest')
    tempRes,tempMod=sncosmo.nest_lc(data,temp_nest_mod,vparam_names=vparam_names,bounds=bounds,guess_amplitude_bound=True,maxiter=None,npoints=200)

    return float(tempMod.get('t0'))


def _nested_wrapper(curves,data,model,vparams,bounds,guess_amplitude_bound,microlensing,zpsys,kernel,maxiter,npoints,nsamples):

    temp=deepcopy(data)
    vparam_names=deepcopy(vparams)


    if microlensing is not None:
        nest_res,nest_fit=sncosmo.nest_lc(temp,model,vparam_names=vparam_names,bounds=bounds,guess_amplitude_bound=guess_amplitude_bound,maxiter=maxiter,npoints=npoints)




        micro,sigma,x_pred,y_pred,samples=fit_micro(curves,nest_res,nest_fit,temp,zpsys,nsamples,micro_type=microlensing,kernel=kernel)


        temp=deepcopy(data)

        t0s=pyParz.foreach(samples.T,_micro_uncertainty,[nest_fit,np.array(temp),temp.colnames,x_pred,vparam_names,bounds])
        mu,sigma=scipy.stats.norm.fit(t0s)

        nest_res.errors['micro']=np.sqrt(np.abs(nest_fit.get('t0')-mu)**2+(3*sigma)**2)
        bestRes=nest_res
        bestMod=nest_fit


    else:
        bestRes,bestMod=sncosmo.nest_lc(data,model,vparam_names=vparam_names,bounds=bounds,guess_amplitude_bound=guess_amplitude_bound,maxiter=maxiter,npoints=npoints)

    return(bestRes,bestMod)

def _maxFromModel(mod,band,zp,zpsys):
    time=np.arange(mod.mintime(),mod.maxtime(),.1)
    flux=mod.bandflux(band,time,zp,zpsys)
    return (time[flux==np.max(flux)],np.max(flux))

def fit_micro(curves,res,fit,dat,zpsys,nsamples,micro_type='achromatic',kernel='RBF'):
    t0=fit.get('t0')
    fit.set(t0=t0)
    data=deepcopy(dat)

    data['time']-=t0
    data=data[data['time']<=40.]
    data=data[data['time']>=-15.]
    achromatic=micro_type.lower()=='achromatic'
    if achromatic:
        allResid=[]
        allErr=[]
        allTime=[]
    else:
        allResid=dict([])
        allErr=dict([])
        allTime=dict([])
    for b in np.unique(data['band']):
        tempData=data[data['band']==b]
        tempData=tempData[tempData['flux']>.1]
        tempTime=tempData['time']
        mod=fit.bandflux(b,tempTime+t0,zpsys=zpsys,zp=tempData['zp'])
        #fig=plt.figure()
        #ax=fig.gca()
        #ax.plot(tempTime,mod)
        #ax.scatter(tempData['time'],tempData['flux'])
        #plt.show()
        tempData=tempData[mod>.1]
        residual=tempData['flux']/mod[mod>.1]
        tempData=tempData[~np.isnan(residual)]
        residual=residual[~np.isnan(residual)]
        tempTime=tempTime[~np.isnan(residual)]

        if achromatic:
            allResid=np.append(allResid,residual)
            allErr=np.append(allErr,residual*tempData['fluxerr']/tempData['flux'])
            allTime=np.append(allTime,tempTime)
        else:
            allResid[b]=residual
            allErr[b]=residual*tempData['fluxerr']/tempData['flux']
            allTime[b]=tempTime

    if kernel=='RBF':
        kernel = RBF(10., (20., 50.))


    if achromatic:

        gp = GaussianProcessRegressor(kernel=kernel, alpha=allErr ** 2,
                                      n_restarts_optimizer=100)

        try:
            gp.fit(np.atleast_2d(allTime).T,allResid.ravel())
        except:
            temp=np.atleast_2d(allTime).T
            temp2=allResid.ravel()
            temp=temp[~np.isnan(temp2)]
            temp2=temp2[~np.isnan(temp2)]
            gp.fit(temp,temp2)



        X=np.atleast_2d(np.linspace(np.min(allTime), np.max(allTime), 1000)).T

        y_pred, sigma = gp.predict(X, return_std=True)
        samples=gp.sample_y(X,nsamples)


        '''
        plt.close()
        fig=plt.figure()
        ax=fig.gca()
        for i in range(samples.shape[1]):
            if i==0:
                ax.plot(X, samples[:,i],alpha=.1,label='Posterior Samples',color='b')
            else:
                ax.plot(X, samples[:,i],alpha=.1,color='b')
        ax.errorbar(allTime.ravel(), allResid, allErr, fmt='r.', markersize=10, label=u'Observations')
        print(len(allResid))
        ax.plot(X, y_pred - 3 * sigma, '--g')
        ax.plot(X, y_pred + 3 * sigma, '--g',label='$3\sigma$ Bounds')
        ax.plot(X,y_pred,'k-.',label="GPR Prediction")

        ax.set_ylabel('Magnification ($\mu$)')
        ax.set_xlabel('Observer Frame Time (Days)')
        ax.plot(X,curves.images['S2'].simMeta['microlensing_params'](X/(1+1.33))/np.median(curves.images['S2'].simMeta['microlensing_params'](X/(1+1.33))),'k',label='True $\mu$-Lensing')
        ax.legend(fontsize=10)
        plt.savefig('microlensing_gpr.pdf',format='pdf',overwrite=True)
        plt.close()
        sys.exit()
        '''

        tempX=X[:,0]
        tempX=np.append([fit._source._phase[0]*(1+fit.get('z'))],np.append(tempX,[fit._source._phase[-1]*(1+fit.get('z'))]))
        y_pred=np.append([1.],np.append(y_pred,[1.]))
        sigma=np.append([0.],np.append(sigma,[0.]))

        result=AchromaticMicrolensing(tempX/(1+fit.get('z')),y_pred,magformat='multiply')

        '''
        fig=plt.figure()
        ax=fig.gca()
        #plt.plot(X, resid, 'r:', label=u'$f(x) = x\,\sin(x)$')
        ax.errorbar(allTime.ravel(), allResid, allErr, fmt='r.', markersize=10, label=u'Observations')
        #for c in np.arange(0,1.01,.1):
        #    for d in np.arange(-np.max(sigma),np.max(sigma),np.max(sigma)/10):
        #        plt.plot(X, 1+(y_pred[1:-1]-1)*c+d, 'b-')#, label=u'Prediction')
        ax.plot(X, y_pred[1:-1], 'b-' ,label=u'Prediction')

        ax.fill(np.concatenate([X, X[::-1]]),
                 np.concatenate([y_pred[1:-1] - 1.9600 * sigma[1:-1],
                                 (y_pred[1:-1] + 1.9600 * sigma[1:-1])[::-1]]),
                 alpha=.5, fc='b', ec='None', label='95% confidence interval')
        ax.set_xlabel('$x$')
        ax.set_ylabel('$f(x)$')
        #plt.xlim(-10, 50)
        #plt.ylim(.8, 1.4)

        ax.legend(loc='upper left')
        #plt.show()
        #plt.show()
        #figures.append(ax)
        #plt.clf()
        plt.close()
        '''


    else:
        pass
        #TODO make chromatic microlensing a thing



    return result,sigma,X[:,0],y_pred[1:-1],samples

def timeDelaysAndMagnifications(curves):
    times=dict([])
    fluxes=dict([])
    time_errors=dict([])
    flux_errors=dict([])

    for d in curves.images.keys():

        b=[x for x in curves.bands if curves.images[d].fits.model.bandoverlap(x)]
        if len(b)==0:
            print('None of your bands overlap your fit?')
            sys.exit()
        b=b[0]
        try:
            timeOfPeak=curves.images[d].fits.model.get('t0')
            band_time_errors=curves.images[d].fits.final_errs['t0']

        except:
            timeOfPeak,peakFlux=_maxFromModel(curves.images[d].fits.model,b,
                                              curves.images[d].table['zp'][curves.images[d]['band']==b],
                                              curves.images[d].zpsys)
            band_time_errors=0
            band_flux_errors=0

        for amp in ['x0','amplitude','A']:
            try:
                peakFlux=curves.images[d].fits.model.get(amp)
                band_flux_errors=curves.images[d].fits.final_errs[amp]
                success=True
                break
            except:
                success=False

        if not success:
            peakFlux=curves.images[d].fits.model.bandflux(b,timeOfPeak,
                                                          zp=curves.images[d].table['zp'][curves.images[d]['band']==b],
                                                          zpsys=curves.images[d].zpsys)
            band_flux_errors=0


        band_times=timeOfPeak
        band_fluxes=peakFlux
        times[d]=band_times
        fluxes[d]=band_fluxes
        time_errors[d]=band_time_errors
        flux_errors[d]=band_flux_errors
    ims=np.sort(list(curves.images.keys()))
    delays=dict([])
    mags=dict([])
    delay_errs=dict([])
    mag_errs=dict([])
    ref1=None
    ref2=None
    ref1_err=None
    ref2_err=None
    for im in ims:
        if not ref1:
            ref1=times[im]
            delays[im]=0
            ref2=fluxes[im]
            mags[im]=1
            ref1_err=time_errors[im]
            ref2_err=flux_errors[im]
            delay_errs[im]=0
            mag_errs[im]=0

        else:
            delays[im]=times[im]-ref1
            mags[im]=fluxes[im]/ref2
            delay_errs[im]=math.sqrt(time_errors[im]**2+ref1_err**2)

            mag_errs[im]=mags[im]*math.sqrt((flux_errors[im]/fluxes[im])**2+(ref2_err/ref2)**2)

    return (delays,delay_errs,mags,mag_errs,times,fluxes,time_errors,flux_errors)



def _plot_marginal_pdfs( res, nbins=101, **kwargs):
    """ plot the results of a classification run
    :return:
    """
    from matplotlib import pyplot as pl
    import numpy as np

    nparam = len(res.vparam_names)
    # nrow = np.sqrt( nparam )
    # ncol = nparam / nrow + 1
    nrow, ncol = 1, nparam

    pdfdict = _get_marginal_pdfs( res, nbins )

    fig = plt.gcf()
    for parname in res.vparam_names :
        iax = res.vparam_names.index( parname )+1
        ax = fig.add_subplot( nrow, ncol, iax )

        parval, pdf, mean, std = pdfdict[parname]
        ax.plot(  parval, pdf, **kwargs )
        if np.abs(std)>=0.1:
            ax.text( 0.95, 0.95, '%s  %.1f +- %.1f'%( parname, np.round(mean,1), np.round(std,1)),
                     ha='right',va='top',transform=ax.transAxes )
        elif np.abs(std)>=0.01:
            ax.text( 0.95, 0.95, '%s  %.2f +- %.2f'%( parname, np.round(mean,2), np.round(std,2)),
                     ha='right',va='top',transform=ax.transAxes )
        elif np.abs(std)>=0.001:
            ax.text( 0.95, 0.95, '%s  %.3f +- %.3f'%( parname, np.round(mean,3), np.round(std,3)),
                     ha='right',va='top',transform=ax.transAxes )
        else :
            ax.text( 0.95, 0.95, '%s  %.3e +- %.3e'%( parname, mean, std),
                     ha='right',va='top',transform=ax.transAxes )

    plt.draw()

def _fit_data_wrap(args):
    try:
        return _fit_data(args)
    except RuntimeError:
        print('There was an issue running model {0}, skipping...'.format(args[0]))
        return None#args[0]

#todo decide about multiple versions of model
def _fit_data(args):
    """
    Helper function that allows parallel processing to occur.
    :param args: All the arguments given from fit_data
    :return: modResults: tuple (Name of current model (including version if exists),name of the sncosmo model,version
                    of sncosmo model,list of tuples containing index and dcurve.fit[modname] for each dcurve in curves)
    """
    warnings.simplefilter("ignore")
    mod=args[0]

    args=args[1]
    if isinstance(mod, tuple):
        version = mod[1]
        mod = mod[0]
    else:
        version = None
    dust_dict={'SFD98Map':sncosmo.SFD98Map,'CCM89Dust':sncosmo.CCM89Dust,'OD94Dust':sncosmo.OD94Dust,'F99Dust':sncosmo.F99Dust}
    if args['dust']:
        dust=dust_dict[args['dust']]()
    else:
        dust=[]
    effect_names=args['effect_names']
    effect_frames=args['effect_frames']
    effects=[dust for i in range(len(effect_names))] if effect_names else []
    effect_names=effect_names if effect_names else []
    effect_frames=effect_frames if effect_frames else []
    if not isinstance(effect_names,(list,tuple)):
        effects=[effect_names]
    if not isinstance(effect_frames,(list,tuple)):
        effects=[effect_frames]
    if isinstance(mod,str):
        modName=mod+'_'+version if version else deepcopy(mod)
    else:
        modName=mod.name+'_'+version if version else deepcopy(mod)

    if isinstance(mod,str) or isinstance(mod,sncosmo.Source):
        source=sncosmo.get_source(mod)
        smod = sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
    else:
        smod=mod
    params=args['params'] if args['params'] else [x for x in smod.param_names]
    fits=newDict()
    dcurve=args['curve']
    #if not np.any([smod.bandoverlap(band) for band in dcurve.bands]):
    #    raise RuntimeError("No band overlap for model %s"%modName)
    fits.method=args['fitting_method']
    fits.bounds=args['bounds'] if args['bounds'] else {}
    fits.ignore=args['ignore'] if args['ignore'] else []
    fits.constants = args['constants'] if args['constants'] else {x: y for x, y in zip(dcurve.meta.keys(),dcurve.meta.values()) if x != 'info'}


    no_bound = {x for x in params if x in _needs_bounds and x not in fits.bounds.keys() and x not in fits.constants.keys()}
    if no_bound:
        params=list(set(params)-no_bound)
    params= [x for x in params if x not in fits.ignore and x not in fits.constants.keys()]
    fits.params = params
    if fits.constants is not None:
        try:
            smod.set(**fits.constants)
        except:
            raise RuntimeError('You may have some parameters in "constants" that are not in your model.')

    if args['doFit']:
        if args['fitting_method']=='mcmc':
            fits.res, fits.model = args['sn_func'][args['fitting_method']](dcurve.table, smod, fits.params, fits.bounds, **args['props'])
        elif args['fitting_method']=='nest':
            fits.res, fits.model = args['sn_func'][args['fitting_method']](dcurve.table, smod, fits.params, fits.bounds,guess_amplitude_bound=True, verbose=False, **args['props'])
        else:
            fits.res, fits.model = args['sn_func'][args['fitting_method']](dcurve.table, smod, fits.params, fits.bounds,verbose=False, **args['props'])
        return(pyParz.parReturn(fits))
    else:
        fits.model=smod
        fits.res=newDict()
        fits.res['vparam_names']=fits.params
        return (fits)


def _joint_likelihood(resList,verbose=False):


    params=[]
    for res in resList.values():
        params=list(set(np.append(params,res.vparam_names)))
    otherParams=[x for x in params if x in __thetaL__ ]
    snparams=[x for x in params if x in __thetaSN__ or x[-1] =='z' or 'ebv' in x or 'r_v' in x]

    outDict={p:dict([]) for p in otherParams}
    for param in params:
        if verbose:
            print(param)
        probsList=[]
        weightsList=[]
        zweights=[]
        testList=[]
        for k in resList.keys():
            if param not in resList[k].vparam_names:
                continue
            res=resList[k]
            pdf=_get_marginal_pdfs(res,nbins=100,verbose=False)
            if param in snparams:
                weightsList=np.append(weightsList,pdf[param][1]*(1/np.abs(pdf[param][4])))
                testList=np.append(testList,pdf[param][1])
                zweights=np.append(zweights,1/np.abs(pdf[param][4]))
                probsList=np.append(probsList,pdf[param][0])
            elif param in otherParams:
                if verbose:
                    print( 'Image %s:  <%s> = %.6e +- %.2e'%(k, param, pdf[param][2], pdf[param][3]) )
                outDict[param][k]=(pdf[param][2],pdf[param][3])
        if param in otherParams:
            continue

        expectedValue=((weightsList*probsList).sum())/(zweights.sum())
        std = np.sqrt( ((weightsList * (probsList-expectedValue)**2 ).sum())/(zweights.sum()) )
        if verbose:
            print( '  <%s> = %.3e +- %.3e'%( param, expectedValue, std) )

        outDict[param]=(expectedValue,std)
    return(outDict)

def _get_marginal_pdfs( res, nbins=51, verbose=True ):
    """ Given the results <res> from a nested sampling chain, determine the
    marginalized posterior probability density functions for each of the
    parameters in the model.
    :param res:  the results of a nestlc run
    :param nbins: number of bins (steps along the x axis) for sampling
       each parameter's marginalized posterior probability
    :return: a dict with an entry for each parameter, giving a 2-tuple containing
       NDarrays of length nbins.  The first array in each pair gives the parameter
       value that defines the left edge of each bin along the parameter axis.
       The second array gives the posterior probability density integrated
       across that bin.
    """
    vparam_names = res.vparam_names
    weights = res.weights
    samples = res.samples

    pdfdict = {}

    for param in vparam_names :
        ipar = vparam_names.index( param )
        paramvals = samples[:,ipar]

        if nbins>1:
            if param in res.bounds :
                parvalmin, parvalmax = res.bounds[param]
            else :
                parvalmin, parvalmax = 0.99*paramvals.min(), 1.01*paramvals.max()
            parambins = np.linspace( parvalmin, parvalmax, nbins, endpoint=True ).flatten()
            binindices = np.digitize( paramvals, parambins )

            # we estimate the marginalized pdf by summing the weights of all points in the bin,
            # where the weight of each point is the prior volume at that point times the
            # likelihood, divided by the total evidence
            pdf = np.array( [ weights[np.where( binindices==ibin )].sum() for ibin in range(len(parambins)) ] )
        else :
            parambins = None
            pdf = None


        mean = (weights  * samples[:,ipar]).sum()
        #print(samples[:,ipar]-mean)
        #print(weights)
        std = np.sqrt( (weights * (samples[:,ipar]-mean)**2 ).sum() )


        pdfdict[param] = (parambins,pdf,mean,std,res.logz)

        if verbose :
            if np.abs(std)>=0.1:
                print( '  <%s> =  %.2f +- %.2f'%( param, np.round(mean,2), np.round(std,2))  )
            elif np.abs(std)>=0.01:
                print( '  <%s> =  %.3f +- %.3f'%( param, np.round(mean,3), np.round(std,3)) )
            elif np.abs(std)>=0.001:
                print( '  <%s> =  %.4f +- %.4f'%( param, np.round(mean,4), np.round(std,4)) )
            else :
                print( '  <%s> = %.3e +- %.3e'%( param, mean, std) )


        if param == 'x0' :
            salt2 = sncosmo.Model( source='salt2')
            salt2.source.set_peakmag( 0., 'bessellb', 'ab' )
            x0_AB0 = salt2.get('x0')
            mBmean = -2.5*np.log10(  mean / x0_AB0 )
            mBstd = 2.5*np.log10( np.e ) *  std / mean
            mBbins = -2.5*np.log10(  parambins / x0_AB0 )

            pdfdict['mB'] = ( mBbins, pdf, mBmean, mBstd )
            if verbose:
                print( '  <%s> =  %.3f +- %.3f'%( 'mB', np.round(mBmean,3), np.round(mBstd,3)) )

    return( pdfdict )


def param_fit(args,modName,fit=False):
    sources={'BazinSource':BazinSource}

    source=sources[modName](args['curve'].table,colorCurve=args['color_curve'])
    mod=sncosmo.Model(source)
    if args['constants']:
        mod.set(**args['constants'])
    if not fit:
        res=sncosmo.utils.Result()
        res.vparam_names=args['params']
    else:
        #res,mod=sncosmo.fit_lc(args['curve'].table,mod,args['params'], bounds=args['bounds'],guess_amplitude=True,guess_t0=True,maxcall=1)

        if 'amplitude' in args['bounds']:
            guess_amp_bound=False
        else:
            guess_amp_bound=True




        res,mod=sncosmo.nest_lc(args['curve'].table,mod,vparam_names=args['params'],bounds=args['bounds'],guess_amplitude_bound=guess_amp_bound,maxiter=1000,npoints=200)
    return({'res':res,'model':mod})

