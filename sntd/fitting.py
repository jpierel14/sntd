import inspect,sncosmo,os,sys,warnings,pyParz,math
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy,copy
from scipy import stats
from scipy.interpolate import interp1d
from astropy.table import Table
from astropy.extern import six
import nestle
from collections import OrderedDict

from .util import *
from .io import _sntd_deepcopy,table_factory
from .spl import spl
from .plotting import display
import models
__all__=['fit_data','colorFit','spline_fit']

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


class fits(dict):
    '''
    This is going to be a fits object that will contain a dictionary
    of models with keys as models names. Each value will then be a fitted models
    class from sncosmo.

    If microlensing is added, I guess run through pycs first to get initial time
    delays and microlensing effects, write out the result, and run through sncosmo
    to get best fit models?
    '''
    def __init__(self):
        super(fits, self).__init__() #init for the super class

        self.method='minuit'
        self.params=[]
        self.bounds={}
        self.ignore=[]
        self.modelFits=newDict()
        #self.constants = constants if constants else {x: y for x, y in zip(dcurve.meta.keys(), dcurve.meta.values()) if
        #                                             x != 'info'}
        self.constants={}
        self.spline=False
        self.poly=False
        self.micro='spline'
        #for band in self.bands:
        #    self[band]=sntd.curve()
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





def fit_data(lcs, snType='Ia',bands=None,method='minuit', models=None, params=None, bounds={}, ignore=None, constants=None,
             spline=False, poly=False,micro='spline',combined_or_separate='separate',combinedGrids=None, **kwargs):
    """
    The main, high-level fitting function.
    :param curves: list of objects containing lightcurves to fit, or a single dcurve
    :type curves: ~io.curveDict or list of ~io.curveDict
    :param bands: The list of bands you'd like to fit, optional (all will be fit if you don't specify)
    :type bands: list
    :param method: The fitting method for sncosmo (minuit,mcmc,nest), optional
    :type method: str
    :param models: The models to fit, optional (if not present, all will be tried)
    :type models: str or ~sncosmo.Models
    :param params: The parameters to be varied in the models,optional (if not present, all default parameters used)
    :type params: list
    :param bounds: Bounds for models parameters (necessary for certain parameters)
    :type bounds: dict, {param:(x0,x1)}
    :param ignore: Parameters to ignore in models, optional
    :type ignore: list
    :param constants: Models parameters to be kept constant (i.e., redshift if known), optional
    :type constants: dict, {param:value}
    :param spline: If you'd like to fit a spline to the data instead of an sncosmo template, optional
    :type spline: Boolean
    :param knotstep: Number of knots you'd like in your spline if spline=True
    :type knotstep: int
    :param poly: If you'd like to fit a polynomial to the data instead of an sncosmo template,optional
    :type poly: Boolean
    :param degree: Degree of polynomial if poly=True
    :type degree: int
    :param micro: If don't want to include microlensing, set to None, otherwise choose poly or spline,optional
    :type micro: str (or None)
    :return: fits object containing all fit information for each lightcurve
    """
    #curves=deepcopy(lcs)
    #print(curves.keys())
    curves=_sntd_deepcopy(deepcopy(lcs))

    #curves=deepcopy(lcs)
    args = locals()
    args['curves'] = curves
    args['bands'] = [bands] if bands and not isinstance(bands,(tuple,list)) else bands
    args['t0_guess']=kwargs.get('t0_guess',None)
    args['refModel']=kwargs.get('refModel',None)
    #sets the bands to user's if defined (set, so that they're unique), otherwise to all the bands that exist in curves
    args['bands'] = set(bands) if bands else curves.bands


    args['constants']=constants
    args['effect_names']=kwargs.get('effect_names',[])
    args['effect_frames']=kwargs.get('effect_frames',[])
    args['dust']=kwargs.get('dust',None)
    args['knots']=kwargs.get('knots',3)
    args['func']=kwargs.get('func','spline')
    args['degree']=kwargs.get('degree',3)
    args['smooth']=kwargs.get('smooth',1)
    args['snType']=snType
    args['flip']=kwargs.get('flip',False)
    args['guess_amp']=kwargs.get('guess_amplitude',True)
    args['combinedError']=kwargs.get('combinedError',None)
    args['showPlots']=kwargs.get('showPlots',False)
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
             x in [y for y in inspect.getargspec(args['sn_func'][method])[0]] and x != 'verbose'}

    if combined_or_separate not in ['separate','combined','both']:
        raise RuntimeError('Parameter "combined_or_separate must be "separate","combined", or "both".')
    if combined_or_separate=='separate':
        curves=_fitSeparate(curves,mods,args,bounds)
    elif combined_or_separate=='combined':
        curves=_fitCombined(curves,mods,args,bounds,combinedGrids)
    else:
        curves=_fitSeparate(curves,mods,args,bounds)
        curves=_fitCombined(curves,mods,args,bounds,combinedGrids)
    return curves



def _fitCombined(curves,mods,args,bounds,grids,guess_amplitude_bound=False,
                 minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
                 maxiter=None, maxcall=None, modelcov=False, rstate=None,
                 verbose=False, warn=True, **kwargs):
    """
    if not grids:
        if curves.combined.table:
            args['curve']=curves.combined
        else:
            curves.combine_curves()
            args['curve']=curves.combined
        tdGrid=[[0] for i in range(len(curves.images.keys()))]
        muGrid=[[1] for i in range(len(curves.images.keys()))]
    else:
        if 'td' not in grids:
            tdGrid=[[0] for i in range(len(curves.images.keys()))]
            if 'mu' not in grids:
                raise RuntimeError('Your grid specification does not contain "td" or "mu"')
        else:
            tds=_guess_time_delays(curves)

            tdGrid=np.linspace(grids['td'][0]+tds[k],grids['td'][1]+tds[k],10)
        if 'mu' not in grids:
            muGrid=[[1] for i in range(len(curves.images.keys()))]
        else:
            mus=_guess_magnifications(curves)
            muGrid=np.linspace(grids['td'][0],grids['td'][1],10)
            muGrid=muGrid[muGrid>0]

    """

    if not curves.combined.table:
        curves.combine_curves()
    args['curve']=curves


    for b in [x for x in np.unique(args['curve'].combined.table['band']) if x not in args['bands']]:
        args['curve'].combined.table=args['curve'].combined.table[args['curve'].combined.table['band']!=b]
    #for k in curves.images.keys():
    #    print(curves.images[k].simMeta)
    if 't0' in args['bounds']:
        t0Bounds=copy(args['bounds']['t0'])
    if 'amplitude' in args['bounds']:
        ampBounds=copy(args['bounds']['amplitude'])
        guess_amp_bound=False
    if 't0' in args['bounds'] and args['t0_guess'] is not None:
        args['bounds']['t0']=(t0Bounds[0]+args['t0_guess'],t0Bounds[1]+args['t0_guess'])



    if 'amplitude' in args['bounds'] and args['guess_amp']:
        args['bounds']['amplitude']=(ampBounds[0]*np.max(args['curve'].table['flux']),ampBounds[1]*np.max(args['curve'].table['flux']))
    #if args['refModel']:
        #bestFit,bestRes=args['refModel']
        #refModel=True
    if not args['refModel']:
        refModel=False
        if len(args['curve'].table)>63 or len(models==1) or args['snType']=='Ia':
            fits=[]
            for mod in mods:
                if mod=='SplineSource':

                    fits.append(spline_fit(args))
                elif mod in ['BazinSource','KarpenkaSource','NewlingSource','PierelSource']:
                    fits.append(param_fit(args,mod,fit=True))


                else:
                    fits.append(_fit_data_wrap((mod,args)))
        else:
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

    #if 'amplitude' not in args['bounds']:
    #    guess_amp_bounds=True
    #else:
    #    guess_amp_bounds=False

    #if 'amplitude' in args['bounds'] and args['guess_amp']:
    #    args['bounds']['amplitude']=(ampBounds[0]*np.max(tempTable['flux']),ampBounds[1]*np.max(tempTable['flux']))
    #sncosmo.plot_lc(args['curve'].combined.table,model=bestFit,errors=bestRes.errors)
    #plt.savefig(list(args['bands'])[0]+'_firstCombinedFit.pdf',format='pdf',overwrite=True)
    #plt.close()
    #print(curves.combined.meta)
    gridBounds=dict([])
    for k in curves.images.keys():
        for par in grids.keys():
            if par=='td' and args['curve'].combined.meta[par][k]==0:
                continue
            if par=='mu' and args['curve'].combined.meta[par][k]==1:
                continue
            gridBounds[k+'_'+par]=grids[par]+args['curve'].combined.meta[par][k]
    args['curve'].combined.time_delays,args['curve'].combined.magnifications,args['curve'].combined.time_delay_errors,args['curve'].combined.magnification_errors=nest_combined_lc(args['curve'],
                                        vparam_names=np.append([k+'_td' for k in args['curve'].images.keys() if args['curve'].combined.meta['td'][k]!=0],[k+'_mu' for k in args['curve'].images.keys() if args['curve'].combined.meta['mu'][k]!=1]),
                                        band=list(args['bands'])[0],refModel=args['refModel'],bounds=gridBounds,snBounds=bounds,guess_amplitude_bound=True,maxiter=None,npoints=100)

    args['curve'].combined.time_delay_errors={k:np.sqrt(args['curve'].combined.time_delay_errors[k]**2+args['combinedError']['td']**2) for k in args['curve'].combined.time_delay_errors.keys() if args['curve'].combined.time_delay_errors[k]!=0}
    args['curve'].combined.magnification_errors={k:np.sqrt(args['curve'].combined.magnification_errors[k]**2+args['combinedError']['mu']**2) for k in args['curve'].combined.magnification_errors.keys() if args['curve'].combined.magnification_errors[k] !=1}


    args['curve'].combine_curves(tds=args['curve'].combined.time_delays,mus=args['curve'].combined.magnifications)


    args['curve'].combined.fits=newDict()
    args['curve'].combined.fits['model']=args['refModel']#finalfit

    return args['curve']

def create_composite_model(curves,band,ref,weight='chisq'):
    minTime=np.max([curves.images[im].fits.model.mintime()-curves.images[im].fits.model.get('t0') for im in curves.images.keys()])#+curves.images[ref].fits.model.get('t0')
    maxTime=np.min([curves.images[im].fits.model.maxtime()-curves.images[im].fits.model.get('t0') for im in curves.images.keys()])

    compTime=np.arange(minTime,maxTime,1)
    fluxes=np.array([curves.images[im].fits.model.bandflux(band,compTime+curves.images[im].fits.model.get('t0'),zp=curves.images[im].zp,zpsys=curves.images[im].zpsys)/curves.magnifications[im] for im in np.sort(curves.images.keys())])
    compTime+=curves.images[ref].fits.model.get('t0')

    try:
        weights=np.array([curves.images[im].fits.res[weight] for im in np.sort(curves.images.keys())])
    except:
        raise RuntimeError("You asked to use %s as the weight, but it isn't in your fit."%weight)
    try:
        error={'td':np.average([curves.images[im].fits.final_errs['t0'] for im in curves.images.keys()],weights=1./weights)/np.sqrt(len(weights)),
               'mu':np.average([curves.images[im].fits.final_errs['amplitude'] for im in curves.images.keys()],weights=1./weights)/np.sqrt(len(weights))}
    except:
        error={'td':np.average([curves.images[im].fits.final_errs['t0'] for im in curves.images.keys()],weights=1./weights)/np.sqrt(len(weights)),
               'mu':np.average([curves.images[im].fits.final_errs['A'] for im in curves.images.keys()],weights=1./weights)/np.sqrt(len(weights))}

    finalFlux=np.average(fluxes,weights=1./weights,axis=0)

    compFunc=interp1d(compTime,finalFlux)

    return(compFunc,error)





def nest_combined_lc(curves,vparam_names,bounds,snBounds,guess_amplitude_bound=False,
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

    ref=[x for x in curves.images.keys() if x not in [z[0:2] for z in vparam_names]][0]
    all_delays=dict([])
    all_mus=dict([])
    all_mu_err=dict([])
    all_delay_err=dict([])
    all_mu_err[ref]=0
    all_delay_err[ref]=0
    all_delays[ref]=0
    all_mus[ref]=1
    '''
    for im in [x for x in curves.images.keys() if x !=ref]:
        temp_vparams=[x for x in vparam_names if x[0:2]==im]
        '''
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

    def chisquare(observed_values,expected_values,errors):
        z = (observed_values - expected_values) / errors
        chi2 = np.sum(z ** 2)

        return chi2


    def loglike(parameters):
        tempCurve=_sntd_deepcopy(deepcopy(curves))

        tempTds=dict([])
        tempMus=dict([])
        for i in range(len(parameters)):
            if iparam_names[i][-2:]=='td':
                tempTds[iparam_names[i][0:2]]=parameters[i]
            elif iparam_names[i][-2:]=='mu':
                tempMus[iparam_names[i][0:2]]=parameters[i]
        #for k in all_delays:
        tempTds[ref]=0#all_delays[k]
        tempMus[ref]=1#all_mus[k]
        #tempCurve.images={k:v for k,v in tempCurve.images.items() if k in tempTds}

                          #tempTds[ref]=0
        #tempMus[ref]=1
        tempCurve.combine_curves(tds=tempTds,mus=tempMus)


        #tempRes,tempMod=sncosmo.nest_lc(tempCurve.combined.table,model,vparam_names=res.vparam_names,bounds=snBounds,guess_amplitude_bound=True,maxiter=1000,npoints=50)
        tempCurve.combined.table=tempCurve.combined.table[tempCurve.combined.table['band']==band]
        tempCurve.combined.table.sort('time')
        chisq=-chisquare(tempCurve.combined.table['flux'],refModel(tempCurve.combined.table['time']),tempCurve.combined.table['fluxerr'])/(len(tempCurve.combined.table)-1)

        return(chisq)




    #t0 = time.time()
    res = nestle.sample(loglike, prior_transform, ndim, npdim=npdim,
                        npoints=npoints, method=method, maxiter=maxiter,
                        maxcall=maxcall, rstate=rstate,
                        callback=(nestle.print_progress if verbose else None))

    #elapsed = time.time() - t0
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



    return all_delays,all_mus,all_delay_err,all_mu_err#(ref,dict(zip(vparam_names,vparameters)),OrderedDict(zip(vparam_names,
            #                                                        np.sqrt(np.diagonal(cov)))))

# def nest_combined_lc(curves,model,res,vparam_names,bounds,snBounds,guess_amplitude_bound=False,
#                      minsnr=5.,refModel=False, priors=None, ppfs=None, npoints=100, method='single',
#                      maxiter=None, maxcall=None, modelcov=False, rstate=None,
#                      verbose=False, warn=True, **kwargs):
#
#     # experimental parameters
#     tied = kwargs.get("tied", None)
#
#
#     model=copy(model)
#     if ppfs is None:
#         ppfs = {}
#     if tied is None:
#         tied = {}
#
#     # Convert bounds/priors combinations into ppfs
#     if bounds is not None:
#         for key, val in six.iteritems(bounds):
#             if key in ppfs:
#                 continue  # ppfs take priority over bounds/priors
#             a, b = val
#             if priors is not None and key in priors:
#                 # solve ppf at discrete points and return interpolating
#                 # function
#                 x_samples = np.linspace(0., 1., 101)
#                 ppf_samples = sncosmo.utils.ppf(priors[key], x_samples, a, b)
#                 f = sncosmo.utils.Interp1D(0., 1., ppf_samples)
#             else:
#                 f = sncosmo.utils.Interp1D(0., 1., np.array([a, b]))
#             ppfs[key] = f
#
#     # NOTE: It is important that iparam_names is in the same order
#     # every time, otherwise results will not be reproducible, even
#     # with same random seed.  This is because iparam_names[i] is
#     # matched to u[i] below and u will be in a reproducible order,
#     # so iparam_names must also be.
#     iparam_names = [key for key in vparam_names if key in ppfs]
#     ppflist = [ppfs[key] for key in iparam_names]
#     npdim = len(iparam_names)  # length of u
#     ndim = len(vparam_names)  # length of v
#
#     # Check that all param_names either have a direct prior or are tied.
#     for name in vparam_names:
#         if name in iparam_names:
#             continue
#         if name in tied:
#             continue
#         raise ValueError("Must supply ppf or bounds or tied for parameter '{}'"
#                          .format(name))
#
#     def prior_transform(u):
#         d = {}
#         for i in range(npdim):
#             d[iparam_names[i]] = ppflist[i](u[i])
#         v = np.empty(ndim, dtype=np.float)
#         for i in range(ndim):
#             key = vparam_names[i]
#             if key in d:
#                 v[i] = d[key]
#             else:
#                 v[i] = tied[key](d)
#         return v
#
#     snParams=[]
#     ref=[x for x in curves.images.keys() if x not in [z[0:2] for z in iparam_names]][0]
#
#     def loglike(parameters):
#         tempCurve=_sntd_deepcopy(deepcopy(curves))
#         tempTds=dict([])
#         tempMus=dict([])
#         for i in range(len(parameters)):
#             if iparam_names[i][-2:]=='td':
#                 tempTds[iparam_names[i][0:2]]=parameters[i]
#             elif iparam_names[i][-2:]=='mu':
#                 tempMus[iparam_names[i][0:2]]=parameters[i]
#         tempTds[ref]=0
#         tempMus[ref]=1
#         tempCurve.combine_curves(tds=tempTds,mus=tempMus)
#         #tempRes,tempMod=sncosmo.nest_lc(tempCurve.combined.table,model,vparam_names=res.vparam_names,bounds=snBounds,guess_amplitude_bound=True,maxiter=1000,npoints=50)
#         if refModel:
#             for b in [x for x in np.unique(tempCurve.combined.table['band']) if not model.bandoverlap(x)]:
#                 tempCurve.combined.table=tempCurve.combined.table[tempCurve.combined.table['band']!=b]
#             return(-0.5*sncosmo.fitting.chisq(tempCurve.combined.table,model))
#
#         else:
#             tempRes,tempMod=sncosmo.fit_lc(tempCurve.combined.table,model,vparam_names=res.vparam_names,bounds=snBounds,maxcall=5)
#             #snParams.append((parameters,tempRes,tempMod))
#             #snParams.append(parameters)
#             #return tempRes.logz
#             return -0.5*tempRes.chisq
#
#
#     #t0 = time.time()
#     res = nestle.sample(loglike, prior_transform, ndim, npdim=npdim,
#                         npoints=npoints, method=method, maxiter=maxiter,
#                         maxcall=maxcall, rstate=rstate,
#                         callback=(nestle.print_progress if verbose else None))
#     #elapsed = time.time() - t0
#
#     # estimate parameters and covariance from samples
#     vparameters, cov = nestle.mean_and_cov(res.samples, res.weights)
#
#
#     return (ref,dict(zip(vparam_names,vparameters)),OrderedDict(zip(vparam_names,
#                                                                     np.sqrt(np.diagonal(cov)))))




def _fitSeparate(curves,mods,args,bounds):
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

        if 'amplitude' in args['bounds'] and args['guess_amp']:
            args['bounds']['amplitude']=(ampBounds[0]*np.max(args['curve'].table['flux']),ampBounds[1]*np.max(args['curve'].table['flux']))



        curves.images[d].fits=newDict()
        if len(args['curve'].table)>63 or len(mods)==1 or curves['snType']=='Ia':
            fits=[]
            for mod in mods:
                if mod=='SplineSource':

                    fits.append(spline_fit(args))
                elif mod in ['BazinSource','KarpenkaSource','NewlingSource','PierelSource']:
                    fits.append(param_fit(args,mod))


                else:
                    if len(mods)==1:
                        doFit=False
                    else:
                        doFit=True
                    args['doFit']=doFit
                    fits.append(_fit_data_wrap((mod,args)))
        else:
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
            #print(bestFit.get('amplitude'))
            #sncosmo.plot_lc(args['curve'].table,model=bestFit,errors=bestRes.errors,zp=30.,zpsys='ab')
            #plt.show()
            #plt.close()
        #for param in [p for p in __thetaSN__ if p in bestRes.vparam_names]:
        #    if param not in newBounds:
        #        newBounds[param]=(.75*bestFit.get(param),1.25*bestFit.get(param))
        #    else:
        #        lower=np.mean([newBounds[param][0],.75*bestFit.get(param)])
        #        upper=np.mean([newBounds[param][1],1.25*bestFit.get(param)])
        #        newBounds[param]=(lower,upper)

        fitDict[d]=[fits,bestFit,bestRes]
    #for param in [p for p in bounds if p not in newBounds.keys()]:
    #    newBounds[param]=bounds[param]
    #if all the best models aren't the same, take the one with min chisq (not the best way to do this)
    if not all([fitDict[d][1]._source.name==fitDict[fitDict.keys()[0]][1]._source.name for d in fitDict.keys()]):
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
        tempTable=copy(curves.images[d].table)
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
        if 'amplitude' in args['bounds'] and args['guess_amp']:
            args['bounds']['amplitude']=(ampBounds[0]*np.max(tempTable['flux']),ampBounds[1]*np.max(tempTable['flux']))
        if 'amplitude' not in args['bounds']:
            guess_amp_bounds=True
        else:
            guess_amp_bounds=False
        nest_res,nest_fit=sncosmo.nest_lc(tempTable,bestFit,vparam_names=bestRes.vparam_names,bounds=args['bounds'],guess_amplitude_bound=guess_amp_bounds,maxiter=1000,npoints=200)
        #print(d,nest_res.h)

        if nest_res.ndof != len(tempTable)- len(nest_res.vparam_names):
            dofs[d]=len(tempTable)- len(nest_res.vparam_names)
        else:
            dofs[d]=nest_res.ndof

        #tempMod=copy(nest_fit)
        #tempMod._source._phase=tempMod._source._phase[tempMod._source._phase>=(np.min(tempTable['time'])-tempMod.get('t0'))]
        #tempMod._source._phase=tempMod._source._phase[tempMod._source._phase<=(np.max(tempTable['time'])-tempMod.get('t0'))]
        #sncosmo.plot_lc(tempTable,model=tempMod,errors=nest_res.errors,zp=tempTable['zp'][0],zpsys=tempTable['zpsys'][0])
        #plt.show()
        #plt.close()

        resList[d]=nest_res
        curves.images[d].fits=newDict()
        curves.images[d].fits['model']=nest_fit
        curves.images[d].fits['res']=nest_res
        #errs[d]=nest_res.errors
        #print(curves.images[d].simMeta)


    joint=_joint_likelihood(resList,verbose=False)
    #for p in joint.keys():

    for d in np.sort(curves.images.keys()):
        errs=dict([])
        tempTable=copy(curves.images[d].table)

        for b in [x for x in np.unique(tempTable['band']) if x not in args['bands']]:
            tempTable=tempTable[tempTable['band']!=b]
        if not np.all([x in __thetaL__ for x in joint.keys()]):
            bds=dict([])

            final_vparams=[]
            for p in joint.keys():
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


                    #bds[p]=(joint[p][d][0]-joint[p][d][1],joint[p][d][0]+joint[p][d][1])
                    #errs[p]=joint[p][d][1]
                else:
                    curves.images[d].fits.model.set(**{p:joint[p][0]})
                    errs[p]=joint[p][1]
                    #final_vparams.append(p)
            #finalRes,finalFit=sncosmo.fit_lc(tempTable,curves.images[d].fits.model,final_vparams,bounds=bds,guess_amplitude=False,guess_t0=False,maxcall=500)
            finalRes,finalFit=sncosmo.nest_lc(tempTable,curves.images[d].fits.model,final_vparams,bounds=bds,guess_amplitude_bound=True,maxiter=500)

            finalRes.ndof=dofs[d]
            #print(d,finalRes.chisq/finalRes.ndof,stats.chi2.sf(finalRes.chisq,finalRes.ndof),finalRes.chisq,finalRes.ndof)
            curves.images[d].fits=newDict()
            curves.images[d].fits['model']=finalFit#nest_fit
            curves.images[d].fits['res']=finalRes

        else:
            for p in joint.keys():
                if isinstance(joint[p],dict):
                    curves.images[d].fits.model.set(**{p:joint[p][d][0]})
                    errs[p]=joint[p][d][1]

                else:
                    curves.images[d].fits.model.set(**{p:joint[p][0]})
                    errs[p]=joint[p][1]
        curves.images[d].fits['final_errs']=errs

    #rint(errs)
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

            tempMod._source._phase=tempMod._source._phase[tempMod._source._phase>=(np.min(tempTable['time'])-tempMod.get('t0'))]
            tempMod._source._phase=tempMod._source._phase[tempMod._source._phase<=(np.max(tempTable['time'])-tempMod.get('t0'))]
            sncosmo.plot_lc(tempTable,model=tempMod,zp=tempTable['zp'][0],zpsys=tempTable['zpsys'][0])
            #plt.savefig(nest_fit._source.name+'_'+tempTable['band'][0]+'_refs_'+d+'.pdf',format='pdf',overwrik4ite=True)
            plt.show()
            plt.clf()
            plt.close()

    return curves

def _maxFromModel(mod,band,zp,zpsys):
      time=np.arange(mod.mintime(),mod.maxtime(),.1)
      flux=mod.bandflux(band,time,zp,zpsys)
      return (time[flux==np.max(flux)],np.max(flux))

def timeDelaysAndMagnifications(curves):
    times=dict([])
    fluxes=dict([])
    time_errors=dict([])
    flux_errors=dict([])
    for d in curves.images.keys():
        band_times=dict([])
        band_fluxes=dict([])
        band_time_errors=dict([])
        band_flux_errors=dict([])
        for b in [x for x in curves.bands if curves.images[d].fits.model.bandoverlap(x)]:
            try:
                timeOfPeak=curves.images[d].fits.model.get('t0')
                band_time_errors[b]=curves.images[d].fits.final_errs['t0']

            except:
                timeOfPeak,peakFlux=_maxFromModel(curves.images[d].fits.model,b,curves.images[d].zp,curves.images[d].zpsys)
                band_time_errors[b]=0
                band_flux_errors[b]=0


            try:
                peakFlux=curves.images[d].fits.model.get('amplitude')
                band_flux_errors[b]=curves.images[d].fits.final_errs['amplitude']
            except:
                try:
                    peakFlux=curves.images[d].fits.model.get('A')
                    band_flux_errors[b]=curves.images[d].fits.final_errs['A']
                except:
                    peakFlux=curves.images[d].fits.model.bandflux(b,timeOfPeak,zp=curves.images[d].zp,zpsys=curves.images[d].zpsys)
                    band_flux_errors[b]=0



            band_times[b]=timeOfPeak
            band_fluxes[b]=peakFlux
        times[d]=band_times
        fluxes[d]=band_fluxes
        time_errors[d]=band_time_errors
        flux_errors[d]=band_flux_errors
    ims=np.sort(curves.images.keys())
    delays=dict([])
    mags=dict([])
    delay_errs=dict([])
    mag_errs=dict([])
    ref1=None
    ref2=None
    ref1_err=None
    ref2_err=None
    for b in [x for x in curves.bands if curves.images[d].fits.model.bandoverlap(x)]:
        for im in ims:
            if not ref1:
                ref1=times[im][b]
                delays[im]=0
                ref2=fluxes[im][b]
                mags[im]=1
                ref1_err=time_errors[im][b]
                ref2_err=flux_errors[im][b]
                delay_errs[im]=0
                mag_errs[im]=0

            else:
                delays[im]=ref1-times[im][b]
                mags[im]=fluxes[im][b]/ref2
                delay_errs[im]=math.sqrt(time_errors[im][b]**2+ref1_err**2)
                mag_errs[im]=mags[im]*math.sqrt((flux_errors[im][b]/fluxes[im][b])**2+(ref2_err/ref2)**2)
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

    source=sncosmo.get_source(mod)

    #print(mod)
    smod = sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
    params=args['params'] if args['params'] else [x for x in smod.param_names]

    fits=newDict()
    dcurve=args['curve']
    if not np.any([smod.bandoverlap(band) for band in dcurve.bands]):
        print('yep')
        raise RuntimeError("No band overlap for model %s"%modName)
    fits.method=args['method']
    fits.bounds=args['bounds'] if args['bounds'] else {}
    fits.ignore=args['ignore'] if args['ignore'] else []
    fits.constants = args['constants'] if args['constants'] else {x: y for x, y in zip(dcurve.meta.keys(),dcurve.meta.values()) if x != 'info'}
    fits.spline = args['spline']
    fits.poly = args['poly']
    fits.micro = args['micro']

    no_bound = {x for x in params if x in _needs_bounds and x not in fits.bounds.keys() and x not in fits.constants.keys()}
    if no_bound:
        params=list(set(params)-no_bound)
    params= [x for x in params if x not in fits.ignore and x not in fits.constants.keys()]
    fits.params = params
    if fits.constants:
        constants = {x: fits.constants[x] for x in fits.constants.keys() if x in smod.param_names}
        smod.set(**constants)

    if args['doFit']:
        if args['method']=='mcmc':
            fits.res, fits.model = args['sn_func'][args['method']](dcurve.table, smod, fits.params, fits.bounds, **args['props'])
        elif args['method']=='nest':
            fits.res, fits.model = args['sn_func'][args['method']](dcurve.table, smod, fits.params, fits.bounds,guess_amplitude_bound=True, verbose=False, **args['props'])
        else:
            fits.res, fits.model = args['sn_func'][args['method']](dcurve.table, smod, fits.params, fits.bounds,verbose=False, **args['props'])
        return(pyParz.parReturn(fits))
    else:

        fits.model=smod
        fits.res=newDict()
        fits.res['vparam_names']=fits.params
        return (fits)




#def _tdMin(delay,time,curves):
#    return(chisquare())

def _joint_likelihood(resList,verbose=False):


    params=[]
    for res in resList.values():
        params=list(set(np.append(params,res.vparam_names)))
    otherParams=[x for x in params if x in __thetaL__]
    snparams=[x for x in params if x in __thetaSN__]
    outDict={p:dict([]) for p in otherParams}
    for param in params:
        if verbose:
            print(param)
        probsList=[]
        weightsList=[]
        zweights=[]
        testList=[]
        for k in resList.keys():
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
            parambins = np.linspace( parvalmin, parvalmax, nbins, endpoint=True )
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
    sources={'NewlingSource':models.NewlingSource,'KarpenkaSource':models.KarpenkaSource,'BazinSource':models.BazinSource,'PierelSource':models.PierelSource}

    source=sources[modName](args['curve'].table)
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
        #return
    #t0,_=sncosmo.fitting.guess_t0_and_amplitude(sncosmo.photdata.photometric_data(args['curve'].table),mod,5.0)
    #if 't0' in args['params']:
    #    args['bounds']['t0']=(t0-10,t0+10)
    #if 'phi' in args['params']:# and 'phi' not in args['bounds'].keys()
    #    args['bounds']['phi']=(np.min(args['curve'].table['time'])-100-t0,np.min(args['curve'].table['time'])-t0)
    #    print(args['bounds']['phi'])
    #print(fit.get('t0'))
    #if modName=='NewlingSource':
        #fit.set(phi=fit.get('t0')-fit.get('k')*fit.get('sigma'))
    #    fit.set(t0=fit.get('t0'))
    #    res.vparam_names=[x for x in res.vparam_names if x != 't0']
    #print(fit.get('t0'),fit.get('fall'),fit.get('rise'),res.chisq/res.ndof)
    #sncosmo.plot_lc(args['curve'].table,model=fit,errors=res)
    #plt.show()
    #plt.clf()
    #sys.exit()
    return({'res':res,'model':mod})

def spline_fit(args):
    source=models.SplineSource(args['curve'].table,knots=args['knots'],smooth=args['smooth'],func=args['func'],degree=args['degree'])
    mod=sncosmo.Model(source)

    if args['constants']:
        mod.set(**args['constants'])
    #mod.set(t0=0)
    res,fit=sncosmo.fit_lc(args['curve'].table,mod,args['params'],bounds=args['bounds'],guess_amplitude=False,guess_t0=False,maxcall=1)
    #sncosmo.plot_lc(args['curve'].table,model=fit,errors=res)
    #plt.show()
    return({'res':res,'model':fit})

'''
def spline_fit(curves):
    shifts=_guess_time_delays(curves)
    lcs=curves.images.values()
    keys=curves.images.keys()
    if 'mag' not in curves.table.colnames:
        for lc in lcs:
            bandDict={x:sncosmo.get_bandpass(x) for x in lc.bands}
            lc.table=flux_to_mag(lc.table,bandDict,zpsys=lc.zpsys)


    for b in curves.bands:
        toFit=[]
        fig=plt.figure()
        ax=fig.gca()
        for i in range(len(lcs)):
            if b in lcs[i].bands:
                tempCurve=copy(lcs[i])
                tempCurve.table=tempCurve.table[tempCurve.table['band']==b]
                tempCurve.fluxes=tempCurve.table['flux']
                tempCurve.fluxerrs=tempCurve.table['fluxerr']
                tempCurve.mags=tempCurve.table['mag']
                tempCurve.magerrs=tempCurve.table['magerr']
                tempCurve.jds=tempCurve.table['time']

                ax.scatter(tempCurve.jds,tempCurve.mags)

                tempCurve.mask=np.array([True for j in range(len(tempCurve.mags))])
                tempCurve.shifttime(shifts[keys[i]])
                tempCurve.object=keys[i]
                pycs.gen.polyml.addtolc(tempCurve,nparams=3)
                toFit.append(tempCurve)
        ax.invert_yaxis()
        plt.savefig(b+'mags.pdf',format='pdf',overwrite=True)
        plt.clf()
        spline=spl(toFit)

        #delay = pycs.gen.lc.gettimeshifts(toFit)[0]


        print(pycs.gen.lc.getnicetimedelays(toFit, separator="\n"))
        display(toFit, [spline],knotsize=0.01,style='posterpdf',filename=b+'_spline.pdf')
    for k in curves.images.keys():
        print(k+': ')
        print(curves.images[k].simMeta)
    sys.exit()
'''

