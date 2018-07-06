import inspect,sncosmo,os,sys,warnings,pyParz,pycs
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy,copy
from scipy.interpolate import splrep,splev
from astropy.table import Table
from astropy.extern import six
import nestle

from .util import *
from .spl import spl
from .plotting import display
import models
__all__=['fit_data','colorFit','spline_fit']

__thetaSN__=['z','hostebv','screenebv','screenz','rise','fall','sigma','k']
__thetaL__=['t0','amplitude','dt0','A','B','t1','psi','phi']


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





def fit_data(curves, snType='Ia',bands=None,method='minuit', models=None, params=None, bounds=None, ignore=None, constants=None,
             spline=False, poly=False, micro='spline',combined_or_separate='separate',combinedGrids=None, **kwargs):
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
    args = locals()
    args['curves'] = curves
    args['bands'] = [bands] if bands and not isinstance(bands,(tuple,list)) else bands
    #sets the bands to user's if defined (set, so that they're unique), otherwise to all the bands that exist in curves
    args['bands'] = set(bands) if bands else curves.bands
    args['constants']=constants
    args['effect_names']=kwargs.get('effect_names',[])
    args['effect_frames']=kwargs.get('effect_frames',[])
    args['dust']=kwargs.get('dust',None)
    args['knots']=kwargs.get('knots',3)
    args['func']=kwargs.get('func','spline')
    args['degree']=kwargs.get('degree',3)
    args['snType']=snType
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
    if curves.combined.table:
        args['curve']=curves.combined
    else:
        curves.combine_curves()
        args['curve']=curves.combined

    for k in curves.images.keys():
        print(curves.images[k].simMeta)

    if len(args['curve'].table)>63:
        fits=[]
        for mod in mods:
            fits.append(_fit_data_wrap((mod,args)))
    else:
        fits=pyParz.foreach(mods,_fit_data,args)

    bestChisq=np.inf
    bestFit=None
    bestRes=None
    for f in fits:
        if f:
            res=f['res']
            mod=f['model']
            if res.chisq <bestChisq:
                bestChisq=res.chisq
                bestFit=mod
                bestRes=res
    sncosmo.plot_lc(curves.combined.table,model=bestFit,errors=bestRes)
    plt.savefig('firstCombinedFit.pdf',format='pdf',overwrite=True)
    plt.close()
    print(curves.combined.meta)
    #TODO don't do this so simply, include in nested sampling or something (scipy optimize next)
    gridBounds=dict([])
    for k in curves.images.keys():
        for par in grids.keys():
            if par=='td' and curves.combined.meta[par][k]==0:
                continue
            if par=='mu' and curves.combined.meta[par][k]==1:
                continue
            gridBounds[k+'_'+par]=grids[par]+curves.combined.meta[par][k]

    ref,names,best,res,params=nest_combined_lc(curves,bestFit,bestRes,vparam_names=np.append([k+'_td' for k in curves.images.keys() if curves.combined.meta['td'][k]!=0],[k+'_mu' for k in curves.images.keys() if curves.combined.meta['mu'][k]!=1]),bounds=gridBounds,snBounds=bounds,guess_amplitude_bound=True,maxiter=100,npoints=10)
    print(res)
    best_nest_fit=None
    best_nest_res=None

    for i in range(len(params)):
        p,nest_res,nest_fit=params[i]
        #testRes=dict(zip(names,samples[i]))
        if np.all(best==p):
            best_nest_fit=nest_fit
            best_nest_res=nest_res
            break
    if not best_nest_fit:
        raise RuntimeError('Huh, no matches from double nestle!')
    mus={k[0:2]:res[k] for k in res.keys() if k[-2:]=='mu'}
    tds={k[0:2]:res[k] for k in res.keys() if k[-2:]=='td'}
    mus[ref]=1
    tds[ref]=0
    curves.combine_curves(tds=tds,mus=mus)
    sncosmo.plot_lc(curves.combined.table,model=best_nest_fit,errors=best_nest_res)
    plt.savefig('nestedCombinedFit.pdf',format='pdf',overwrite=True)
    plt.close()
    print(curves.combined.meta)
    sys.exit()

    curves.combined.fit=newDict()
    curves.combined.fit['model']=bestFit
    curves.combined.fit['res']=bestRes
    return curves



def nest_combined_lc(curves,model,res,vparam_names,bounds,snBounds,guess_amplitude_bound=False,
                     minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
                     maxiter=None, maxcall=None, modelcov=False, rstate=None,
                     verbose=False, warn=True, **kwargs):

    # experimental parameters
    tied = kwargs.get("tied", None)


    model=copy(model)
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

    snParams=[]
    ref=[x for x in ['S'+str(i+1) for i in range(len(curves.images.keys()))] if x not in [z[0:2] for z in iparam_names]][0]
    def loglike(parameters):
        tempCurve=copy(curves)
        tempTds=dict([])
        tempMus=dict([])
        for i in range(len(parameters)):
            if iparam_names[i][-2:]=='td':
                tempTds[iparam_names[i][0:2]]=parameters[i]
            elif iparam_names[i][-2:]=='mu':
                tempMus[iparam_names[i][0:2]]=parameters[i]
        tempTds[ref]=0
        tempMus[ref]=1
        tempCurve.combine_curves(tds=tempTds,mus=tempMus)
        #tempRes,tempMod=sncosmo.nest_lc(tempCurve.combined.table,model,vparam_names=res.vparam_names,bounds=snBounds,guess_amplitude_bound=True,maxiter=1000,npoints=50)
        tempRes,tempMod=sncosmo.fit_lc(tempCurve.combined.table,model,vparam_names=res.vparam_names,bounds=snBounds,maxcall=50)
        snParams.append((parameters,tempRes,tempMod))
        #return tempRes.logz
        return -0.5*tempRes.chisq


    #t0 = time.time()
    res = nestle.sample(loglike, prior_transform, ndim, npdim=npdim,
                        npoints=npoints, method=method, maxiter=maxiter,
                        maxcall=maxcall, rstate=rstate,
                        callback=(nestle.print_progress if verbose else None))
    #elapsed = time.time() - t0

    # estimate parameters and covariance from samples
    vparameters, cov = nestle.mean_and_cov(res.samples, res.weights)


    return (ref,iparam_names,vparameters,dict(zip(vparam_names,vparameters)),snParams)



def _fitSeparate(curves,mods,args,bounds):
    resList=dict([])
    fitDict=dict([])
    #newBounds=dict([])
    for d in curves.images.keys():
        #print(curves.images[d].simMeta)
        args['curve']=curves.images[d]
        curves.images[d].fits=newDict()
        if len(args['curve'].table)>63 or len(mods)==1 or curves['snType']=='Ia':
            fits=[]
            for mod in mods:
                if mod=='SplineSource':

                    fits.append(spline_fit(args))
                elif mod in ['BazinSource','KarpenkaSource','NewlingSource']:
                    fit,bounds=param_fit(args,mod)
                    fits.append(fit)

                else:
                    fits.append(_fit_data_wrap((mod,args)))
        else:
            fits=pyParz.foreach(mods,_fit_data,args)

        bestChisq=np.inf
        for f in fits:
            if f:
                res=f['res']
                mod=f['model']
                if res.chisq <bestChisq:
                    bestChisq=res.chisq
                    bestFit=mod
                    bestRes=res
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

    for d in fitDict.keys():
        _,bestFit,bestMod=fitDict[d]
        nest_res,nest_fit=sncosmo.nest_lc(curves.images[d].table,bestFit,vparam_names=bestRes.vparam_names,bounds=bounds,guess_amplitude_bound=False,maxiter=50,npoints=10)
        #sncosmo.plot_lc(data=curves.images[d].table,model=nest_fit,errors=nest_res.errors)
        #plt.show()
        #plt.close()
        resList[d]=nest_res
        curves.images[d].fits=newDict()
        curves.images[d].fits['model']=nest_fit
        curves.images[d].fits['res']=nest_res
        #print(curves.images[d].simMeta)


    joint=_joint_likelihood(resList,verbose=True)
    for p in joint.keys():
        for d in curves.images.keys():
            if isinstance(joint[p],dict):
                curves.images[d].fits.model.set(**{p:joint[p][d][0]})
            else:
                curves.images[d].fits.model.set(**{p:joint[p][0]})
            #print(curves.images[d].fits.model._source.t0)
    for d in curves.images.keys():
        sncosmo.plot_lc(curves.images[d].table,model=curves.images[d].fits.model,errors=curves.images[d].fits.res)
        plt.show()
        plt.close()

    return curves

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

    if args['method']=='mcmc':
        fits.res, fits.model = args['sn_func'][args['method']](dcurve.table, smod, fits.params, fits.bounds, **args['props'])
    elif args['method']=='nest':
        fits.res, fits.model = args['sn_func'][args['method']](dcurve.table, smod, fits.params, fits.bounds,guess_amplitude_bound=True, verbose=False, **args['props'])
    else:
        fits.res, fits.model = args['sn_func'][args['method']](dcurve.table, smod, fits.params, fits.bounds,verbose=False, **args['props'])


    return(pyParz.parReturn(fits))

#def _tdMin(delay,time,curves):
#    return(chisquare())

def _joint_likelihood(resList,verbose=True):


    params=[]
    for res in resList.values():
        params=list(set(np.append(params,res.vparam_names)))
    otherParams=[x for x in params if x in __thetaL__]
    snparams=[x for x in params if x in __thetaSN__]
    outDict={p:dict([]) for p in otherParams}
    for param in params:
        print(param)
        probsList=[]
        weightsList=[]
        for k in resList.keys():
            res=resList[k]
            pdf=_get_marginal_pdfs(res,nbins=100,verbose=False)
            if param in snparams:
                weightsList=np.append(weightsList,pdf[param][1])
                probsList=np.append(probsList,pdf[param][0])
            elif param in otherParams:
                if verbose:
                    print( 'Image %s:  <%s> = %.6e +- %.2e'%(k, param, pdf[param][2], pdf[param][3]) )
                outDict[param][k]=(pdf[param][2],pdf[param][3])
        if param in otherParams:
            continue
        expectedValue=(weightsList*probsList).sum()
        expectedValue/=len(resList)
        std = np.sqrt( ((weightsList * (probsList-expectedValue)**2 ).sum())/4 )
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
    import numpy as np
    import sncosmo
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
        std = np.sqrt( (weights * (samples[:,ipar]-mean)**2 ).sum() )

        pdfdict[param] = (parambins,pdf,mean,std)

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
            print( '  <%s> =  %.3f +- %.3f'%( 'mB', np.round(mBmean,3), np.round(mBstd,3)) )

    return( pdfdict )


def param_fit(args,modName):
    sources={'NewlingSource':models.NewlingSource,'KarpenkaSource':models.KarpenkaSource,'BazinSource':models.BazinSource}
    source=sources[modName](args['curve'].table)
    mod=sncosmo.Model(source)
    if args['constants']:
        mod.set(**args['constants'])
    t0,_=sncosmo.fitting.guess_t0_and_amplitude(sncosmo.photdata.photometric_data(args['curve'].table),mod,5.0)
    if 't0' in args['params']:
        args['bounds']['t0']=(t0-10,t0+10)
    if 'phi' in args['params']:# and 'phi' not in args['bounds'].keys()
        args['bounds']['phi']=(np.min(args['curve'].table['time'])-100-t0,np.min(args['curve'].table['time'])-t0)
        print(args['bounds']['phi'])
    res,fit=sncosmo.fit_lc(args['curve'].table,mod,args['params'], bounds=args['bounds'],guess_amplitude=False,guess_t0=True,maxcall=1)
    #print(fit.get('t0'))
    #if modName=='NewlingSource':
        #fit.set(phi=fit.get('t0')-fit.get('k')*fit.get('sigma'))
    #    fit.set(t0=fit.get('t0'))
    #    res.vparam_names=[x for x in res.vparam_names if x != 't0']
    #print(fit.get('t0'),fit.get('fall'),fit.get('rise'),res.chisq/res.ndof)
    #sncosmo.plot_lc(args['curve'].table,model=fit,errors=res)
    #plt.show()
    #sys.exit()
    return({'res':res,'model':fit},args['bounds'])

def spline_fit(args):
    source=models.SplineSource(args['curve'].table,knots=args['knots'],func=args['func'],degree=args['degree'])
    mod=sncosmo.Model(source)

    if args['constants']:
        mod.set(**args['constants'])
    mod.set(t0=0)
    res,fit=sncosmo.fit_lc(args['curve'].table,mod,['dt0','amplitude'],bounds=args['bounds'],guess_amplitude=False,guess_t0=False,maxcall=20)
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

