import inspect,sncosmo,os,sys,warnings,pyParz
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from copy import deepcopy
from scipy.interpolate import splrep,splev
from scipy.optimize import curve_fit
from scipy.stats import chisquare

from .io import _get_default_prop_name
from .util import __dir__
from astropy.table import Table

__all__=['fit_data','colorFit']

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
             spline=False, poly=False, micro='spline', **kwargs):
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

    #currently sets the models to run through to all if user doesn't define, will probably change that
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
    #mods = [x for x in sncosmo.models._SOURCES._loaders.keys() if 'snana' in x[0] or 'salt2' in x[0]] if not models else models
    #mods=[mods] if not isinstance(mods,(tuple,list)) else mods
    args['sn_func'] = {'minuit': sncosmo.fit_lc, 'mcmc': sncosmo.mcmc_lc, 'nest': sncosmo.nest_lc}
    #get any properties set in kwargs that exist for the defined fitting function
    args['props'] = {x: kwargs[x] for x in kwargs.keys() if
             x in [y for y in inspect.getargspec(args['sn_func'][method])[0]] and x != 'verbose'}
    #get a unique set of the models to run (not multiple versions of the same model)

    #if not models:
    #    mods = {x[0] if isinstance(x,(tuple,list)) else x for x in mods}
    #elif  not isinstance(models,(tuple,list)):


    for d in curves.images.keys():
        print(d)
        args['curve']=curves.images[d]
        curves.images[d].fits=newDict()
        if len(args['curve'].table)>63:
            fits=[]
            for mod in mods:
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
        curves.images[d].fits=newDict()
        curves.images[d].fits['model']=bestFit
        curves.images[d].fits['res']=bestRes

        #print(curves.images[d].fits.model._source.name)
        #print(curves.images[d].fits.res.param_names,curves.images[d].fits.res.parameters,curves.images[d].fits.res.errors)
        #print(curves.images[d].simMeta)
    #for img in args['curves'].images:
    #        args['curves'].images[img].fits=dict((mod,newDict()))
    #set up parallel processing
    #run each model in parallel and keep track using _pool_results_to_dict

    #print(args['curves'][1].fit['salt2'].res.errors)
    #fig=plt.figure()
    #plt.plot(args['curves'][0].fits.modelFits['salt2']['sdssi'].time,args['curves'][0].fits.modelFits['salt2']['sdssi'].fluxes)
    #plt.show()

    #sncosmo.plot_lc(curves.images['S1'].table,model=curves.images['S1'].fits.model,errors=curves.images['S1'].fits.res.errors)
    #plt.show()

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

    modName=mod+'_'+version if version else deepcopy(mod)
    source=sncosmo.get_source(mod)
    #print(mod)
    smod = sncosmo.Model(source=source)
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
        constants = {x: dcurve.fits.constants[x] for x in fits.constants.keys() if x in smod.param_names}
        smod.set(**constants)

    if args['method']=='mcmc':
        fits.res, fits.model = args['sn_func'][args['method']](dcurve.table, smod, fits.params, fits.bounds, **args['props'])
    elif args['method']=='nest':
        fits.res, fits.model = args['sn_func'][args['method']](dcurve.table, smod, fits.params, fits.bounds,guess_amplitude_bound=True, verbose=False, **args['props'])
    else:
        fits.res, fits.model = args['sn_func'][args['method']](dcurve.table, smod, fits.params, fits.bounds,verbose=False, **args['props'])


    '''
    try:
        dcurve=_snmodel_to_flux(dcurve,modName)
    except RuntimeError:
        print('woops')
        continue
    '''



    #dcurve = _snmodel_to_flux(dcurv e, modName)
    #fits.model._source=None
    #print(fits.res)
    return(pyParz.parReturn(fits))

#def _tdMin(delay,time,curves):
#    return(chisquare())

def _findMax(time,curve):
    #TODO check edge cases
    t0=np.where(curve==np.max(curve))[0][0]

    if t0==0:
        #return(time[0])
        return None

    elif t0==len(time)-1:
        #return(time[-1])
        return None

    else:
        fit=splrep(time[t0-1:t0+2],curve[t0-1:t0+2],k=2)

    interptime=np.linspace(time[t0-1],time[t0+1],100)
    flux=splev(interptime,fit)
    return(interptime[flux==np.max(flux)])


def _findMin(time,curve):
    #TODO check edge cases
    t0=np.where(curve==np.min(curve))[0][0]

    if t0==0:
        #return(time[0])
        return None

    elif t0==len(time)-1:
        #return(time[-1])
        return None

    else:
        fit=splrep(time[t0-1:t0+2],curve[t0-1:t0+2],k=2)

    interptime=np.linspace(time[t0-1],time[t0+1],100)
    flux=splev(interptime,fit)
    return(interptime[flux==np.min(flux)])


def colorFit(lcs):
    #TODO make this work for arbitrary timescales
    colors=combinations(lcs.bands,2)
    figure=plt.figure()
    ax=figure.gca()
    allDelays=[]
    for col in colors:
        curves=dict([])
        for d in lcs.images.keys():
            dcurve=lcs.images[d]
            spl1=splrep(dcurve.table['time'][dcurve.table['band']==col[0]],dcurve.table['flux'][dcurve.table['band']==col[0]])
            spl2=splrep(dcurve.table['time'][dcurve.table['band']==col[1]],dcurve.table['flux'][dcurve.table['band']==col[1]])
            time=np.linspace(max(np.min(dcurve.table['time'][dcurve.table['band']==col[0]]),np.min(dcurve.table['time'][dcurve.table['band']==col[1]])),
                        min(np.max(dcurve.table['time'][dcurve.table['band']==col[0]]),np.max(dcurve.table['time'][dcurve.table['band']==col[1]])),50)
            ccurve=splev(time,spl1)/splev(time,spl2)
            #curves.append(dcurve.table['flux'][dcurve.table['band']==col[0]]-dcurve.table['flux'][dcurve.table['band']==col[1]])
            curves[d]=(time,ccurve)
        ref=False


        delays=dict([])
        for k in curves.keys():
            time,curve=curves[k]
            maxValue=_findMax(time,curve)
            minValue=_findMin(time,curve)
            if not minValue and not maxValue:
                sys.exit()

            if not ref:
                ref=True
                refMax=maxValue
                refMin=minValue
                refName=k
                delays[k]=0

            else:

                #print(maxValue-refMax,minValue-refMin)
                if refMax and maxValue:
                    if refMin and minValue:
                        delays[k]=np.mean([maxValue-refMax,minValue-refMin])
                    else:
                        delays[k]=maxValue-refMax
                elif refMin and minValue:
                    delays[k]=minValue-refMin
                else:
                    sys.exit()



        allDelays.append(delays)
    finalDelays=dict([])
    for k in allDelays[0].keys():

        est=np.mean([x[k] for x in allDelays])
        est=est[0] if isinstance(est,list) else est
        finalDelays[k]=est
        print('True: '+str(lcs.images[k].simMeta['td']-lcs.images[refName].simMeta['td']),'Estimate: '+str(finalDelays[k]))
    for k in curves.keys():
        time,curve=curves[k]
        #print(np.min(time-delays[k]-curves[refName][0]))
        ax.scatter(time,curve)
        #ax.scatter(time-finalDelays[k]-curves[refName][0],curve)
    plt.show()
    return(finalDelays)



            #res=minimize(_tdMin,np.zeros(1),args=(time,curves[i]))





#def maxFit(lcs):
'''
#todo figure out how to deal with negative flux from model flux to mag
def _snmodel_to_flux(dcurve,modName):
    warnings.simplefilter("ignore")
    for band in dcurve.bands:
        if not dcurve.fits.modelFits[modName].model.bandoverlap(band):
            continue
        dcurve.fits.modelFits[modName][band]=sntd.curve(band=band,zp=dcurve[band].zp,zpsys=dcurve[band].zpsys)
        tmin = []
        tmax = []
        tmin.append(np.min(dcurve[band].table[_get_default_prop_name('time')]) - 10)
        tmax.append(np.max(dcurve[band].table[_get_default_prop_name('time')]) + 10)
        tmin.append(dcurve.fits.modelFits[modName].model.mintime())
        tmax.append(dcurve.fits.modelFits[modName].model.maxtime())
        tmin = min(tmin)
        tmax = max(tmax)
        tgrid = np.linspace(tmin, tmax, int(tmax - tmin) + 1)
        mflux = dcurve.fits.modelFits[modName].model.bandflux(band, tgrid, zp=dcurve[band].zp, zpsys=dcurve[band].zpsys)
        mmag = -2.5 * np.log10(mflux) + dcurve[band].zp
        # todo do real error thing

        dcurve.fits.modelFits[modName][band].mags=mmag
        dcurve.fits.modelFits[modName][band].fluxes = mflux
        dcurve.fits.modelFits[modName][band].jds = tgrid

        dcurve.fits.modelFits[modName][band].magerrs = np.abs(mmag)*.1
        dcurve.fits.modelFits[modName][band].fluxerrs = np.abs(mflux)*.1
        #dcurve.fits.modelFits[modName][band]=sntd.io.factory(tgrid,mmag,mflux,band,dcurve[band].zp,dcurve[band].zpsys,np.abs(mmag)*.1,np.abs(mflux)*.1,dcurve[band].telescopename,dcurve[band].object)
    return dcurve
'''