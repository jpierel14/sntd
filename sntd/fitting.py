import inspect
import sncosmo
import sntd
import os,sys
import warnings
import numpy as np
import functools
import multiprocessing
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.stats import norm
from copy import deepcopy
from .io import _get_default_prop_name
from astropy.table import Table

__all__=['fit_data']

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

def fit_data(curves, bands=None,method='minuit', models=None, params=None, bounds=None, ignore=None, constants=None,
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
    args['curves'] = [curves] if not isinstance(curves, (tuple, list)) else curves
    args['bands'] = [bands] if bands and not isinstance(bands,(tuple,list)) else bands
    #sets the bands to user's if defined (set, so that they're unique), otherwise to all the bands that exist in curves
    args['bands'] = set(bands) if bands else set(
        np.hstack(np.array([x for x in [list(y.bands) for y in args['curves']]])))
    if not args['bands']:
        raise (RuntimeError, "You don't have any bands to analyze!")
    #currently sets the models to run through to all if user doesn't define, will probably change that
    mods = [x for x in sncosmo.models._SOURCES._loaders.keys() if 'snana' in x[0] or 'salt2' in x[0]] if not models else models
    #mods=[mods] if not isinstance(mods,(tuple,list)) else mods
    args['sn_func'] = {'minuit': sncosmo.fit_lc, 'mcmc': sncosmo.mcmc_lc, 'nest': sncosmo.nest_lc}
    #get any properties set in kwargs that exist for the defined fitting function
    args['props'] = {x: kwargs[x] for x in kwargs.keys() if
             x in [y for y in inspect.getargspec(args['sn_func'][method])[0]] and x != 'verbose'}
    #get a unique set of the models to run (not multiple versions of the same model)

    if not models:
        mods = {x[0] if isinstance(x,(tuple,list)) else x for x in mods}
    elif  not isinstance(models,(tuple,list)):
        mods=[models]

    def _pool_results_to_dict(modResults):
        """
        This function is just used in the parallel processing in order to "share" a data structure between threads.
        :param modResults: The object returned by each call to _fit_data during multiprocessing.
        :type modResults: tuple (Name of current model (including version if exists),name of the sncosmo model,version
                    of sncosmo model,list of tuples containing index and dcurve.fit[modname] for each dcurve in curves)
        :return:None, but updates args
        """
        if modResults:
            modName, source, version, results = modResults
            for i, tempFit in results:
                args['curves'][i].fits.modelFits[modName] = tempFit
                #args['curves'][i].fits.modelFits[modName].model._source = sncosmo.get_source(source, version=version)
                #todo: right now there is not astropy table in the model curve object because of pickling issues, decide if this is a bad idea

    for curve in args['curves']:
        curve.fits.modelFits=dict((mod,newDict()) if not isinstance(mod,(list,tuple)) else (mod[0]+'_'+mod[1],newDict()) for mod in mods)
    #set up parallel processing
    p = Pool(processes=multiprocessing.cpu_count())
    #run each model in parallel and keep track using _pool_results_to_dict
    fits=[]
    for x in p.imap_unordered(_fit_data_wrap,[(x,args) for x in mods]):
        _pool_results_to_dict(x)
        #fits.append(x)
    p.close()
    print('done')
    #print(args['curves'][1].fit['salt2'].res.errors)
    #fig=plt.figure()
    #plt.plot(args['curves'][0].fits.modelFits['salt2']['sdssi'].time,args['curves'][0].fits.modelFits['salt2']['sdssi'].fluxes)
    #plt.show()
    sncosmo.plot_lc(args['curves'][0].table,model=args['curves'][0].fits.modelFits['salt2'].model,errors=args['curves'][0].fits.modelFits['salt2'].res.errors)
    plt.show()

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
    print(mod)
    smod = sncosmo.Model(source=source)
    params=args['params'] if args['params'] else [x for x in smod.param_names]


    for i,dcurve in enumerate(args['curves']):
        if not np.any([smod.bandoverlap(band) for band in dcurve.bands]):
            continue
        dcurve.fits.method=args['method']
        dcurve.fits.bounds=args['bounds'] if args['bounds'] else {}
        dcurve.fits.ignore=args['ignore'] if args['ignore'] else []
        dcurve.fits.constants = args['constants'] if args['constants'] else {x: y for x, y in zip(dcurve.meta.keys(),dcurve.meta.values()) if x != 'info'}
        dcurve.fits.spline = args['spline']
        dcurve.fits.poly = args['poly']
        dcurve.fits.micro = args['micro']
        no_bound = {x for x in params if x in _needs_bounds and x not in dcurve.fits.bounds.keys() and x not in dcurve.fits.constants.keys()}
        if no_bound:
            params=list(set(params)-no_bound)
        params= [x for x in params if x not in dcurve.fits.ignore and x not in dcurve.fits.constants.keys()]
        dcurve.fits.params = params
        if dcurve.fits.constants:
            constants = {x: dcurve.fits.constants[x] for x in dcurve.fits.constants.keys() if x in smod.param_names}
            smod.set(**constants)

        dcurve.fits.modelFits[modName].res, dcurve.fits.modelFits[modName].model = args['sn_func'][args['method']](dcurve.table, smod, dcurve.fits.params, dcurve.fits.bounds,verbose=False, **args['props'])


        '''
        try:
            dcurve=_snmodel_to_flux(dcurve,modName)
        except RuntimeError:
            print('woops')
            continue
        '''



        dcurve = _snmodel_to_flux(dcurve, modName)
        #dcurve.fits.modelFits[modName].model._source = None
    return((modName,mod,version,[(i,dcurve.fits.modelFits[modName]) for i,dcurve in enumerate(args['curves'])]))

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
