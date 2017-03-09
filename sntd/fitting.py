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
from contextlib import contextmanager
from copy import deepcopy
from .io import _get_default_prop_name

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


class fit(dict):
    '''
    This is going to be a fits object that will contain a dictionary
    of models with keys as models names. Each value will then be a fitted models
    class from sncosmo.

    If microlensing is added, I guess run through pycs first to get initial time
    delays and microlensing effects, write out the result, and run through sncosmo
    to get best fit models?
    '''
    def __init__(self,dcurve,bands=None,method='minuit', models=None, params=None, bounds=None, ignore=None, constants=None,
             spline=False, poly=False, micro='spline', **kwargs):
        super(fit, self).__init__() #init for the super class
        self.bands={band for band in bands} if bands else bands
        self.method=method
        self.params=params if params else []
        self.bounds=bounds if bounds else {}
        self.ignore=ignore if ignore else []
        self.constants = constants if constants else {x: y for x, y in zip(dcurve.meta.keys(), dcurve.meta.values()) if
                                                      x != 'info'}
        if not self.constants:
            self.constants={}
        self.spline=spline
        self.poly=poly
        self.micro=micro
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
    args['bands'] = {x for x in bands} if bands else set(
        np.hstack(np.array([x for x in [list(y.bands) for y in args['curves']]])))
    if not args['bands']:
        raise (RuntimeError, "You don't have any bands to analyze!")
    mods = sncosmo.models._SOURCES._loaders.keys() if not models else models
    args['sn_func'] = {'minuit': sncosmo.fit_lc, 'mcmc': sncosmo.mcmc_lc, 'nest': sncosmo.nest_lc}
    args['props'] = {x: kwargs[x] for x in kwargs.keys() if
             x in [y for y in inspect.getargspec(args['sn_func'][method])[0]] and x != 'verbose'}

    def _pool_results_to_dict(modResults):
        modName,source,modResults=modResults
        for i,tempFit in modResults:
            if not hasattr(args['curves'][i],'fit'):
                args['curves'][i].fit=fit(args['curves'][i],**args)
            args['curves'][i].fit[modName] = tempFit
            args['curves'][i].fit[modName].model._source=sncosmo.get_source(source)

    mods={x[0] for x in mods}
    p=Pool(processes=multiprocessing.cpu_count())
    for x in p.imap(_fit_data,[(x,args) for x in mods]):
        _pool_results_to_dict(x)


#todo decide about multiple versions of model
def _fit_data(args):
    warnings.simplefilter("ignore")
    mod=args[0]
    args=args[1]
    if isinstance(mod, tuple):
        version = mod[1]
        mod = mod[0]
    else:
        version = None
    modName=mod+'_'+version if version and version != '1.0' else deepcopy(mod)
    source=sncosmo.get_source(mod)
    smod = sncosmo.Model(source=source)
    params=args['params'] if args['params'] else [x for x in smod.param_names]
    bands=set()
    for dcurve in args.get('curves'):
        dcurve.fit = fit(dcurve, **args)
        dcurve.fit[modName]=newDict()
        no_bound = {x for x in params if
                    x in _needs_bounds and x not in dcurve.fit.bounds.keys() and x not in dcurve.fit.constants.keys()}
        if no_bound:
            params=list(params-no_bound)
        params= [x for x in params if x not in dcurve.fit.ignore and x not in dcurve.fit.constants.keys()]
        dcurve.fit.params = params
        if dcurve.fit.constants:
            constants = {x: dcurve.fit.constants[x] for x in dcurve.fit.constants.keys() if x in smod.param_names}
            smod.set(**constants)
        dcurve.fit[modName].res, dcurve.fit[modName].model = args.get('sn_func')[args.get('method')](dcurve.table, smod, dcurve.fit.params, dcurve.fit.bounds,verbose=False, **args.get('props'))
        try:
            dcurve=_snmodel_to_flux(dcurve,modName)
        except:
            continue
        bands.update(dcurve.fit.bands)
        dcurve.fit[modName].model._source = None
    return((modName,mod,[(i,dcurve.fit[modName]) for i,dcurve in enumerate(args['curves'])]))

#todo figure out how to deal with negative flux from model flux to mag
def _snmodel_to_flux(dcurve,modName):
    warnings.simplefilter("ignore")
    for band in dcurve.fit.bands & dcurve.bands:
        if not dcurve.fit[modName].model.bandoverlap(band):
            continue
        dcurve.fit[modName][band]=sntd.curve(band=band,zp=dcurve[band].zp,zpsys=dcurve[band].zpsys)
        tmin = []
        tmax = []
        tmin.append(np.min(dcurve[band].table[_get_default_prop_name('time')]) - 10)
        tmax.append(np.max(dcurve[band].table[_get_default_prop_name('time')]) + 10)
        tmin.append(dcurve.fit[modName].model.mintime())
        tmax.append(dcurve.fit[modName].model.maxtime())
        tmin = min(tmin)
        tmax = max(tmax)
        tgrid = np.linspace(tmin, tmax, int(tmax - tmin) + 1)
        mflux = dcurve.fit[modName].model.bandflux(band, tgrid, zp=dcurve[band].zp, zpsys=dcurve[band].zpsys)
        mmag = -2.5 * np.log10(mflux) + dcurve[band].zp
        dcurve.fit[modName][band].mags=mmag
        dcurve.fit[modName][band].fluxes = mflux
        dcurve.fit[modName][band].time = tgrid
        dcurve.fit[modName][band].magerrs = mmag*.1
        dcurve.fit[modName][band].fluxerrs = mflux*.1
    return dcurve
