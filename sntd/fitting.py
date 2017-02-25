import inspect
import sncosmo
import sntd
import sys

import numpy as np

from .io import _get_default_prop_name

__all__=['fit_data']

_needs_bounds={'z'}


class fitDict(dict):
    #todo document this class
    def __init__(self,curves,bands=None,method='minuit', models=None, params=None, bounds=None, ignore=None, constants=None,
             spline=False, poly=False, micro=True):
        """
        Constructor for curveDict class. Inherits from the dictionary class, and is the main object of organization used by SNTD.
        :param curves: Object containing list of lightcurves to fit
        :type curves: list of ~io.curveDict
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
        :param poly: If you'd like to fit a polynomial to the data instead of an sncosmo template,optional
        :type poly: Boolean
        :param micro: If you'd like to include microlensing, optional
        :type micro: Boolean
        """
        super(fitDict, self).__init__() #init for the super class
        self.curves=curves
        self.bands = bands if bands else set(np.ndarray.flatten(np.asarray([x for x in [y.bands for y in self.curves]])))
        for band in self.bands:
            self[band]=[]
            for dcurve in self.curves:
                self[band].append(fit())

        self.meta = {k: v if k != 'info' else '' for d in [x.meta for x in self.curves] for k, v in d.items()}
        for d in [x.meta for x in self.curves]:
            self.meta['info']+=d['info']

        """@type: dict
            @ivar: The metadata for the fitDict object, intialized with the curveDict metadata. It's
            populated by added _metachar__ characters into the header of your data file.
        """

        self.method=method
        self.models=[models] if models and not isinstance(models,(tuple,list)) else models
        self.params=[] if not params else params
        self.bounds={} if not bounds else bounds
        self.ignore=[] if not ignore else ignore
        self.constants=constants if constants else {x:y for x,y in zip(self.meta.keys(),self.meta.values()) if x!='info'}
        self.spline=spline
        self.poly=poly
        self.micro=micro
        self.telescopename=curves.telescopename
        """
        @type: str
        @ivar: Name of the telescope that the data were gathered from
        """
        self.object=curves.object
        """
        @type: str
        @ivar: Object of interest
        """

    #these three functions allow you to access the curveDict via "dot" notation
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

    def __str__(self):
        """
        A replacement for the print method of the class, so that when you run print(curveDict()), this is how it shows
        up.
        """
        print('Telescope: %s' % self.telescopename)
        print('Object: %s' % self.object)
        print('Number of fitted bands: %d' % len(self.bands))
        print('Method: %s'%self.method)
        if self.spline:
            print('Type: Spline')
        elif self.poly:
            print('Type: Polynomial')
        else:
            if self.models:
                print('Models: {0}'.format(', '.join(self.models)))
            else:
                print('Running all modelss.')
        print('')
        print('Metadata:')
        print('\n'.join('{}:{}'.format(*t) for t in zip([x for x in self.meta.keys() if x not in self.constants.keys()],
                                                        [y for y in self.meta.values() if
                                                         y not in self.constants.values()])))
        print('')
        if self.params:
            print('Models Parameters: {0}'.format(', '.join(self.params)))
        else:
            print('Varying all non-constant default models parameters.')
        print('')
        print('Constants:')
        if self.constants:
            print('\n'.join('{}={}'.format(*t) for t in zip(self.constants.keys(), self.constants.values())))
        else:
            print('None')
        if self.ignore:
            if isinstance(self.ignore,list):
                print('Ignoring: {0}'.format(', '.join(self.ignore)))
            else:
                print('Ignoring: %s'%self.ignore)
        print('')
        #todo print all of the fit information once that's done
        '''
        for band in self.bands:
            print('------------------')
            print('Band: %s' % band)
            print('Date Range: %.5f->%.5f' % (
                min(self.table[self.table[_get_default_prop_name('band')] == band][_get_default_prop_name('time')]),
                max(self.table[self.table[_get_default_prop_name('band')] == band][_get_default_prop_name('time')])))
            print('Number of points: %d' % len(self[band].table))
        '''
        return '------------------'


class fit(dict):
    '''
    This is going to be a fits object that will contain a dictionary
    of modelss with keys as models names. Each value will then be a fitted models
    class from sncosmo.

    If microlensing is added, I guess run through pycs first to get initial time
    delays and microlensing effects, write out the result, and run through sncosmo
    to get best fit models?
    '''
    def __init__(self,dcurve,bands=None,method='minuit', models=None, params=None, bounds={}, ignore=None, constants={},
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

        for band in self.bands:
            self[band]=sntd.curve()
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
    args=locals()
    curves=[curves] if not isinstance(curves,(tuple,list)) else curves
    args['bands'] = {x for x in bands} if bands else set(np.hstack(np.array([x for x in [list(y.bands) for y in curves]])))
    if not args['bands']:
        raise(RuntimeError,"You don't have any bands to analyze!")
    #fits=fitDict(**{x:args[x] for x in args.keys() if x != 'kwargs'})
    mods = sncosmo.models._SOURCES._loaders.keys() if not models else models
    sn_func = {'minuit': sncosmo.fit_lc, 'mcmc': sncosmo.mcmc_lc, 'nest': sncosmo.nest_lc}
    props = {x: kwargs[x] for x in kwargs.keys() if
             x in [y for y in inspect.getargspec(sn_func[method])[0]] and x != 'verbose'}
    for dcurve in curves:  # here dcurve=curveDict
        dcurve.fit = fit(dcurve, **args)
    printed=False
    for model in mods:
        if isinstance(model, tuple):
            model = model[0]
        source = sncosmo.get_source(model)
        mod = sncosmo.Model(source=source)
        """
        tmin=mod.mintime()
        tmax=mod.maxtime()
        tgrid = np.linspace(tmin, tmax, int(tmax - tmin) + 1)
        mflux = mod.bandflux('F160W', tgrid)

        fig = plt.figure()
        ax = plt.gca()
        ax.plot(tgrid, mflux)
        plt.show()
        """
        if not params:
            if len(mods) == 1:
                print(
                "Did not supply params, using default model params. (These are: '{1}' for model {0!r})".format(model,
                                                                                                               "', '".join(
                                                                                                                   mod.param_names)))
            params = {x for x in mod.param_names}
        else:
            params=set(params)
        bands=set()
        for dcurve in curves:
            no_bound = {x for x in params if
                        x in _needs_bounds and x not in dcurve.fit.bounds.keys() and x not in dcurve.fit.constants.keys()}
            if no_bound:
                if not printed:
                    print("Ignoring following parameter(s), didn't have required bounds: {0}".format(', '.join(no_bound)))
                    print('')
                params=list(params-no_bound)
                printed=True
            params= [x for x in params if x not in dcurve.fit.ignore and x not in dcurve.fit.constants.keys()]
            dcurve.fit.params = params
            if dcurve.fit.constants:
                constants = {x: dcurve.fit.constants[x] for x in dcurve.fit.constants.keys() if x in mod.param_names}
                mod.set(**constants)
            dcurve.fit.res, dcurve.fit.model = sn_func[method](dcurve.table, mod, dcurve.fit.params, dcurve.fit.bounds,verbose=False, **props)
            dcurve=_snmodel_to_flux(dcurve)
            bands.update(dcurve.fit.bands)
        """
        This is the point where all light curves have their own fit object for each band containing the flux/mag data for the model template,
        ready to be sent to pycs if that's what we want to do. bands contains all of the bands for all of the curves in a set.
        """
        #todo figure out if we need to do microlensing before time delay measurements, then continue the modeling process
        """
        fig = plt.figure()
        ax = plt.gca()
        try:
            ax.plot(dcurve.fit['F160W'].curve.time-dcurve.fit.model.parameters[1], dcurve.fit['F160W'].curve.fluxes)
            plt.show()
            sncosmo.plot_lc(dcurve.table,dcurve.fit.model)
            plt.show()
        except:
            ax.plot(dcurve.fit['sdssr'].curve.time-dcurve.fit.model.parameters[1], dcurve.fit['sdssr'].curve.fluxes)
            plt.show()
            sncosmo.plot_lc(dcurve.table,dcurve.fit.model)
            plt.show()
        """
        sys.exit()


def _snmodel_to_flux(dcurve):
    for band in dcurve.fit.bands & dcurve.bands:
        dcurve.fit[band]=sntd.curve(band=band,zp=dcurve[band].zp,zpsys=dcurve[band].zpsys)
        tmin = []
        tmax = []
        tmin.append(np.min(dcurve[band].table[_get_default_prop_name('time')]) - 10)
        tmax.append(np.max(dcurve[band].table[_get_default_prop_name('time')]) + 10)
        tmin.append(dcurve.fit.model.mintime())
        tmax.append(dcurve.fit.model.maxtime())
        tmin = min(tmin)
        tmax = max(tmax)
        tgrid = np.linspace(tmin, tmax, int(tmax - tmin) + 1)
        mflux = dcurve.fit.model.bandflux(band, tgrid, zp=dcurve[band].zp, zpsys=dcurve[band].zpsys)
        mmag = -2.5 * np.log10(mflux) + dcurve[band].zp
        dcurve.fit[band].mags=mmag
        dcurve.fit[band].fluxes = mflux
        dcurve.fit[band].time = tgrid
        dcurve.fit[band].magerrs = mmag*.1
        dcurve.fit[band].fluxerrs = mflux*.1
    return dcurve




    """
    list(filter(lambda x: mod.bandoverlap(x),fits.bands)):
    Okay I think it makes sense to do it this way. Go through each curveDict object and run sncosmo on all bands
    at the same time. Each time, keep track of the model result for each band in a fits object. Then send the fits
    object to the pycs fitting function. It will take all of the model results for each band one at a time and send them
    to pycs, which will try to fit them to the data instead of a spline. This means I'll need a fits object with a
    hierarchy looking something like:
    fitsDict-->band keys, each band has a fit object--> each fit object contains model results for that band?

    sn_func={'minuit':sncosmo.fit_lc,'mcmc':sncosmo.mcmc_lc,'nest':sncosmo.nest_lc}
    props = {x: kwargs[x] for x in kwargs.keys() if x in [y for y in inspect.getargspec(sn_func[method])[0]] and x !='verbose'}

    models= sncosmo.models._SOURCES._loaders.keys()
    for dcurve in curves:
        printed=False
        for model in models:
            if isinstance(model,tuple):
                model=model[0]
            source = sncosmo.get_source(model)
            mod = sncosmo.Model(source=source)

            if not params:
                if len(models)==1:
                    print("Did not supply params, using default model params. (These are: '{1}' for model {0!r})".format(model,
                                                                                                                     "', '".join(
                                                                                                                         mod.param_names)))
                params = [x for x in mod.param_names if x != 'z']
                if 'z' in mod.param_names:
                    if args.get('bounds') and 'z' in args.get('bounds').keys():
                        params.append('z')
                    elif 'z' not in constants and 'z' not in ignore and not printed:
                        print("Ignoring 'z' parameter because no bounds were given (required)")
                        printed=True

            fits.params = [x for x in params if x not in ignore and x not in constants.keys()]
            if fits.constants:
                fits.constants={x:fits.constants[x] for x in constants.keys() if x in mod.param_names}
                mod.set(**fits.constants)
            fits[band].res, fits[band].model = sn_func[method](fits[band].dcurve.table, mod, fits.params,fits.bounds,verbose=False, **props)
            #sncosmo.plot_lc(fits[band].dcurve.table,model=fits[band].model)


            _py_fit(fits[band],band,micro,**kwargs)
            '''
            fig = plt.figure()
            ax = plt.gca()
            ax.plot(tgrid-fits[band].model.parameters[1], mflux)
            plt.show()
            '''

def _py_fit(fit,band,micro,**kwargs):
    tmin = []
    tmax = []
    tmin.append(np.min(fit.dcurve.table[_get_default_prop_name('time')]) - 10)
    tmax.append(np.max(fit.dcurve.table[_get_default_prop_name('time')]) + 10)
    tmin.append(fit.model.mintime())
    tmax.append(fit.model.maxtime())
    tmin = min(tmin)
    tmax = max(tmax)
    tgrid = np.linspace(tmin, tmax, int(tmax - tmin) + 1)
    mflux = fit.model.bandflux(band, tgrid, zp=fit[band].dcurve.zp, zpsys=fit[band].dcurve.zpsys)
    mmag = -2.5 * np.log10(mflux) + fit[band].dcurve.zp

"""