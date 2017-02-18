import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from .io import _get_default_prop_name,_props
import sys,inspect,sncosmo

__all__=['fit_data']
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
        self.bands = bands if bands else list({x for x in [y.bands for y in curves]})
        for band in self.bands:
            self[band]=[]
            for curve in curves:
                self[band].append(fit(curve[band],meta=curve.meta,))


        self.meta = curves.meta
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


class fit:
    '''
    This is going to be a fits object that will contain a dictionary
    of modelss with keys as models names. Each value will then be a fitted models
    class from sncosmo.

    If microlensing is added, I guess run through pycs first to get initial time
    delays and microlensing effects, write out the result, and run through sncosmo
    to get best fit models?
    '''
    def __init__(self,curve):
        self.curve=curve


def fit_data(curves, bands=None,method='minuit', models=None, params=None, bounds=None, ignore=None, constants=None,
             spline=False, poly=False, micro='spline', **kwargs):
    """
    The main, high-level fitting function.
    :param curves: list of objects containing lightcurves to fit, or a single curve
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

    fits = fitDict(curves, bands, method, models, params, bounds, ignore, constants, spline, poly, micro) if isinstance(
        curves, (tuple, list)) else fitDict([curves], bands, method, models, params, bounds, ignore, constants, spline,
                                            poly, micro)
    if spline:
        if poly:
            raise(RuntimeError,"Can't fit spline and polynomial at the same time!")
        #todo implement a normal pycs spline fit to each band
        return None
    elif poly:
        #todo implement a normal pycs polyfit to each band
        return None
    sn_func={'minuit':sncosmo.fit_lc,'mcmc':sncosmo.mcmc_lc,'nest':sncosmo.nest_lc}
    props = {x: kwargs[x] for x in kwargs.keys() if x in [y for y in inspect.getargspec(sn_func[method])[0]] and x !='verbose'}

    models=fits.models if fits.models else sncosmo.models._SOURCES._loaders.keys()
    for band in fits.bands:
        printed=False
        for model in models:
            if isinstance(model,tuple):
                model=model[0]
            source = sncosmo.get_source(model)
            mod = sncosmo.Model(source=source)

            if not fits.params:
                if len(models)==1:
                    print("Did not supply params, using default model params. (These are: '{1}' for model {0!r})".format(model,
                                                                                                                     "', '".join(
                                                                                                                         mod.param_names)))
                params = [x for x in mod.param_names if x != 'z']
                if 'z' in mod.param_names:
                    if args.get('bounds') and 'z' in args.get('bounds').keys():
                        params.append('z')
                    elif 'z' not in fits.constants and 'z' not in fits.ignore and not printed:
                        print("Ignoring 'z' parameter because no bounds were given (required)")
                        printed=True

            fits.params = [x for x in params if x not in fits.ignore and x not in fits.constants.keys()]
            if fits.constants:
                fits.constants={x:fits.constants[x] for x in fits.constants.keys() if x in mod.param_names}
                mod.set(**fits.constants)
            fits[band].res, fits[band].model = sn_func[method](fits[band].curve.table, mod, fits.params,fits.bounds,verbose=False, **props)
            #sncosmo.plot_lc(fits[band].curve.table,model=fits[band].model)


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
    tmin.append(np.min(fit.curve.table[_get_default_prop_name('time')]) - 10)
    tmax.append(np.max(fit.curve.table[_get_default_prop_name('time')]) + 10)
    tmin.append(fit.model.mintime())
    tmax.append(fit.model.maxtime())
    tmin = min(tmin)
    tmax = max(tmax)
    tgrid = np.linspace(tmin, tmax, int(tmax - tmin) + 1)
    mflux = fit.model.bandflux(band, tgrid, zp=fit[band].curve.zp, zpsys=fit[band].curve.zpsys)
    mmag = -2.5 * np.log10(mflux) + fit[band].curve.zp

