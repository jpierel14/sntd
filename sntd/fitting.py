import inspect,sncosmo,os,sys,warnings,pyParz
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from copy import deepcopy
from scipy.interpolate import splrep,splev
from scipy.optimize import curve_fit
from scipy.stats import chisquare

from .io import _get_default_prop_name
from .simulation import _getAbsFromDist,_getAbsoluteDist
from .util import __dir__
from astropy.table import Table

__thetaSN__=['z','hostebv','screenebv','screenz']
__thetaL__=['t0','amplitude']

__all__=['fit_data_separately','fit_combined_data','colorFit']

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



def fit_combined_data():
    pass

def fit_data_separately(curves, snType='Ia',bands=None,method='minuit', models=None, params=None, bounds=None, ignore=None, constants=None,
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
    args['constants']=constants
    args['effect_names']=kwargs.get('effect_names',[])
    args['effect_frames']=kwargs.get('effect_frames',[])
    args['dust']=kwargs.get('dust',None)
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

    resList=dict([])
    fitDict=dict([])
    for d in curves.images.keys():
        #print(curves.images[d].simMeta)
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
        fitDict[d]=[fits,bestFit,bestRes]
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
        nest_res,nest_fit=sncosmo.nest_lc(curves.images[d].table,bestFit,vparam_names=bestRes.vparam_names,bounds=bounds,guess_amplitude_bound=True,maxiter=1000,npoints=100)
        #sncosmo.plot_lc(data=curves.images[d].table,model=nest_fit,errors=nest_res.errors)
        #plt.show()
        #plt.close()
        resList[d]=nest_res
        curves.images[d].fits=newDict()
        curves.images[d].fits['model']=nest_fit
        curves.images[d].fits['res']=nest_res
        print(curves.images[d].simMeta)


    joint=_joint_likelihood(resList,verbose=True)
    for p in joint.keys():
        for d in curves.images.keys():
            if isinstance(joint[p],dict):
                curves.images[d].fits.model.set(**{p:joint[p][d][0]})
            else:
                curves.images[d].fits.model.set(**{p:joint[p][0]})

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
    modName=mod+'_'+version if version else deepcopy(mod)
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