import inspect,sncosmo,os,sys,warnings,pyParz,math,multiprocessing
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy,copy
from scipy import stats
from astropy.table import Table
import nestle
from collections import OrderedDict
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
import scipy
from sncosmo import nest_lc
from scipy.stats import rv_continuous

from .util import *
from .curve_io import _sntd_deepcopy
from .models import *
from .ml import *
#from .sncosmo_fitting import nest_lc

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
			 method='parallel',t0_guess=None,refModel=None,effect_names=[],effect_frames=[],fitting_method='nest',
			 dust=None,flip=False,guess_amplitude=True,seriesError=None,showPlots=False,microlensing=None,fitOrder=None,
			 fit_prior=None,
			 kernel='RBF',refImage='image_1',nMicroSamples=100,color_curve=None,verbose=True,**kwargs):

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
		Needs to be 'parallel', 'series', or 'color'
	t0_guess: :class:`dict`
		Dictionary with image names (i.e. 'image_1','image_2') as keys and a guess for time of peak as values
	refModel: :class:`scipy.interpolate.interp1d` or :class:`sncosmo.Model`
		If doing a series or color fit, a reference model that will be used to fit all the data at once is required.
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
	seriesError: :class:`scipy.interpolate.interp1d`
		If doing a series or color fit, this is the uncertainty on the reference model.
	showPlots: Boolean
		If true, :func:`sncosmo.plot_lc` function is called during the fitting
	microlensing: str
		If None microlensing is ignored, otherwise should be str (e.g. achromatic, chromatic)
	kernel: str
		The kernel to use for microlensing GPR
	seriesGrids: :class:`dict`
		A dictionary with 'td' or 'mu' as keys and tuples with additive bounds as values
	refImage: str
		The name of the image you want to be the reference image (i.e. image_1,image_2, etc.)
	nMicroSamples: int
		The number of pulls from the GPR posterior you want to use for microlensing uncertainty estimation
	color_curve: :class:`astropy.Table`
		A color curve to define the relationship between bands for parameterized light curve model.
	verbose: Boolean
		Turns on/off the verbosity flag
	Returns
	-------
	fitted_curveDict: :class:`~sntd.curve_io.curveDict`
		The same curveDict that was passed to fit_data, but with new fits and time delay measurements included
	Examples
	--------
	>>> fitCurves=sntd.fit_data(myMISN2,snType='Ia', models='salt2-extended',bands=['F110W','F125W'],
		params=['x0','x1','t0','c'],constants={'z':1.33},bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},
		method='parallel',microlensing=None)

	"""

	#get together user arguments
	args = locals()
	for k in kwargs.keys():
		args[k]=kwargs[k]

	if isinstance(curves,(list,tuple)):
		args['curves']=[]
		for i in range(len(curves)):
			temp=_sntd_deepcopy(curves[i])
			temp.nsn=i+1
			args['curves'].append(temp)
		args['parlist']=True
	else:
		args['curves']=_sntd_deepcopy(curves)
		args['parlist']=False

	if method !='color':
		args['bands'] = [bands] if bands is not None and not isinstance(bands,(tuple,list,np.ndarray)) else bands
		#sets the bands to user's if defined (set, so that they're unique), otherwise to all the bands that exist in curves

		args['bands'] = list(set(bands)) if bands is not None else None

		args['bands'] = list(curves.bands) if not isinstance(curves,(list,tuple,np.ndarray)) else list(curves[0].bands)
	elif len(args['bands'])!=2:
		print('Must provide exactly 2 bands for color curve fitting.')
		sys.exit(1)

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
	args['mods']=mods
	#sncosmo fitting function for model determination
	args['sn_func'] = {'minuit': sncosmo.fit_lc, 'mcmc': sncosmo.mcmc_lc, 'nest': nest_lc}

	#get any properties set in kwargs that exist for the defined fitting function
	args['props'] = {x: kwargs[x] for x in kwargs.keys() if
					 x in [y for y in inspect.signature(args['sn_func'][fitting_method]).parameters] and x != 'verbose'}

	if method not in ['parallel','series','color']:
		raise RuntimeError('Parameter "method" must be "parallel","series", or "color".')
	if microlensing is not None and method !='parallel':
		print('Microlensing uncertainty only set up for parallel right now, switching to parallel method...')
		method='parallel'
	
	if method=='parallel':
		if args['parlist']:
			par_arg_vals=[]
			for i in range(len(args['curves'])):
				temp_args={}
				for par_key in ['snType','bounds','constants','t0_guess','refModel','color_curve','seriesGrids']:
					if isinstance(args[par_key],(list,tuple,np.ndarray)):
						temp_args[par_key]=args[par_key][i] 
				for par_key in ['bands','models','ignore','params']:
					if isinstance(args[par_key],(list,tuple,np.ndarray)) and np.any([isinstance(x,(list,tuple,np.ndarray)) for x in args[par_key]]):
						temp_args[par_key]=args[par_key][i] 
				par_arg_vals.append([args['curves'][i],temp_args])
			curves=pyParz.foreach(par_arg_vals,_fitparallel,[args])

		else:
			curves=_fitparallel(args)
	elif method=='series':
		if args['parlist']:

			par_arg_vals=[]
			for i in range(len(args['curves'])):
				temp_args={}
				for par_key in ['snType','bounds','constants','t0_guess','refModel','color_curve','seriesGrids']:
					if isinstance(args[par_key],(list,tuple,np.ndarray)):
						temp_args[par_key]=args[par_key][i] 
				for par_key in ['bands','models','ignore','params']:
					if isinstance(args[par_key],(list,tuple,np.ndarray)) and np.any([isinstance(x,(list,tuple,np.ndarray)) for x in args[par_key]]):
						temp_args[par_key]=args[par_key][i] 
				par_arg_vals.append([args['curves'][i],temp_args])
			curves=pyParz.foreach(par_arg_vals,_fitseries,[args])
		else:
			curves=_fitseries(args)

	else:
		if args['parlist']:
			par_arg_vals=[]
			for i in range(len(args['curves'])):
				temp_args={}
				for par_key in ['snType','bounds','constants','t0_guess','refModel','color_curve','seriesGrids']:
					if isinstance(args[par_key],(list,tuple,np.ndarray)):
						temp_args[par_key]=args[par_key][i] 
				for par_key in ['bands','models','ignore','params']:
					if isinstance(args[par_key],(list,tuple,np.ndarray)) and np.any([isinstance(x,(list,tuple,np.ndarray)) for x in args[par_key]]):
						temp_args[par_key]=args[par_key][i] 
				par_arg_vals.append([args['curves'][i],temp_args])
			curves=pyParz.foreach(par_arg_vals,_fitColor,[args])
		else:
			curves=_fitColor(args)

	return curves

def _fitColor(all_args):
	if isinstance(all_args,(list,tuple,np.ndarray)):
		curves,args=all_args
		if isinstance(args,list):
			args=args[0]
		if isinstance(curves,list):
			curves,single_par_vars=curves
			for key in single_par_vars:
				args[key]=single_par_vars[key]

		args['curves']=curves
		if args['verbose']:
			print('Fitting MISN number %i...'%curves.nsn)
	else:
		args=all_args
	args['bands']=list(args['bands'])
	if len(args['bands'])!=2:
		raise RuntimeError("If you want to analyze color curves, you need two bands!")

	imnums=[x[-1] for x in args['curves'].images.keys()]
	nimage=len(imnums)
	snParams=['t_%s'%i for i in imnums]
	all_vparam_names=np.append(args['params'],
							   snParams).flatten()

	ims=list(args['curves'].images.keys())


	for param in all_vparam_names:
		if param not in args['bounds'].keys():
			args['bounds'][param]=np.array(args['bounds']['td'])#+time_delays[im]

	finallogz=-np.inf
	if args['dust'] is not None:
		if isinstance(args['dust'],str):
			dust_dict={'CCM89Dust':sncosmo.CCM89Dust,'OD94Dust':sncosmo.OD94Dust,'F99Dust':sncosmo.F99Dust}
			dust=dust_dict[args['dust']]()
		else:
			dust=args['dust']
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


	for mod in np.array(args['models']).flatten():



		source=sncosmo.get_source(mod)
		tempMod = sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
		tempMod.set(**args['constants'])
		if not args['curves'].color.table:
			args['curves'].color_table(args['bands'][0],args['bands'][1],referenceImage=args['refImage'],static=False,model=tempMod)

			guess_amplitude=False
		else:
			guess_amplitude=True
		tempMod.set(t0=args['curves'].color.meta['reft0'])

		args['bounds'][tempMod.param_names[2]]=(.1*tempMod.parameters[2],10*tempMod.parameters[2])

		res,td_res,td_err,model=nest_color_lc(args['curves'].color.table,tempMod,nimage,color=args['bands'],
											  				bounds=args['bounds'],
															 vparam_names=all_vparam_names,ref=args['refImage'],
															 refModel=args['refModel'],
															 guess_amplitude_bound=guess_amplitude,
															 minsnr=args.get('minsnr',5.),priors=args.get('priors',None),ppfs=args.get('ppfs',None),
															 method=args.get('nest_method','single'),maxcall=args.get('maxcall',None),
															 modelcov=args.get('modelcov',None),rstate=args.get('rstate',None),
															 maxiter=args.get('maxiter',None),npoints=args.get('npoints',100))

		if finallogz<res.logz:
			finalres,finaltd_res,finaltd_err,finalmodel=res,td_res,td_err,model
			time_delays=args['curves'].color.meta['td']


	args['curves'].color.time_delays=dict([])
	args['curves'].color.magnifications=dict([])
	args['curves'].color.magnification_errors=dict([])
	args['curves'].color.time_delay_errors=dict([])

	#sys.exit()
	for k in args['curves'].images.keys():
		if k==args['refImage']:
			args['curves'].color.time_delays[k]=0
			args['curves'].color.time_delay_errors[k]=0
		else:
			args['curves'].color.time_delays[k]=time_delays[k]+(finaltd_res['t_'+k[-1]]-finaltd_res['t_'+args['refImage'][-1]])
			args['curves'].color.time_delay_errors[k]=np.sqrt(finaltd_err['t_'+k[-1]]**2+finaltd_err['t_'+args['refImage'][-1]]**2)


	args['curves'].color_table(args['bands'][0],args['bands'][1],time_delays=args['curves'].color.time_delays)


	args['curves'].color.fits=newDict()
	args['curves'].color.fits['model']=finalmodel
	args['curves'].color.fits['res']=finalres

	return args['curves']

def nest_color_lc(data,model,nimage,color, vparam_names,bounds,guess_amplitude_bound=False,ref='image_1',
				   minsnr=5.,refModel=False,band=None, priors=None, ppfs=None, npoints=100, method='single',
				   maxiter=None, maxcall=None, modelcov=False, rstate=None,
				   verbose=False, warn=True,**kwargs):

	# experimental parameters
	tied = kwargs.get("tied", None)

	vparam_names=list(vparam_names)
	if ppfs is None:
		ppfs = {}
	if tied is None:
		tied = {}

	if guess_amplitude_bound:
		guess_t0,guess_amp=sncosmo.fitting.guess_t0_and_amplitude(sncosmo.photdata.photometric_data(data[data['image']==ref]),
																  model,minsnr)

		#print(guess_t0,guess_amp)
		model.set(t0=guess_t0)
	#bounds['x0']=(.1*guess_amp,10*guess_amp)

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


	model_param_names=[x for x in vparam_names[:len(vparam_names)-nimage]]
	model_idx = np.array([vparam_names.index(name) for name in model_param_names])
	td_params=[x for x in vparam_names[len(vparam_names)-nimage:] if x[0]=='t']
	td_idx=np.array([vparam_names.index(name) for name in td_params])


	im_indices=[np.where(data['image']==i)[0] for i in np.unique(data['image'])]



	def chisq_likelihood(parameters):
		model.set(**{model_param_names[k]:parameters[model_idx[k]] for k in range(len(model_idx))})
		all_data=deepcopy(data)

		for i in range(len(im_indices)):
			all_data['time'][im_indices[i]]-=parameters[td_idx[i]]

		model_observations = model.color(color[0],color[1],all_data['zpsys'][0],all_data['time'])


		if modelcov:
			all_cov=None
			cov = np.diag(all_data['%s-%s_err'%(color[0],color[1])]**2)
			for i in range(2):

				_, mcov = model.bandfluxcov(color[i],
											all_data['time'],
										zp=all_data['zp_%s'%color[i]],
										zpsys=all_data['zpsys'])


				if all_cov is None:
					all_cov=copy(mcov)**2
				else:
					all_cov+=copy(mcov)**2


			all_cov=np.sqrt(all_cov)

			cov = cov + all_cov
			invcov = np.linalg.pinv(cov)

			diff = all_data['%s-%s'%(color[0],color[1])]-model_observations
			chisq=np.dot(np.dot(diff, invcov), diff)

		else:
			chisq=np.sum((all_data['%s-%s'%(color[0],color[1])]-model_observations)**2/\
						 (all_data['%s-%s_err'%(color[0],color[1])]**2))
		return chisq

	def loglike(parameters):
		chisq=chisq_likelihood(parameters)
		if not np.isfinite(chisq):
			return -np.inf

		return(-.5*chisq)






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



	params=[weighted_quantile(res.samples[:,i],[.16,.5,.84],res.weights) for i in range(len(vparam_names))]

	model.set(**{model_param_names[k]:params[model_idx[k]][1] for k in range(len(model_idx))})
	model.set(t0=model.get('t0')+params[td_idx[td_params.index('t_'+ref[-1])]][1])
	td_res={}
	td_err={}
	for i in range(len(td_params)):
		td_res[td_params[i]]=params[td_idx[i]][1]
		td_err[td_params[i]]=np.array([params[td_idx[i]][1]-params[td_idx[i]][0],params[td_idx[i]][2]-params[td_idx[i]][1]])

	return res,td_res,td_err,model


def _fitseries(all_args):
	if isinstance(all_args,(list,tuple,np.ndarray)):
		curves,args=all_args
		if isinstance(args,list):
			args=args[0]
		if isinstance(curves,list):
			curves,single_par_vars=curves
			for key in single_par_vars:
				args[key]=single_par_vars[key]

		args['curves']=curves
		if args['verbose']:
			print('Fitting MISN number %i...'%curves.nsn)
	else:
		args=all_args
	args['bands']=list(args['bands'])


	#TODO remove amplitude and t0 parameters if in params

	imnums=[x[-1] for x in args['curves'].images.keys()]
	nimage=len(imnums)
	snParams=[['t_%s'%i,'a_%s'%i] for i in imnums]
	all_vparam_names=np.append(args['params'],
							   snParams).flatten()

	ims=list(args['curves'].images.keys())


	for param in all_vparam_names:
		if param not in args['bounds'].keys():
			if param[:2]=='t_':
				if args['fit_prior'] is not None:
					im=[x for x in ims if x[-1]==param[-1]][0]
					args['bounds'][param]=2*np.array([args['fit_prior'].images[im].param_quantiles['t0'][0]- \
													  args['fit_prior'].images[im].param_quantiles['t0'][1],
													  args['fit_prior'].images[im].param_quantiles['t0'][2]- \
													  args['fit_prior'].images[im].param_quantiles['t0'][1]])
				else:
					args['bounds'][param]=np.array(args['bounds']['td'])#+time_delays[im]
			elif param[:2]=='a_':
				if args['fit_prior'] is not None:

					im=[x for x in ims if x[-1]==param[-1]][0]
					amp=args['fit_prior'].images[im].fits.model.param_names[2]
					args['bounds'][param]=2*np.array([args['fit_prior'].images[im].param_quantiles[amp][0]/\
													  args['fit_prior'].images[im].param_quantiles[amp][1],
													  args['fit_prior'].images[im].param_quantiles[amp][2]/ \
													  args['fit_prior'].images[im].param_quantiles[amp][1]])-1

				else:
					args['bounds'][param]=np.array(args['bounds']['mu'])#*magnifications[im]
		elif args['fit_prior'] is not None:
			args['bounds'][param]=np.array([args['fit_prior'].images[args['refImage']].param_quantiles[param][0]- \
											  args['fit_prior'].images[args['refImage']].param_quantiles[param][1],
											  args['fit_prior'].images[args['refImage']].param_quantiles[param][2]- \
											  args['fit_prior'].images[args['refImage']].param_quantiles[param][1]])+ \
								  args['fit_prior'].images[args['refImage']].param_quantiles[param][1]
	finallogz=-np.inf
	if args['dust'] is not None:
		if isinstance(args['dust'],str):
			dust_dict={'CCM89Dust':sncosmo.CCM89Dust,'OD94Dust':sncosmo.OD94Dust,'F99Dust':sncosmo.F99Dust}
			dust=dust_dict[args['dust']]()
		else:
			dust=args['dust']
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


	for mod in np.array(args['models']).flatten():



		source=sncosmo.get_source(mod)
		tempMod = sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
		tempMod.set(**args['constants'])
		if not args['curves'].series.table:
			if args['fit_prior'] is not None:
				args['curves'].combine_curves(time_delays=args['fit_prior'].time_delays,
					magnifications=args['fit_prior'].magnifications,referenceImage=args['refImage'])
				args['curves'].series.meta['reft0']=args['fit_prior'].images[args['refImage']].fits.model.get('t0')
				args['curves'].series.meta['refamp']=\
					args['fit_prior'].images[args['refImage']].fits.model.get(tempMod.param_names[2])
			else:
				args['curves'].combine_curves(referenceImage=args['refImage'],static=False,model=tempMod)
			guess_amplitude=False
		else:
			guess_amplitude=True
		tempMod.set(t0=args['curves'].series.meta['reft0'])
		tempMod.parameters[2]=args['curves'].series.meta['refamp']

		#args['bounds'][tempMod.param_names[2]]=(.1*tempMod.parameters[2],10*tempMod.parameters[2])
		res,td_res,mu_res,td_err,mu_err,model=nest_series_lc(args['curves'].series.table,tempMod,nimage,bounds=args['bounds'],
									  vparam_names=all_vparam_names,ref=args['refImage'],
									  refModel=args['refModel'],
									  guess_amplitude_bound=guess_amplitude,
									  minsnr=args.get('minsnr',5.),priors=args.get('priors',None),ppfs=args.get('ppfs',None),
									  method=args.get('nest_method','single'),maxcall=args.get('maxcall',None),
									  modelcov=args.get('modelcov',None),rstate=args.get('rstate',None),
									  maxiter=args.get('maxiter',None),npoints=args.get('npoints',100))
		if finallogz<res.logz:
			finalres,finaltd_res,finalmu_res,finaltd_err,finalmu_err,finalmodel=res,td_res,mu_res,td_err,mu_err,model
			time_delays=args['curves'].series.meta['td']
			magnifications=args['curves'].series.meta['mu']


	args['curves'].series.time_delays=dict([])
	args['curves'].series.magnifications=dict([])
	args['curves'].series.magnification_errors=dict([])
	args['curves'].series.time_delay_errors=dict([])

	#sys.exit()
	for k in args['curves'].images.keys():
		if k==args['refImage']:
			args['curves'].series.time_delays[k]=0
			args['curves'].series.magnifications[k]=1
			args['curves'].series.time_delay_errors[k]=0
			args['curves'].series.magnification_errors[k]=0
		else:
			args['curves'].series.time_delays[k]=time_delays[k]+(finaltd_res['t_'+k[-1]]-finaltd_res['t_'+args['refImage'][-1]])
			args['curves'].series.magnifications[k]=magnifications[k]*(finalmu_res['a_'+k[-1]]/finalmu_res['a_'+args['refImage'][-1]])
			args['curves'].series.time_delay_errors[k]=np.sqrt(finaltd_err['t_'+k[-1]]**2+finaltd_err['t_'+args['refImage'][-1]]**2)
			args['curves'].series.magnification_errors[k]=args['curves'].series.magnifications[k]*\
				np.sqrt((finalmu_err['a_'+k[-1]]/finalmu_res['a_'+k[-1]])**2+(finalmu_err['a_'+args['refImage'][-1]]/\
																			  finalmu_res['a_'+args['refImage'][-1]])**2)



	args['curves'].combine_curves(time_delays=args['curves'].series.time_delays,magnifications=args['curves'].series.magnifications,referenceImage=args['refImage'])


	args['curves'].series.fits=newDict()
	args['curves'].series.fits['model']=finalmodel
	args['curves'].series.fits['res']=finalres

	return args['curves']




def nest_series_lc(data,model,nimage,vparam_names,bounds,guess_amplitude_bound=False,ref='image_1',
					 minsnr=5.,refModel=False,band=None, priors=None, ppfs=None, npoints=100, method='single',
					 maxiter=None, maxcall=None, modelcov=False, rstate=None,
					 verbose=False, warn=True,**kwargs):

	# experimental parameters
	tied = kwargs.get("tied", None)

	vparam_names=list(vparam_names)
	if ppfs is None:
		ppfs = {}
	if tied is None:
		tied = {}

	if guess_amplitude_bound:
		guess_t0,guess_amp=sncosmo.fitting.guess_t0_and_amplitude(sncosmo.photdata.photometric_data(data[data['image']==ref]),
															  model,minsnr)

	#print(guess_t0,guess_amp)
		model.set(t0=guess_t0)
		model.parameters[2]=guess_amp
	#bounds['x0']=(.1*guess_amp,10*guess_amp)

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


	model_param_names=[x for x in vparam_names[:len(vparam_names)-nimage*2]]
	model_idx = np.array([vparam_names.index(name) for name in model_param_names])
	td_params=[x for x in vparam_names[len(vparam_names)-nimage*2:] if x[0]=='t']
	td_idx=np.array([vparam_names.index(name) for name in td_params])
	amp_params=[x for x in vparam_names[len(vparam_names)-nimage*2:] if x[0]=='a']
	amp_idx=np.array([vparam_names.index(name) for name in amp_params])

	im_indices=[np.where(data['image']==i)[0] for i in np.unique(data['image'])]

	def chisq_likelihood(parameters):
		model.set(**{model_param_names[k]:parameters[model_idx[k]] for k in range(len(model_idx))})
		all_data=deepcopy(data)

		for i in range(len(im_indices)):
			all_data['time'][im_indices[i]]-=parameters[td_idx[i]]
			all_data['flux'][im_indices[i]]/=parameters[amp_idx[i]]
		model_observations = model.bandflux(all_data['band'],all_data['time'],
											zp=all_data['zp'],zpsys=all_data['zpsys'])

		if modelcov:
			cov = np.diag(all_data['fluxerr']**2)
			_, mcov = model.bandfluxcov(all_data['band'], all_data['time'],
										zp=all_data['zp'], zpsys=all_data['zpsys'])

			cov = cov + mcov
			invcov = np.linalg.pinv(cov)

			diff = all_data['flux']-model_observations
			chisq=np.dot(np.dot(diff, invcov), diff)

		else:
			chisq=np.sum((all_data['flux']-model_observations)**2/(all_data['fluxerr']**2))
		return chisq

	def loglike(parameters):
		chisq=chisq_likelihood(parameters)
		if not np.isfinite(chisq):
			return -np.inf

		return(-.5*chisq)






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



	params=[weighted_quantile(res.samples[:,i],[.16,.5,.84],res.weights) for i in range(len(vparam_names))]
	model.set(**{model_param_names[k]:params[model_idx[k]][1] for k in range(len(model_idx))})
	model.set(t0=model.get('t0')+params[td_idx[td_params.index('t_'+ref[-1])]][1])
	model.parameters[2]=model.parameters[2]*params[amp_idx[amp_params.index('a_'+ref[-1])]][1]
	td_res={}
	mu_res={}
	td_err={}
	mu_err={}
	for i in range(len(td_params)):
		td_res[td_params[i]]=params[td_idx[i]][1]
		td_err[td_params[i]]=np.array([params[td_idx[i]][1]-params[td_idx[i]][0],params[td_idx[i]][2]-params[td_idx[i]][1]])
		mu_res[amp_params[i]]=params[amp_idx[i]][1]
		mu_err[amp_params[i]]=np.array([params[amp_idx[i]][1]-params[amp_idx[i]][0],params[amp_idx[i]][2]-params[amp_idx[i]][1]])


	return res,td_res,mu_res,td_err,mu_err,model

def weighted_quantile(values, quantiles, sample_weight=None,
					  values_sorted=False, old_style=False):
	""" Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
	values = np.array(values)
	quantiles = np.array(quantiles)
	if sample_weight is None:
		sample_weight = np.ones(len(values))
	sample_weight = np.array(sample_weight)
	assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
		'quantiles should be in [0, 1]'

	if not values_sorted:
		sorter = np.argsort(values)
		values = values[sorter]
		sample_weight = sample_weight[sorter]

	weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
	if old_style:
		# To be convenient with numpy.percentile
		weighted_quantiles -= weighted_quantiles[0]
		weighted_quantiles /= weighted_quantiles[-1]
	else:
		weighted_quantiles /= np.sum(sample_weight)
	return np.interp(quantiles, weighted_quantiles, values)

def par_fit_parallel(all_args):
	d,fitDict,args,bestFit,bestRes=all_args
	_,bestFit,bestMod,bounds=fitDict[d]
	tempTable=deepcopy(args['curves'].images[d].table)

	for b in [x for x in np.unique(tempTable['band']) if x not in args['bands']]:
		tempTable=tempTable[tempTable['band']!=b]
	if args['flip']:
		tempTable['flux']=np.flip(tempTable['flux'],axis=0)


	if 'amplitude' not in bounds.keys():
		guess_amp_bounds=True
	else:
		guess_amp_bounds=False



	nest_res,nest_fit=_nested_wrapper(args['curves'],tempTable,bestFit,vparams=bestRes.vparam_names,bounds=bounds,
									  priors=args.get('priors',None), ppfs=args.get('None'), method=args.get('nest_method','single'),
									  maxcall=args.get('maxcall',None), modelcov=args.get('modelcov',False),
									  rstate=args.get('rstate',None),
									  guess_amplitude_bound=guess_amp_bounds,microlensing=args['microlensing'],
									  zpsys=args['curves'].images[d].zpsys,kernel=args['kernel'],
									  maxiter=args.get('maxiter',None),npoints=args.get('npoints',100),nsamples=args['nMicroSamples'])








	return([d,nest_fit,nest_res,len(tempTable)- len(nest_res.vparam_names)])


def _fitparallel(all_args):
	if isinstance(all_args,(list,tuple,np.ndarray)):
		curves,args=all_args
		if isinstance(args,list):
			args=args[0]
		if isinstance(curves,list):
			curves,single_par_vars=curves
			for key in single_par_vars:
				args[key]=single_par_vars[key]

		args['curves']=curves
		if args['verbose']:
			print('Fitting MISN number %i...'%curves.nsn)
	else:
		args=all_args

	fitDict=dict([])
	if 't0' in args['bounds']:
		t0Bounds=copy(args['bounds']['t0'])
	if 'amplitude' in args['bounds']:
		ampBounds=copy(args['bounds']['amplitude'])
	for d in args['curves'].images.keys():
		#print(curves.images[d].simMeta)
		args['curve']=copy(args['curves'].images[d])
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



		args['curves'].images[d].fits=newDict()
		if True or len(args['curve'].table)>63 or len(args['mods'])==1 or args['snType']=='Ia':
			fits=[]
			for mod in args['mods']:
				if mod =='BazinSource' or isinstance(mod,BazinSource):
					fits.append(param_fit(args,mod))


				else:
					if len(args['mods'])==1:
						doFit=False
					else:
						doFit=True
					args['doFit']=doFit
					fits.append(_fit_data_wrap((mod,args)))
		else:
			args['doFit']=True
			fits=pyParz.foreach(args['mods'],_fit_data,args)

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


		fitDict[d]=[fits,bestFit,bestRes,copy(args['bounds'])]

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


	if args['fitOrder'] is None:
		all_SNR=[np.sum(args['curves'].images[d].table['flux']/args['curves'].images[d].table['fluxerr']) \
					for d in np.sort(list(args['curves'].images.keys()))]
		sorted=np.flip(np.argsort(all_SNR))
		args['fitOrder']=np.sort(list(args['curves'].images.keys()))[sorted]


	initial_bounds=deepcopy(args['bounds'])
	first_res=par_fit_parallel([args['fitOrder'][0],fitDict,args,bestFit,bestRes])
	first_params=[weighted_quantile(first_res[2].samples[:,i],[.16,.5,.84],first_res[2].weights)\
				  for i in range(len(first_res[2].vparam_names))]
	args['curves'].images[args['fitOrder'][0]].fits=newDict()
	args['curves'].images[args['fitOrder'][0]].fits['model']=first_res[1]
	args['curves'].images[args['fitOrder'][0]].fits['res']=first_res[2]

	args['curves'].time_delays={args['fitOrder'][0]:0}
	args['curves'].magnifications={args['fitOrder'][0]:1}
	args['curves'].time_delay_errors={args['fitOrder'][0]:0}
	args['curves'].magnification_errors={args['fitOrder'][0]:0}


	t0ind=first_res[2].vparam_names.index('t0')
	ampind=first_res[2].vparam_names.index(first_res[1].param_names[2])

	args['curves'].images[args['fitOrder'][0]].param_quantiles={k:first_params[first_res[2].vparam_names.index(k)] for\
																 k in first_res[2].vparam_names}
	for d in args['fitOrder'][1:]:
		args['curves'].images[d].fits=newDict()
		initial_bounds['t0']=t0Bounds

		params,args['curves'].images[d].fits['model'],args['curves'].images[d].fits['res']\
			=nest_parallel_lc(args['curves'].images[d].table,first_res[1],first_res[2],initial_bounds,
						 	guess_amplitude_bound=True,priors=args.get('priors',None), ppfs=args.get('None'),
						 method=args.get('nest_method','single'),
						 maxcall=args.get('maxcall',None), modelcov=args.get('modelcov',False),
						 rstate=args.get('rstate',None),
						 maxiter=args.get('maxiter',None),npoints=args.get('npoints',1000))

		args['curves'].time_delays[d]=params[t0ind][1]-first_params[t0ind][1]
		args['curves'].magnifications[d]=params[ampind][1]/first_params[ampind][1]
		args['curves'].time_delay_errors[d]=np.sqrt([(params[t0ind][1]-params[t0ind][0])**2+\
													 (first_params[t0ind][1]-first_params[t0ind][0])**2,
										(params[t0ind][2]-params[t0ind][1])**2+ \
										(first_params[t0ind][2]-first_params[t0ind][1])**2])
		args['curves'].magnification_errors[d]=args['curves'].magnifications[d]*\
			np.sqrt([((params[ampind][1]-params[ampind][0])/params[ampind][1])**2+\
					((first_params[ampind][1]-first_params[ampind][0])/first_params[ampind][1])**2,
									 ((params[ampind][2]-params[ampind][1])/params[ampind][1])**2+ \
								   ((first_params[ampind][2]-first_params[ampind][1])/first_params[ampind][1])**2])

		args['curves'].images[d].param_quantiles={k:params[args['curves'].images[d].fits['res'].vparam_names.index(k)]\
												  for k in args['curves'].images[d].fits['res'].vparam_names}
	if args['showPlots']:
		for d in args['curves'].images.keys():
			tempTable=copy(args['curves'].images[d].table)
			for b in [x for x in np.unique(tempTable['band']) if x not in args['bands']]:
				tempTable=tempTable[tempTable['band']!=b]
			tempMod=copy(args['curves'].images[d].fits.model)


			sncosmo.plot_lc(tempTable,model=tempMod)

			#plt.savefig(nest_fit._source.name+'_'+tempTable['band'][0]+'_refs_'+d+'.pdf',format='pdf',overwrik4ite=True)
			#plt.savefig('example_plot_dust_image_'+str(d[-1])+'.png',format='png',overwrite=True)
			plt.show()
			plt.clf()
			plt.close()

	return args['curves']

def nest_parallel_lc(data,model,prev_res,bounds,guess_amplitude_bound=False,
				   minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
				   maxiter=None, maxcall=None, modelcov=False, rstate=None,
				   verbose=False, warn=True,**kwargs):

	# experimental parameters
	tied = kwargs.get("tied", None)

	vparam_names=list(prev_res.vparam_names)
	if ppfs is None:
		ppfs = {}
	if tied is None:
		tied = {}

	model=copy(model)
	if guess_amplitude_bound:
		guess_t0,guess_amp=sncosmo.fitting.guess_t0_and_amplitude(sncosmo.photdata.photometric_data(data),
																  model,minsnr)

		#print(guess_t0,guess_amp)
		model.set(t0=guess_t0)
		model.parameters[2]=guess_amp

		bounds[model.param_names[2]]=(0,10*guess_amp)
		bounds['t0']=np.array(bounds['t0'])+guess_t0

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

	final_priors=[]
	for p in vparam_names:
		if p==model.param_names[2] or p=='t0':
			final_priors.append(lambda x:0)
			continue

		ind=prev_res.vparam_names.index(p)
		temp=posterior('temp',np.min(prev_res.samples[:,ind]),np.max(prev_res.samples[:,ind]))
		samples=np.linspace(np.min(prev_res.samples[:,ind]),
							np.max(prev_res.samples[:,ind]),1000)
		final_priors.append(scipy.interpolate.interp1d(samples,
													   np.log(temp._pdf(samples,prev_res.samples[:,ind],
																		prev_res.weights)),fill_value=-np.inf,
													   bounds_error=False))

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



	def chisq_likelihood(parameters):



		model.set(**{vparam_names[k]:parameters[k] for k in range(len(vparam_names))})
		model_observations = model.bandflux(data['band'],data['time'],
											zp=data['zp'],zpsys=data['zpsys'])

		if modelcov:
			cov = np.diag(data['fluxerr']**2)
			_, mcov = model.bandfluxcov(data['band'], data['time'],
										zp=data['zp'], zpsys=data['zpsys'])

			cov = cov + mcov
			invcov = np.linalg.pinv(cov)

			diff = data['flux']-model_observations
			chisq=np.dot(np.dot(diff, invcov), diff)

		else:
			chisq=np.sum((data['flux']-model_observations)**2/(data['fluxerr']**2))
		return chisq

	def loglike(parameters):
		prior_val=0
		for i in range(len(parameters)):
			temp_prior=final_priors[i](parameters[i])
			if not np.isfinite(temp_prior):
				return -np.inf
			prior_val+=temp_prior

		chisq=chisq_likelihood(parameters)
		if not np.isfinite(chisq):
			return -np.inf

		return(prior_val-.5*chisq)






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



	params=[weighted_quantile(res.samples[:,i],[.16,.5,.84],res.weights) for i in range(len(vparam_names))]

	model.set(**{vparam_names[k]:params[k][1] for k in range(len(vparam_names))})


	return params,model,res

class posterior(rv_continuous):
	"Skewed Normal Distribution"
	def _pdf(self,x,samples,weights):
		pdf,edges=np.histogram(samples,weights=weights,
							   bins=20,density=True)
		func=scipy.interpolate.interp1d([(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)],
										pdf,fill_value=0,bounds_error=False)
		return(func(x))

	def _argcheck(self,*args):
		return True

def _micro_uncertainty(args):
	sample,other=args
	nest_fit,data,colnames,x_pred,vparam_names,bounds,priors=other
	data=Table(data,names=colnames)

	temp_nest_mod=deepcopy(nest_fit)
	tempMicro=AchromaticMicrolensing(x_pred/(1+nest_fit.get('z')),sample,magformat='multiply')
	temp_nest_mod.add_effect(tempMicro,'microlensing','rest')
	tempRes,tempMod=nest_lc(data,temp_nest_mod,vparam_names=vparam_names,bounds=bounds,
							guess_amplitude_bound=True,maxiter=None,npoints=200,priors=priors)

	return float(tempMod.get('t0'))


def _nested_wrapper(curves,data,model,vparams,bounds,priors,guess_amplitude_bound,
					microlensing,zpsys,kernel,maxiter,npoints,nsamples,ppfs, method,
					maxcall, modelcov,
					rstate):

	temp=deepcopy(data)
	vparam_names=deepcopy(vparams)


	if microlensing is not None:
		nest_res,nest_fit=nest_lc(temp,model,vparam_names=vparam_names,bounds=bounds,ppfs=ppfs,
								  guess_amplitude_bound=guess_amplitude_bound,maxiter=maxiter,npoints=npoints,
								  priors=priors,method=method,maxcall=maxcall,modelcov=modelcov,rstate=rstate)




		micro,sigma,x_pred,y_pred,samples=fit_micro(curves,nest_res,nest_fit,temp,zpsys,nsamples,
													micro_type=microlensing,kernel=kernel)


		#temp=deepcopy(data)
		#print(samples.shape)
		#return(nest_res,nest_fit)
		t0s=pyParz.foreach(samples.T,_micro_uncertainty,[nest_fit,np.array(temp),temp.colnames,x_pred,vparam_names,bounds,priors])
		mu,sigma=scipy.stats.norm.fit(t0s)

		nest_res.errors['micro']=np.sqrt(np.abs(nest_fit.get('t0')-mu)**2+(3*sigma)**2)
		bestRes=nest_res
		bestMod=nest_fit


	else:
		bestRes,bestMod=nest_lc(data,model,vparam_names=vparam_names,bounds=bounds,ppfs=ppfs,
								guess_amplitude_bound=guess_amplitude_bound,maxiter=maxiter,npoints=npoints,
								priors=priors,method=method,maxcall=maxcall,modelcov=modelcov,rstate=rstate)

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
		#tempData=tempData[mod>.1]
		residual=tempData['flux']/mod
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


		if False:
			plt.close()
			fig=plt.figure()
			ax=fig.gca()
			for i in range(samples.shape[1]):
				if i==0:
					ax.plot(X, samples[:,i],alpha=.1,label='Posterior Samples',color='b')
				else:
					ax.plot(X, samples[:,i],alpha=.1,color='b')
			ax.errorbar(allTime.ravel(), allResid, allErr, fmt='r.', markersize=10, label=u'Observations')

			ax.plot(X, y_pred - 3 * sigma, '--g')
			ax.plot(X, y_pred + 3 * sigma, '--g',label='$3\sigma$ Bounds')
			ax.plot(X,y_pred,'k-.',label="GPR Prediction")

			ax.set_ylabel('Magnification ($\mu$)')
			ax.set_xlabel('Observer Frame Time (Days)')
			ax.plot(X,curves.images['image_2'].simMeta['microlensing_params'](X/(1+1.33))/np.median(curves.images['image_2'].simMeta['microlensing_params'](X/(1+1.33))),'k',label='True $\mu$-Lensing')
			ax.legend(fontsize=10)
			plt.show()
			#plt.savefig('microlensing_gpr.pdf',format='pdf',overwrite=True)
			plt.close()
			#sys.exit()


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
	dust_dict={'CCM89Dust':sncosmo.CCM89Dust,'OD94Dust':sncosmo.OD94Dust,'F99Dust':sncosmo.F99Dust}
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
	params=args['params'] if args['params'] else [x for x in smod.vparam_names]
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




		res,mod=nest_lc(args['curve'].table,mod,vparam_names=args['params'],bounds=args['bounds'],guess_amplitude_bound=guess_amp_bound,maxiter=1000,npoints=200)
	return({'res':res,'model':mod})

