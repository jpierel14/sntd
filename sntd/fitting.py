import warnings,sncosmo,os,sys,pyParz,pickle,subprocess,glob,math,time,tarfile
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy,copy
from scipy import stats
from astropy.table import Table
import nestle
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
import scipy
import itertools
from sncosmo import nest_lc
from itertools import combinations


from .util import *
from .util import __filedir__,__current_dir__
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



def fit_data(curves=None, snType='Ia',bands=None, models=None, params=None, bounds={}, ignore=None, constants={},ignore_models=[],
			 method='parallel',t0_guess=None,effect_names=[],effect_frames=[],batch_init=None,cut_time=None,force_positive_param=[],
			 dust=None,microlensing=None,fitOrder=None,color_bands=None,color_param_ignore=[],min_points_per_band=3,identify_micro=False,
			 min_n_bands=1,max_n_bands=None,n_cores_per_node=1,npar_cores=4,max_batch_jobs=199,max_cadence=None,fit_colors=None,
			 fit_prior=None,par_or_batch='parallel',batch_partition=None,nbatch_jobs=None,batch_python_path=None,n_per_node=None,fast_model_selection=True,
			 wait_for_batch=False,band_order=None,set_from_simMeta={},guess_amplitude=True,trial_fit=True,clip_data=False,
			 kernel='RBF',refImage='image_1',nMicroSamples=100,color_curve=None,warning_supress=True,
			 verbose=True,**kwargs):

	"""The main high-level fitting function.

	Parameters
	----------
	curves: :class:`~sntd.curve_io.curveDict`
		The curveDict object containing the multiple images to fit.
	snType: str
		The supernova classification
	bands: :class:`~list` of :class:`~sncosmo.Bandpass` or :class:`~str`, or :class:`~sncosmo.Bandpass` or :class:`~str`
		The band(s) to be fit
	models: :class:`~list` of :class:`~sncosmo.Model` or str, or :class:`~sncosmo.Model` or :class:`~str`
		The model(s) to be used for fitting to the data
	params: :class:`~list` of :class:`~str`
		The parameters to be fit for the models inside of the parameter models
	bounds: :class:`dict`
		A dictionary with parameters in params as keys and a tuple of bounds as values
	ignore: :class:`~list` of :class:`~str`
		List of parameters to ignore
	constants: :class:`dict`
		Dictionary with parameters as keys and the constant value you want to set them to as values
	ignore_models: class:`~list`
		List of model names to ignore, usually used if you did not specify the "models" parameter
		and let all models for a given SN type be chosen, but you want to ignore one or more.
	method: :class:`~str` or :class:`~list`
		Needs to be 'parallel', 'series', or 'color', or a list containting one or more of these
	t0_guess: :class:`dict`
		Dictionary with image names (i.e. 'image_1','image_2') as keys and a guess for time of peak as values
	effect_names: :class:`~list` of :class:`~str`
		List of effect names if model contains a :class:`~sncosmo.PropagationEffect`.
	effect_frames: :class:`~list` of :class:`~str`
		List of the frames (e.g. obs or rest) that correspond to the effects in effect_names
	batch_init: :class:`~str`
		A string to be pasted into the batch python file (e.g. extra imports or filters added to sncosmo.)
	cut_time: :class:`~list`
		The start and end (rest frame) phase that you want to fit in, default accept all phases. 
	force_positive_param: :class:`~list`
		Optional list of parameters to always make positive.
	dust: :class:`sncosmo.PropagationEffect`
		An sncosmo dust propagation effect to include in the model
	microlensing: str
		If None microlensing is ignored, otherwise should be str (e.g. achromatic, chromatic)
	fitOrder: :class:`~list`
		The order you want to fit the images if using parallel method (default chooses by npoints/SNR)
	color_bands: :class:`~list`
		If using multiple methods (in batch mode), the subset of bands to use for color fitting.
	color_param_ignore: :class:`~list`
		If using multiple methods, parameters you may want to fit for one method but not 
		for color method (e.g. stretch)
	min_points_per_band: int
		Only accept bands to fit with this number of points fitting other criterion (e.g. minsnr)
	identify_micro: bool
		If True, function is run to attempt to identify bands where microlensing is least problematic.
	min_n_bands: int
		Checks the SN to make sure it has this number of bands (with min_points_per_band in each)
	max_n_bands: int
		The best n bands are chosen from the data. 
	n_cores_per_node: int
		The number of cores to run parallelization on per node
	npar_cores: int
		The number of cores to devote to parallelization
	max_batch_jobs: int 
		The maximum number of jobs allowed by your slurm task manager. 
	max_cadence: int
		To clip each image of a MISN to this cadence
	fit_colors: list
		List of colors to use in color fitting (e.g. ['bessellb-bessellv','bessellb-bessellr'])
	fit_prior: :class:`~sntd.curve_io.curveDict` or bool
		if implementing parallel method alongside others and fit_prior is True, will use output of parallel as prior
		for series/color. If SNTD curveDict object, used as prior for series or color.
	par_or_batch: str
		if providing a list of SNe, par means multiprocessing and batch means sbatch. Must supply other batch
		parameters if batch is chosen, so parallel is default.
	batch_partition: str
		The name of the partition for sbatch command
	nbatch_jobs: int
		number of jobs (10 jobs for 100 light curves is 10 light curves per job)
	batch_python_path: str
		path to python you want to use for batch mode (if different from current)
	n_per_node: int
		Number of SNe to fit per node (in series) in batch mode. If none, just distributes all SNe across the number
		of jobs you have by default. 
	fast_model_selection: bool
		If you are providing a list of models and want the best fit, turning this on will make the fitter choose based
		on a simple minuit fit before moving to the full sntd fitting. If false, each model will be fitted with the full
		sntd fitting and the best will be chosen. 
	wait_for_batch: bool
		if false, submits job in the background. If true, waits for job to finish (shows progress bar) and returns output.
	band_order: :class:`~list`
		If you want colors to be fit in a specific order (e.g. B-V instead of V-B depending on band order)
	set_from_simMeta: :class:`~dict`
		Dictionary where keys are model parameters and values are the corresponding key in the 
		:class:`~sntd.curve_io.curveDict`.images.simMeta dictionary (e.g. {'z':'sim_redshift'} if you want to set the model
		redshift based on a simulated redshift in simMeta called 'sim_redshfit')
	guess_amplitude: bool
		If True, the amplitude parameter for the model is estimated, as well as its bounds
	trial_fit: bool
		If true, a simple minuit fit is performed to locate the parameter space for nestle fits, otherwise the full parameter
		range in bounds is used. 
	clip_data: bool
		If true, criterion like minsnr and cut_time actually will remove data from the light curve, as opposed to simply not
		fitting those data points. 
	kernel: str
		The kernel to use for microlensing GPR
	refImage: str
		The name of the image you want to be the reference image (i.e. image_1,image_2, etc.)
	nMicroSamples: int
		The number of pulls from the GPR posterior you want to use for microlensing uncertainty estimation
	color_curve: :class:`astropy.Table`
		A color curve to define the relationship between bands for parameterized light curve model.
	warning_supress: bool
		Turns on or off warnings
	verbose: bool
		Turns on/off the verbosity flag
	Returns
	-------
	fitted_curveDict: :class:`~sntd.curve_io.curveDict` or :class:`~list`
		The same curveDict that was passed to fit_data, but with new fits and time delay measurements included. List
		if list was provided.
	Examples
	--------
	>>> fitCurves=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F125W'],
		params=['x0','x1','t0','c'],constants={'z':1.33},bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},
		method='parallel',microlensing=None)

	"""

	#get together user arguments
	locs = locals()
	args = copy(locs)
	for k in kwargs.keys():
		args[k]=kwargs[k]

	if isinstance(curves,(list,tuple,np.ndarray)):

		if isinstance(curves[0],str):#then its a filename list
			filelist=True
		else:
			filelist=False
			args['curves']=[]
			for i in range(len(curves)):
				temp=_sntd_deepcopy(curves[i])
				temp.nsn=i+1
				args['curves'].append(temp)
		args['parlist']=True
	else:
		args['curves']=_sntd_deepcopy(curves)
		args['parlist']=False

	if method !='color' or identify_micro:
		args['bands'] = [bands] if bands is not None and not isinstance(bands,(tuple,list,np.ndarray)) else bands
		args['bands'] = list(set(bands)) if bands is not None else None
		#sets the bands to user's if defined (set, so that they're unique), otherwise to all the bands that exist in curves
		if args['bands']is None:
			args['bands'] = list(curves.bands) if not isinstance(curves,(list,tuple,np.ndarray)) and not isinstance(args['curves'][0],str) else None
	

	models=[models] if models is not None and not isinstance(models,(tuple,list,np.ndarray)) else models
	if models is None:
		mod,types=np.loadtxt(os.path.join(__filedir__,'data','sncosmo','models.ref'),dtype='str',unpack=True)
		modDict={mod[i]:types[i] for i in range(len(mod))}
		if isinstance(snType,str):
			if snType!='Ia':
				mods = [x[0] for x in sncosmo.models._SOURCES._loaders.keys() if x[0] in modDict.keys() and modDict[x[0]][:len(snType)]==snType]
			elif snType=='Ia':
				mods = [x[0] for x in sncosmo.models._SOURCES._loaders.keys() if 'salt2' in x[0]]
		else:
			mods=[]
			for t in snType:
				if t!='Ia':
					mods = np.append(mods,[x[0] for x in sncosmo.models._SOURCES._loaders.keys() if x[0] in modDict.keys() and modDict[x[0]][:len(t)]==t])
				elif t=='Ia':
					mods = np.append(mods,[x[0] for x in sncosmo.models._SOURCES._loaders.keys() if 'salt2' in x[0]])


	else:
		mods=models
	mods=np.unique(mods)
	for ig_mod in ignore_models:
		if ig_mod not in mods:
			temp=snana_to_sncosmo(ig_mod)
			if temp is not None:
				mods=[x for x in mods if x!=temp[1]]
		else:
			mods=[x for x in mods if x!=ig_mod]
	args['models']=mods

	if warning_supress:
		warnings.simplefilter('ignore')
	if identify_micro and not args['parlist']:
		all_bands,color_bands=identify_micro_func(args)
		args['color_bands']=color_bands
		args['bands']=all_bands
		args['curves'].micro_bands=all_bands
		args['curves'].micro_color_bands=color_bands

	if fit_prior is False:
		args['fit_prior']=None
	
	if args['parlist'] and n_per_node is None and par_or_batch=='batch':
		if nbatch_jobs is None:
			print('Must set n_per_node node and/or nbatch_jobs')
		n_per_node = math.ceil(len(args['curves'])/nbatch_jobs)



	if isinstance(method,(list,np.ndarray,tuple)):
		if len(method)==1:
			method=method[0]
		elif 'parallel' in method and fit_prior==True:
			method=np.append(['parallel'],[x for x in method if x!='parallel'])
		if args['parlist']:
			if par_or_batch=='parallel':
				print('Have not yet set up parallelized multi-fit processing')
				sys.exit(1)
			else:
				
				if n_cores_per_node>1:
					parallelize=n_cores_per_node
					n_per_node=max(n_per_node,n_cores_per_node)
					micro_par=None
				elif microlensing is not None:
					parallelize=None
					micro_par=npar_cores
				else:
					parallelize=None
					micro_par=None
				total_jobs=math.ceil(len(args['curves'])/n_per_node)
				if nbatch_jobs is None:
					nbatch_jobs=min(total_jobs,max_batch_jobs)
				script_name_init,folder_name=make_sbatch(partition=batch_partition,
												   njobs=nbatch_jobs,njobstotal=min(total_jobs,max_batch_jobs),python_path=batch_python_path,init=True,parallelize=parallelize,microlensing_cores=micro_par)
				script_name,folder_name=make_sbatch(partition=batch_partition,folder=folder_name,
												  njobs=nbatch_jobs,python_path=batch_python_path,init=False,parallelize=parallelize,microlensing_cores=micro_par)

				pickle.dump(constants,open(os.path.join(folder_name,'sntd_constants.pkl'),'wb'))
				pickle.dump(args['curves'],open(os.path.join(folder_name,'sntd_data.pkl'),'wb'))
				pyfiles=['run_sntd_init.py','run_sntd.py'] if parallelize is None else ['run_sntd_init_par.py','run_sntd_par.py']
				for pyfile in pyfiles:
					with open(os.path.join(__filedir__,'batch',pyfile)) as f:
						batch_py=f.read()
					if 'init' in pyfile:
						batch_py=batch_py.replace('nlcsreplace',str(min(int(n_per_node*nbatch_jobs),len(args['curves']))))
						batch_py=batch_py.replace('njobsreplace',str(nbatch_jobs))
					else:
						batch_py=batch_py.replace('nlcsreplace',str(n_per_node))
					if batch_init is None:
						batch_py=batch_py.replace('batchinitreplace','print("Nothing to initialize...")')
					else:
						batch_py=batch_py.replace('batchinitreplace',batch_init)
					batch_py=batch_py.replace('ncores',str(n_cores_per_node))

					indent1=batch_py.find('fitCurves=')
					indent=batch_py.find('try:')+len('try:')+1


					sntd_command=''
					for i in range(len(method)):
						fit_method=method[i]
						sntd_command+='sntd.fit_data('
						for par,val in locs.items():
							if par =='curves':
								if i==0:
									if parallelize is None:
										sntd_command+='curves=all_dat[i],'
									else:
										sntd_command+='curves=all_input,'
								else:
									sntd_command+='curves=fitCurves,'
							elif par=='constants':
								if parallelize is None:
									sntd_command+='constants=all_dat[i].constants,'
								else:
									sntd_command+='constants={'+'},'
							elif par=='batch_init':
								sntd_command+='batch_init=None,'

							elif par=='identify_micro' and identify_micro:
								if i>0:
									sntd_command+='identify_micro=False,'
								else:
									sntd_command+='identify_micro=True,'
							elif par=='bands' and identify_micro:
								if i>0:
									if parallelize is None:
										if fit_method!='color':
											sntd_command+='bands=fitCurves.micro_bands,'
										else:
											sntd_command+='bands=fitCurves.micro_color_bands,'
									else:
										print('Have not implemented this yet.')
										sys.exit(1)

								else:
									sntd_command+='bands=None,'
							elif fit_method=='color' and par=='bands':
								if color_bands is not None:
									sntd_command+='bands='+str(color_bands)+','

								else:
									sntd_command+='bands='+str(val)+','

							elif par=='method':
								sntd_command+='method="'+fit_method+'",'
							elif par=='fit_prior' and fit_method!='parallel' and (fit_prior is not None and fit_prior is not False):
								if parallelize is None:
									sntd_command+='fit_prior=fitCurves,'
								else:
									sntd_command+='fit_prior=True,'
							elif par=='par_or_batch' and parallelize is not None:
								sntd_command+='par_or_batch="parallel",'
							elif par=='npar_cores' and parallelize is not None:
								sntd_command+='npar_cores=%i,'%n_cores_per_node
							elif isinstance(val,str):
								sntd_command+=str(par)+'="'+str(val)+'",'
							elif par=='kwargs':

								for par2,val2 in val.items():
									if isinstance(val,str):
										sntd_command+=str(par2)+'="'+str(val2)+'",'
									else:
										sntd_command+=str(par2)+'='+str(val2)+','
							else:
								sntd_command+=str(par)+'='+str(val)+','

						sntd_command=sntd_command[:-1]+')\n'
						if i<len(method)-1:
							sntd_command+=' '*(indent1-indent)+'fitCurves='
				



					batch_py=batch_py.replace('sntdcommandreplace',sntd_command)

					with open(os.path.join(os.path.abspath(folder_name),pyfile),'w') as f:
						f.write(batch_py)

				return run_sbatch(folder_name,script_name_init,script_name,total_jobs,max_batch_jobs,n_per_node,wait_for_batch,parallelize,len(args['curves']))


		else:
			initBounds=deepcopy(args['bounds'])
			if 'parallel' in method:
				if verbose:
					print('Starting parallel method...')
				curves=_fitparallel(args)
				if args['fit_prior']==True:
					args['fit_prior']=curves
				args['curves']=curves
				args['bounds']=copy(initBounds)
			if 'series' in method:
				if verbose:
					print('Starting series method...')
				if 'td' not in args['bounds']:
					if verbose:
						print('td not in bounds for series method, choosing based on parallel bounds...')
					args['bounds']['td']=args['bounds']['t0']
				if 'mu' not in args['bounds']:
					if verbose:
						print('mu not in bounds for series method, choosing defaults...')
					args['bounds']['mu']=[0,10]

				curves=_fitseries(args)
				args['curves']=curves
				args['bounds']=copy(initBounds)
			if 'color' in method:
				if verbose:
					print('Starting color method...')
				if 'td' not in args['bounds']:
					if verbose:
						print('td not in bounds for color method, choosing based on parallel bounds...')
					args['bounds']['td']=args['bounds']['t0']
				
				curves=_fitColor(args)
	elif method not in ['parallel','series','color']:
			raise RuntimeError('Parameter "method" must be "parallel","series", or "color".')

	elif method=='parallel':
		if args['parlist']:
			if par_or_batch=='parallel':
				par_arg_vals=[]
				for i in range(len(args['curves'])):
					temp_args={}
					
					for par_key in ['snType','bounds','constants','t0_guess']:
						if isinstance(args[par_key],(list,tuple,np.ndarray)):
							try:
								temp_args[par_key]=args[par_key][i]
							except:
								pass
					for par_key in ['bands','models','ignore','params']:
						if isinstance(args[par_key],(list,tuple,np.ndarray)) and np.any([isinstance(x,(list,tuple,np.ndarray)) for x in args[par_key]]):
							try:
								temp_args[par_key]=args[par_key][i]
							except:
								pass
					par_arg_vals.append([args['curves'][i],temp_args])

				curves=pyParz.foreach(par_arg_vals,_fitparallel,[args],numThreads=min(npar_cores,len(par_arg_vals)))
			else:
				if n_cores_per_node>1:
					parallelize=n_cores_per_node
					n_per_node=max(n_per_node,n_cores_per_node)
					micro_par=None
				elif microlensing is not None:
					parallelize=None
					micro_par=npar_cores
				else:
					parallelize=None
					micro_par=None
				total_jobs=math.ceil(len(args['curves'])/n_per_node)
				if nbatch_jobs is None:
					nbatch_jobs=min(total_jobs,max_batch_jobs)
				script_name_init,folder_name=make_sbatch(partition=batch_partition,
												   njobs=nbatch_jobs,njobstotal=min(total_jobs,max_batch_jobs),python_path=batch_python_path,init=True,parallelize=parallelize,microlensing_cores=micro_par)
				script_name,folder_name=make_sbatch(partition=batch_partition,folder=folder_name,
												  njobs=nbatch_jobs,python_path=batch_python_path,init=False,parallelize=parallelize,microlensing_cores=micro_par)

				pickle.dump(constants,open(os.path.join(folder_name,'sntd_constants.pkl'),'wb'))
				pickle.dump(args['curves'],open(os.path.join(folder_name,'sntd_data.pkl'),'wb'))
				pyfiles=['run_sntd_init.py','run_sntd.py'] if parallelize is None else ['run_sntd_init_par.py','run_sntd_par.py']
				for pyfile in pyfiles:
					with open(os.path.join(__filedir__,'batch',pyfile)) as f:
						batch_py=f.read()
					if 'init' in pyfile:
						batch_py=batch_py.replace('nlcsreplace',str(min(int(n_per_node*nbatch_jobs),len(args['curves']))))
						batch_py=batch_py.replace('njobsreplace',str(nbatch_jobs))
					else:
						batch_py=batch_py.replace('nlcsreplace',str(n_per_node))
					if batch_init is None:
						batch_py=batch_py.replace('batchinitreplace','print("Nothing to initialize...")')
					else:
						batch_py=batch_py.replace('batchinitreplace',batch_init)
					batch_py=batch_py.replace('ncores',str(n_cores_per_node))

					indent1=batch_py.find('fitCurves=')
					indent=batch_py.find('try:')+len('try:')+1


					sntd_command='sntd.fit_data('
					for par,val in locs.items():

						if par =='curves':
							if parallelize is None:
								sntd_command+='curves=all_dat[i],'
							else:
								sntd_command+='curves=all_input,'
						elif par=='batch_init':
							sntd_command+='batch_init=None,'
						elif par=='constants':
							if parallelize is None:
									sntd_command+='constants=all_dat[i].constants,'
							else:
								sntd_command+='constants=const_list,'
						elif par=='method':
							sntd_command+='method="parallel",'
						elif par=='par_or_batch' and parallelize is not None:
								sntd_command+='par_or_batch="parallel",'
						elif par=='npar_cores' and parallelize is not None:
							sntd_command+='npar_cores=%i,'%n_cores_per_node
						elif isinstance(val,str):
							sntd_command+=str(par)+'="'+str(val)+'",'
						elif par=='kwargs':

							for par2,val2 in val.items():
								if isinstance(val,str):
									sntd_command+=str(par2)+'="'+str(val2)+'",'
								else:
									sntd_command+=str(par2)+'='+str(val2)+','

						else:
							sntd_command+=str(par)+'='+str(val)+','

					sntd_command=sntd_command[:-1]+')'

					batch_py=batch_py.replace('sntdcommandreplace',sntd_command)

					with open(os.path.join(os.path.abspath(folder_name),pyfile),'w') as f:
						f.write(batch_py)

				

				
				return run_sbatch(folder_name,script_name_init,script_name,total_jobs,max_batch_jobs,n_per_node,wait_for_batch,parallelize,len(args['curves']))



		else:
			curves=_fitparallel(args)
	elif method=='series':
		if args['parlist']:
			if par_or_batch=='parallel':
				par_arg_vals=[]
				for i in range(len(args['curves'])):
					temp_args={}
					try:
						for par_key in ['snType','bounds','constants','t0_guess']:
							if isinstance(args[par_key],(list,tuple,np.ndarray)):
								temp_args[par_key]=args[par_key][i]
						for par_key in ['bands','models','ignore','params']:
							if isinstance(args[par_key],(list,tuple,np.ndarray)) and np.any([isinstance(x,(list,tuple,np.ndarray)) for x in args[par_key]]):
								temp_args[par_key]=args[par_key][i]
					except:
						pass
					par_arg_vals.append([args['curves'][i],temp_args])
				curves=pyParz.foreach(par_arg_vals,_fitseries,[args],numThreads=min(npar_cores,len(par_arg_vals)))
			else:
				if n_cores_per_node>1:
					parallelize=n_cores_per_node
					n_per_node=max(n_per_node,n_cores_per_node)
					micro_par=None
				elif microlensing is not None:
					parallelize=None
					micro_par=npar_cores
				else:
					parallelize=None
					micro_par=None
				total_jobs=math.ceil(len(args['curves'])/n_per_node)
				if nbatch_jobs is None:
					nbatch_jobs=min(total_jobs,max_batch_jobs)
				script_name_init,folder_name=make_sbatch(partition=batch_partition,
												   njobs=nbatch_jobs,njobstotal=min(total_jobs,max_batch_jobs),python_path=batch_python_path,init=True,parallelize=parallelize,microlensing_cores=micro_par)
				script_name,folder_name=make_sbatch(partition=batch_partition,folder=folder_name,
												  njobs=nbatch_jobs,python_path=batch_python_path,init=False,parallelize=parallelize,microlensing_cores=micro_par)

				pickle.dump(constants,open(os.path.join(folder_name,'sntd_constants.pkl'),'wb'))
				pickle.dump(args['curves'],open(os.path.join(folder_name,'sntd_data.pkl'),'wb'))
				pyfiles=['run_sntd_init.py','run_sntd.py'] if parallelize is None else ['run_sntd_init_par.py','run_sntd_par.py']
				for pyfile in pyfiles:
					with open(os.path.join(__filedir__,'batch',pyfile)) as f:
						batch_py=f.read()
					if 'init' in pyfile:
						batch_py=batch_py.replace('nlcsreplace',str(min(int(n_per_node*nbatch_jobs),len(args['curves']))))
						batch_py=batch_py.replace('njobsreplace',str(nbatch_jobs))
					else:
						batch_py=batch_py.replace('nlcsreplace',str(n_per_node))
					if batch_init is None:
						batch_py=batch_py.replace('batchinitreplace','print("Nothing to initialize...")')
					else:
						batch_py=batch_py.replace('batchinitreplace',batch_init)
					batch_py=batch_py.replace('ncores',str(n_cores_per_node))

					indent1=batch_py.find('fitCurves=')
					indent=batch_py.find('try:')+len('try:')+1
					sntd_command='sntd.fit_data('
					for par,val in locs.items():
						if par =='curves':
							if parallelize is None:
								sntd_command+='curves=all_dat[i],'
							else:
								sntd_command+='curves=all_input,'
						elif par=='batch_init':
							sntd_command+='batch_init=None,'
						elif par=='constants':
							if parallelize is None:
									sntd_command+='constants=all_dat[i].constants,'
							else:
								sntd_command+='constants={'+'},'
						elif par=='method':
							sntd_command+='method="series",'
						elif par=='par_or_batch' and parallelize is not None:
								sntd_command+='par_or_batch="parallel",'
						elif par=='npar_cores' and parallelize is not None:
							sntd_command+='npar_cores=%i,'%n_cores_per_node
						elif isinstance(val,str):
							sntd_command+=str(par)+'="'+str(val)+'",'
						elif par=='kwargs':

							for par2,val2 in val.items():
								if isinstance(val,str):
									sntd_command+=str(par2)+'="'+str(val2)+'",'
								else:
									sntd_command+=str(par2)+'='+str(val2)+','
						else:
							sntd_command+=str(par)+'='+str(val)+','

					sntd_command=sntd_command[:-1]+')'

					batch_py=batch_py.replace('sntdcommandreplace',sntd_command)

					with open(os.path.join(os.path.abspath(folder_name),pyfile),'w') as f:
						f.write(batch_py)
				
				return run_sbatch(folder_name,script_name_init,script_name,total_jobs,max_batch_jobs,n_per_node,wait_for_batch,parallelize,len(args['curves']))
		else:
			curves=_fitseries(args)

	elif method=='color':
		if args['parlist']:
			if par_or_batch=='parallel':
				par_arg_vals=[]
				for i in range(len(args['curves'])):
					temp_args={}
					try:
						for par_key in ['snType','bounds','constants','t0_guess']:
							if isinstance(args[par_key],(list,tuple,np.ndarray)):
								temp_args[par_key]=args[par_key][i]
						for par_key in ['bands','models','ignore','params']:
							if isinstance(args[par_key],(list,tuple,np.ndarray)) and np.any([isinstance(x,(list,tuple,np.ndarray)) for x in args[par_key]]):
								temp_args[par_key]=args[par_key][i]
					except:
						pass
					par_arg_vals.append([args['curves'][i],temp_args])
				curves=pyParz.foreach(par_arg_vals,_fitColor,[args],numThreads=min(npar_cores,len(par_arg_vals)))
			else:
				if n_cores_per_node>1:
					parallelize=n_cores_per_node
					n_per_node=max(n_per_node,n_cores_per_node)
					micro_par=None
				elif microlensing is not None:
					parallelize=None
					micro_par=npar_cores
				else:
					parallelize=None
					micro_par=None
				total_jobs=math.ceil(len(args['curves'])/n_per_node)
				if nbatch_jobs is None:
					nbatch_jobs=min(total_jobs,max_batch_jobs)
				script_name_init,folder_name=make_sbatch(partition=batch_partition,
												   njobs=nbatch_jobs,njobstotal=min(total_jobs,max_batch_jobs),python_path=batch_python_path,init=True,parallelize=parallelize,microlensing_cores=micro_par)
				script_name,folder_name=make_sbatch(partition=batch_partition,folder=folder_name,
												  njobs=nbatch_jobs,python_path=batch_python_path,init=False,parallelize=parallelize,microlensing_cores=micro_par)

				pickle.dump(constants,open(os.path.join(folder_name,'sntd_constants.pkl'),'wb'))
				pickle.dump(args['curves'],open(os.path.join(folder_name,'sntd_data.pkl'),'wb'))
				pyfiles=['run_sntd_init.py','run_sntd.py'] if parallelize is None else ['run_sntd_init_par.py','run_sntd_par.py']
				for pyfile in pyfiles:
					with open(os.path.join(__filedir__,'batch',pyfile)) as f:
						batch_py=f.read()
					if 'init' in pyfile:
						batch_py=batch_py.replace('nlcsreplace',str(min(int(n_per_node*nbatch_jobs),len(args['curves']))))
						batch_py=batch_py.replace('njobsreplace',str(nbatch_jobs))
					else:
						batch_py=batch_py.replace('nlcsreplace',str(n_per_node))
					if batch_init is None:
						batch_py=batch_py.replace('batchinitreplace','print("Nothing to initialize...")')
					else:
						batch_py=batch_py.replace('batchinitreplace',batch_init)
					batch_py=batch_py.replace('ncores',str(n_cores_per_node))

					indent1=batch_py.find('fitCurves=')
					indent=batch_py.find('try:')+len('try:')+1
					sntd_command='sntd.fit_data('
					for par,val in locs.items():
						if par =='curves':
							if parallelize is None:
								sntd_command+='curves=all_dat[i],'
							else:
								sntd_command+='curves=all_input,'
						elif par=='batch_init':
							sntd_command+='batch_init=None,'
						elif par=='constants':
							if parallelize is None:
									sntd_command+='constants=all_dat[i].constants,'
							else:
								sntd_command+='constants={'+'},'
						elif par=='method':
							sntd_command+='method="color",'
						elif par=='par_or_batch' and parallelize is not None:
								sntd_command+='par_or_batch="parallel",'
						elif par=='npar_cores' and parallelize is not None:
							sntd_command+='npar_cores=%i,'%n_cores_per_node
						elif isinstance(val,str):
							sntd_command+=str(par)+'="'+str(val)+'",'
						elif par=='kwargs':

							for par2,val2 in val.items():
								if isinstance(val,str):
									sntd_command+=str(par2)+'="'+str(val2)+'",'
								else:
									sntd_command+=str(par2)+'='+str(val2)+','
						else:
							sntd_command+=str(par)+'='+str(val)+','

					sntd_command=sntd_command[:-1]+')'

					batch_py=batch_py.replace('sntdcommandreplace',sntd_command)

					with open(os.path.join(os.path.abspath(folder_name),pyfile),'w') as f:
						f.write(batch_py)

				

				return run_sbatch(folder_name,script_name_init,script_name,total_jobs,max_batch_jobs,n_per_node,wait_for_batch,parallelize,len(args['curves']))
		else:
			if args['color_bands'] is not None:
				args['bands']=args['color_bands']

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
		if isinstance(curves,str):
			args['curves']=pickle.load(open(curves,'rb'))
		else:
			args['curves']=curves

			if args['verbose']:
				print('Fitting MISN number %i...'%curves.nsn)
	else:
		args=all_args

	for p in args['curves'].constants.keys():
		if p not in args['constants'].keys():
			args['constants'][p]=args['curves'].constants[p]

	if args['clip_data']:
		for im in args['curves'].images.keys():
			args['curves'].clip_data(im=im,minsnr=args.get('minsnr',0),max_cadence=args['max_cadence'])

	args['bands']=list(args['bands'])
	_,band_SNR,_=getBandSNR(args['curves'],args['bands'],args['min_points_per_band'])
	
	if len(args['bands'])<2:
		raise RuntimeError("If you want to analyze color curves, you need two bands!")
	else:
		if args['fit_colors'] is None:
			final_bands=[]
			for band in np.unique(args['curves'].images[args['refImage']].table['band']):
				to_add=True
				for im in args['curves'].images.keys():
					if len(np.where(args['curves'].images[im].table['band']==band)[0])<args['min_points_per_band']:
						to_add=False
				if to_add:
					final_bands.append(band)
			if np.any([x not in final_bands for x in args['bands']]):
				all_SNR=[]
				for band in final_bands:
					ims=[]
					for d in args['curves'].images.keys():
						inds=np.where(args['curves'].images[d].table['band']==band)[0]
						if len(inds)==0:
							ims.append(0)
						else:
							ims.append(np.sum(args['curves'].images[d].table['flux'][inds]/args['curves'].images[d].table['fluxerr'][inds])*\
								 np.sqrt(len(inds)))
					all_SNR.append(np.sum(ims))
				sorted=np.flip(np.argsort(all_SNR))
				args['bands']=np.array(final_bands)[sorted]
				if args['max_n_bands'] is not None:
					args['bands']=args['bands'][:args['max_n_bands']]
			colors_to_fit=[x for x in combinations(args['bands'],2)]
			if args['color_bands'] is not None:
				for i in range(len(colors_to_fit)):
					colors_to_fit[i]=[x for x in args['color_bands'] if x in colors_to_fit[i]]

		else:
			colors_to_fit=[x.split('-') for x in args['fit_colors']]


	imnums=[x[-1] for x in args['curves'].images.keys()]
	if args['fit_prior'] is not None:
		if args['fit_prior']==True:
			args['fit_prior']=args['curves']
		ref=args['fit_prior'].parallel.fitOrder[0]
		refnum=ref[-1]
	else:
		ref=args['refImage']
		refnum=ref[-1]
	inds=np.arange(0,len(args['curves'].images[ref].table),1).astype(int)
	nimage=len(imnums)
	snParams=['dt_%s'%i for i in imnums if i!=refnum]
	all_vparam_names=np.append(args['params'],
							   snParams).flatten()
	
	ims=list(args['curves'].images.keys())

	
	for param in all_vparam_names:
		if param in args['color_param_ignore'] and args['fit_prior'] is not None and param not in args['constants']:
			par_ref=args['fit_prior'].parallel.fitOrder[0]
			args['constants'][param]=args['fit_prior'].images[par_ref].param_quantiles[param][1]
			if param in all_vparam_names:
				all_vparam_names=np.array([x for x in all_vparam_names if x !=param])
		if param not in args['bounds'].keys():
			if param.startswith('dt_'):
				if args['fit_prior'] is not None:
					im=[x for x in ims if x[-1]==param[-1]][0]
					args['bounds'][param]=np.array([-1,1])*3*np.sqrt(args['fit_prior'].parallel.time_delay_errors[im]**2+ \
													args['fit_prior'].parallel.time_delay_errors[ref]**2)+\
							(args['fit_prior'].parallel.time_delays[im]-args['fit_prior'].parallel.time_delays[ref])
					
				else:
					args['bounds'][param]=np.array(args['bounds']['td'])#+time_delays[im]
			
			

		elif args['fit_prior'] is not None:
			par_ref=args['fit_prior'].parallel.fitOrder[0]
			if param not in args['fit_prior'].images[par_ref].param_quantiles.keys():
				continue
			args['bounds'][param]=3*np.array([args['fit_prior'].images[par_ref].param_quantiles[param][0]- \
											  args['fit_prior'].images[par_ref].param_quantiles[param][1],
											  args['fit_prior'].images[par_ref].param_quantiles[param][2]- \
											  args['fit_prior'].images[par_ref].param_quantiles[param][1]])+ \
								  args['fit_prior'].images[par_ref].param_quantiles[param][1]



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

	if 'ignore_models' in args['set_from_simMeta'].keys():
		to_ignore=args['curves'].images[ref].simMeta[args['set_from_simMeta']['ignore_models']]
		if isinstance(to_ignore,str):
			to_ignore=[to_ignore]
		args['models']=[x for x in np.array(args['models']).flatten() if x not in to_ignore]

	if args['fit_prior'] is not None and args['fit_prior'].images[args['fit_prior'].parallel.fitOrder[0]].fits.model._source.name not in args['models']:
		print('Wanted to use a fit prior but do not have the same model as an option.')
		raise RuntimeError
	elif args['fit_prior'] is not None:
		args['models']=args['fit_prior'].images[args['fit_prior'].parallel.fitOrder[0]].fits.model._source.name

	if not args['curves'].quality_check(min_n_bands=2,
						min_n_points_per_band=args['min_points_per_band'],clip=False,method='parallel'):
		return
	all_fit_dict={}
	if args['fast_model_selection'] and len(np.array(args['models']).flatten())>1:
		for b in args['force_positive_param']:
			if b in args['bounds'].keys():
				args['bounds'][b]=np.array([max([args['bounds'][b][0],0]),max([args['bounds'][b][1],0])])
			else:
				args['bounds'][b]=np.array([0,np.inf])
		
		minchisq=np.inf
		init_inds=deepcopy(inds)
		for mod in np.array(args['models']).flatten():
			inds=deepcopy(init_inds)
			if isinstance(mod,str):
				if mod.upper() in ['BAZIN','BAZINSOURCE']:
					mod='BAZINSOURCE'
					if len(np.unique(args['curves'].images[ref].table['band']))>1 and args['color_curve'] is None:
						best_band=band_SNR[args['fitOrder'][0]][0]
						inds=np.where(args['curves'].images[ref].table['band']==best_band)[0]
						

					source=BazinSource(data=args['curves'].images[ref].table[inds],colorCurve=args['color_curve'])
				else:
					source=sncosmo.get_source(mod)
				tempMod = sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
			else:
				tempMod=copy(mod)
			
			tempMod.set(**{k:args['constants'][k] for k in args['constants'].keys() if k in tempMod.param_names})
			tempMod.set(**{k:args['curves'].images[args['refImage']].simMeta[args['set_from_simMeta'][k]] for k in args['set_from_simMeta'].keys() if k in tempMod.param_names})

			if mod=='BAZINSOURCE':
				tempMod.set(z=0)
			try:
				res,fit=sncosmo.fit_lc(args['curves'].images[ref].table[inds],tempMod,[x for x in args['params'] if x in tempMod.param_names],
										bounds={b:args['bounds'][b] for b in args['bounds'] if b not in ['t0',tempMod.param_names[2]]},
										minsnr=args.get('minsnr',0))
			except:
				if args['verbose']:
					print('Issue with %s, skipping...'%mod)
				continue
			tempchisq=res.chisq/(len(inds)+len([x for x in args['params'] if x in tempMod.param_names])-1)
			if tempchisq<minchisq:
				minchisq=tempchisq
				bestres=copy(res)
				bestfit=copy(fit)
				bestmodname=copy(mod)
			all_fit_dict[mod]=[copy(fit),copy(res)]
		try:
			args['models']=[bestmodname]
		except:
			print('Every model had an error.')
			sys.exit(1)
	finallogz=-np.inf
	for mod in np.array(args['models']).flatten():


		if isinstance(mod,str):
			source=sncosmo.get_source(mod)
			tempMod = sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
		else:
			tempMod=copy(mod)
		tempMod.set(**{k:args['constants'][k] for k in args['constants'].keys() if k in tempMod.param_names})
		tempMod.set(**{k:args['curves'].images[args['refImage']].simMeta[args['set_from_simMeta'][k]] for k in args['set_from_simMeta'].keys() if k in tempMod.param_names})

		if args['fit_prior'] is not None:
			par_ref=args['fit_prior'].parallel.fitOrder[0]

			if mod!=args['fit_prior'].images[par_ref].fits.model._source.name:
				continue
			temp_delays={k:args['fit_prior'].parallel.time_delays[k]-args['fit_prior'].parallel.time_delays[par_ref]\
						 for k in args['fit_prior'].parallel.fitOrder}
			
			args['curves'].color_table([x[0] for x in colors_to_fit],[x[1] for x in colors_to_fit],time_delays={im:0 for im in args['curves'].images.keys()},
										minsnr=args.get('minsnr',0))
			args['curves'].color.meta['reft0']=args['fit_prior'].images[par_ref].fits.model.get('t0')
			args['curves'].color.meta['td']=temp_delays
		else:
			par_ref=args['refImage']
			im_name=args['refImage'][:-1]
			if args['trial_fit']:
				for b in args['force_positive_param']:
					if b in args['bounds'].keys():
						args['bounds'][b]=np.array([max([args['bounds'][b][0],0]),max([args['bounds'][b][1],0])])
					else:
						args['bounds'][b]=np.array([0,np.inf])
				best_bands=band_SNR[args['refImage']][:min(len(band_SNR[args['refImage']]),2)]
				temp_delays={}
				temp_mags={}
				fit_order=np.flip(args['fitOrder']) if args['fitOrder'] is not None else \
										[x for x in args['curves'].images.keys() if x != args['refImage']]+[args['refImage']]
				for im in fit_order:
					temp_bands=[]
					for b in best_bands:
						temp_bands=np.append(temp_bands,np.where(args['curves'].images[im].table['band']==b)[0])
					temp_inds=temp_bands.astype(int)
					
					res,fit=sncosmo.fit_lc(deepcopy(args['curves'].images[im].table[temp_inds]),tempMod,
											[x for x in args['params'] if x in tempMod.param_names]+[tempMod.param_names[2]]+\
											[x for x in tempMod.param_names if x in args['bounds'].keys()],
											bounds={b:args['bounds'][b] for b in args['bounds'].keys() if b not in ['t0',tempMod.param_names[2]]},
											minsnr=args.get('minsnr',0))
					temp_delays[im]=fit.get('t0')
				for param in args['color_param_ignore']:
					if param not in args['constants']:
						args['constants'][param]=fit.get(param)
					if param in all_vparam_names:
						all_vparam_names=np.array([x for x in all_vparam_names if x !=param])

				tempMod.set(**{k:args['constants'][k] for k in args['constants'].keys() if k in tempMod.param_names})
				args['curves'].color.meta['reft0']=temp_delays[args['refImage']]

				
				temp_delays={im:temp_delays[im]-temp_delays[args['refImage']] for im in temp_delays.keys()}
				for b in args['bounds']:
					if b in list(res.errors.keys()):
						if b not in all_vparam_names:
							tempMod.set(**{b:fit.get(b)})
						elif b !='t0':
							args['bounds'][b]=(np.array(args['bounds'][b])-np.median(args['bounds'][b]))/2+fit.get(b)
						else:
							args['bounds'][b]=(np.array(args['bounds'][b])-np.median(args['bounds'][b]))/2+args['curves'].color.meta['reft0']
					 

					elif b.startswith('dt_'):
						args['bounds'][b]=np.array(args['bounds']['td'])/2+temp_delays[im_name+b[-1]]
					
				
				if 't0' not in args['bounds'].keys():
					args['bounds']['t0']=np.array(args['bounds']['td'])/2+args['curves'].color.meta['reft0']

				args['curves'].color_table([x[0] for x in colors_to_fit],[x[1] for x in colors_to_fit],time_delays={im:0 for im in args['curves'].images.keys()},
												minsnr=args.get('minsnr',0))
				args['curves'].color.meta['td']=temp_delays
				
			else:
				args['curves'].color_table([x[0] for x in colors_to_fit],[x[1] for x in colors_to_fit],referenceImage=args['refImage'],static=True,model=tempMod,
																	minsnr=args.get('minsnr',0))
				for b in args['bounds']:
					if b.startswith('dt_'):
						args['bounds'][b]=np.array(args['bounds']['td'])+args['curves'].color.meta['td'][im_name+b[-1]]
					elif b=='t0':
						args['bounds'][b]=np.array(args['bounds'][b])+args['curves'].color.meta['reft0']

				if 't0' not in args['bounds'].keys():
					args['bounds']['t0']=np.array(args['bounds']['td'])+args['curves'].color.meta['reft0']


		

		

		if args['cut_time'] is not None:
			
			for im in args['curves'].images.keys():
				args['curves'].color.table=args['curves'].color.table[np.where(np.logical_or(args['curves'].color.table['image']!=im,
																	np.logical_and(args['curves'].color.table['time']>= \
																	args['cut_time'][0]*(1+tempMod.get('z'))+args['curves'].color.meta['reft0']+\
																	args['curves'].color.meta['td'][im],
																	args['curves'].color.table['time']<= \
																	args['cut_time'][1]*(1+tempMod.get('z'))+args['curves'].color.meta['reft0']+\
																	args['curves'].color.meta['td'][im])))[0]]
		
		all_vparam_names=np.array([x for x in all_vparam_names if x!=tempMod.param_names[2]])
		if args['band_order'] is not None:
			args['bands']=[x for x in args['band_order'] if x in args['bands']]
		for b in args['force_positive_param']:
			if b in args['bounds'].keys():
				args['bounds'][b]=np.array([max([args['bounds'][b][0],0]),max([args['bounds'][b][1],0])])
			else:
				args['bounds'][b]=np.array([0,np.inf])
		if not args['curves'].quality_check(min_n_bands=args['min_n_bands'],
						min_n_points_per_band=args['min_points_per_band'],clip=args['clip_data'],method='color'):
			return
		params,res,model=nest_color_lc(args['curves'].color.table,tempMod,nimage,colors=colors_to_fit,
											bounds=args['bounds'],
											 vparam_names=[x for x in all_vparam_names if x in tempMod.param_names or x in snParams],ref=par_ref,
											 minsnr=args.get('minsnr',5.),priors=args.get('priors',None),ppfs=args.get('ppfs',None),
											 method=args.get('nest_method','single'),maxcall=args.get('maxcall',None),
											 modelcov=args.get('modelcov',None),rstate=args.get('rstate',None),
											 maxiter=args.get('maxiter',None),npoints=args.get('npoints',100))
		if finallogz<res.logz:
			finallogz=res.logz
			finalres,finalmodel=res,model
			time_delays=args['curves'].color.meta['td']
			final_param_quantiles=params

	args['curves'].color.time_delays=dict([])
	args['curves'].color.time_delay_errors=dict([])

	args['curves'].color.t_peaks=dict([])

	
	if par_ref==args['refImage']:
		args['curves'].color.time_delays[par_ref]=0
		args['curves'].color.time_delay_errors[par_ref]=np.array([0,0])
		args['curves'].color.t_peaks[par_ref]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('t0')],.5,finalres.weights)
		for k in args['curves'].images.keys():
			if k==par_ref:
				continue
			else:
				args['curves'].color.t_peaks[k]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('dt_'+k[-1])]+\
												finalres.samples[:,finalres.vparam_names.index('t0')],
												.5,finalres.weights)
				dt_quant=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('dt_'+k[-1])],[.16,.5,.84],finalres.weights)
				args['curves'].color.time_delays[k]=dt_quant[1]
				args['curves'].color.time_delay_errors[k]=np.array([dt_quant[0]-dt_quant[1],dt_quant[2]-dt_quant[1]])
				
	else:
		args['curves'].color.time_delays[args['refImage']]=0
		args['curves'].color.time_delay_errors[args['refImage']]=np.array([0,0])
		
		trefSamples=finalres.samples[:,finalres.vparam_names.index('dt_'+args['refImage'][-1])]
		args['curves'].color.t_peaks[args['refImage']]=weighted_quantile(trefSamples+finalres.samples[:,finalres.vparam_names.index('t0')],.5,finalres.weights)
		for k in args['curves'].images.keys():
			if k==args['refImage']:
				continue
			elif k==par_ref:
				args['curves'].color.t_peaks[k]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('t0')],.5,finalres.weights)
				
				dt_quant=weighted_quantile(-1*trefSamples,[.16,.5,.84],finalres.weights)
				args['curves'].color.time_delays[k]=dt_quant[1]
				args['curves'].color.time_delay_errors[k]=np.array([dt_quant[0]-dt_quant[1],dt_quant[2]-dt_quant[1]])
				
			else:
				args['curves'].color.t_peaks[k]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('dt_'+k[-1])]+\
												finalres.samples[:,finalres.vparam_names.index('t0')],
												.5,finalres.weights)
				
				dt_quant=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('dt_'+k[-1])]-trefSamples,[.16,.5,.84],finalres.weights)
				args['curves'].color.time_delays[k]=dt_quant[1]
				args['curves'].color.time_delay_errors[k]=np.array([dt_quant[0]-dt_quant[1],dt_quant[2]-dt_quant[1]])
				
		

	finalmodel.set(t0=args['curves'].color.t_peaks[args['refImage']])

	args['curves'].color_table([x[0] for x in colors_to_fit],[x[1] for x in colors_to_fit],time_delays=args['curves'].color.time_delays,minsnr=args.get('minsnr',0))
	args['curves'].color.meta['td']=time_delays
	args['curves'].color.meta['fit_colors']=colors_to_fit
	args['curves'].color.refImage=args['refImage']
	args['curves'].color.priorImage=par_ref
	args['curves'].color.bands=args['bands']

	args['curves'].color.fits=newDict()
	args['curves'].color.fits['model']=finalmodel
	args['curves'].color.fits['res']=finalres

	return args['curves']

def nest_color_lc(data,model,nimage,colors, vparam_names,bounds,ref='image_1',
				   minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
				   maxiter=None, maxcall=None, modelcov=False, rstate=None,
				   verbose=False, warn=True,**kwargs):
	####Taken from SNCosmo nest_lc

	# experimental parameters
	tied = kwargs.get("tied", None)

	vparam_names=list(vparam_names)
	if ppfs is None:
		ppfs = {}
	if tied is None:
		tied = {}


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


	model_param_names=[x for x in vparam_names[:len(vparam_names)-(nimage-1)]]

	model_idx = np.array([vparam_names.index(name) for name in model_param_names])
	td_params=[x for x in vparam_names[len(vparam_names)-nimage:] if x.startswith('dt')]
	td_idx=np.array([vparam_names.index(name) for name in td_params])



	im_indices=[np.where(data['image']==i)[0] for i in np.unique(data['image']) if i !=ref]
	
	nonan_dict={band:np.where(~np.isnan(data['flux_%s'%band]))[0] for band in np.unique(np.array(colors).flatten())}
	color_ind_dict={color[0]+'-'+color[1]:np.where(~np.isnan(data[color[0]+'-'+color[1]]))[0] for color in colors}

	obs_dict={}
	err_dict={}
	zp_dict={}
	for color in colors:
		col_inds=color_ind_dict[color[0]+'-'+color[1]]
		obs_dict[color[0]+'-'+color[1]]=data['flux_%s'%color[0]][col_inds]/data['flux_%s'%color[1]][col_inds]
		err_dict[color[0]+'-'+color[1]]=(data['flux_%s'%color[0]][col_inds]/data['flux_%s'%color[1]][col_inds])*\
							np.sqrt((data['fluxerr_%s'%color[0]][col_inds]/data['flux_%s'%color[0]][col_inds])**2+\
																(data['fluxerr_%s'%color[1]][col_inds]/data['flux_%s'%color[1]][col_inds])**2)
	

	unique_bands=np.unique(np.array(colors).flatten())
	zp_dict={b:data['zp_%s'%b][nonan_dict[b][0]] for b in unique_bands}
	zpsys=data['zpsys'][0]

	
	def chisq_likelihood(parameters):
		model.set(**{model_param_names[k]:parameters[model_idx[k]] for k in range(len(model_idx))})
		all_data=deepcopy(data)

		for i in range(len(im_indices)):
			all_data['time'][im_indices[i]]-=parameters[td_idx[i]]
		
		
		mod_dict={}
		cov_dict={}
		for b in unique_bands:
			time=all_data[nonan_dict[b]]['time']
			mod_dict[b]=model.bandflux(b,time,zpsys=zpsys,zp=zp_dict[b])
			if modelcov:
				

				_, mcov = model.bandfluxcov(b,
											time,
										zp=zp_dict[b],
										zpsys=zpsys)

				
				cov_dict[b]=mcov

		chisq=0
		for color in colors:
			col_inds=color_ind_dict[color[0]+'-'+color[1]]
			mod_flux1=mod_dict[color[0]]
			mod_flux2=mod_dict[color[1]]
			color_inds1=[i for i in range(len(nonan_dict[color[0]])) if nonan_dict[color[0]][i] in col_inds]
			color_inds2=[i for i in range(len(nonan_dict[color[1]])) if nonan_dict[color[1]][i] in col_inds]
			model_observations=mod_flux1[color_inds1]/mod_flux2[color_inds2]
			

			
			obs=obs_dict[color[0]+'-'+color[1]]
			err=err_dict[color[0]+'-'+color[1]]
			
			colzp1=data['zp_'+color[0]][col_inds[0]]
			colzp2=data['zp_'+color[1]][col_inds[0]]
			
			good_obs=np.where(~np.isnan(model_observations))[0]
			model_observations=model_observations[good_obs]
			obs=obs[good_obs]
			err=err[good_obs]
			
			if modelcov:
				cov=np.diag(err)
				
				mcov1=cov_dict[color[0]][:,np.array(color_inds1)[good_obs]]
				mcov1=mcov1[np.array(color_inds1)[good_obs],:]
				mcov2=cov_dict[color[1]][:,np.array(color_inds2)[good_obs]]
				mcov2=mcov2[np.array(color_inds2)[good_obs],:]

				
				cov = cov + np.sqrt(mcov1**2+mcov2**2)
				invcov = np.linalg.pinv(cov)
				diff = obs-model_observations
				chisq+=np.dot(np.dot(diff, invcov), diff)

			else:
				chisq+=np.sum((obs-model_observations)**2/\
							 err**2)

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

	return params,res,model


def _fitseries(all_args):
	if isinstance(all_args,(list,tuple,np.ndarray)):
		curves,args=all_args
		if isinstance(args,list):
			args=args[0]
		if isinstance(curves,list):
			curves,single_par_vars=curves
			for key in single_par_vars:
				args[key]=single_par_vars[key]

		if isinstance(curves,str):
			args['curves']=pickle.load(open(curves,'rb'))
		else:
			args['curves']=curves
			if args['verbose']:
				print('Fitting MISN number %i...'%curves.nsn)
	else:
		args=all_args


	for p in args['curves'].constants.keys():
		if p not in args['constants'].keys():
			args['constants'][p]=args['curves'].constants[p]

	if args['clip_data']:
		for im in args['curves'].images.keys():
			args['curves'].clip_data(im=im,minsnr=args.get('minsnr',0),max_cadence=args['max_cadence'])

	args['bands'],band_SNR,_=getBandSNR(args['curves'],args['bands'],args['min_points_per_band'])
	args['curves'].series.bands=args['bands'][:args['max_n_bands']]if args['max_n_bands'] is not None else args['bands']

	

	imnums=[x[-1] for x in args['curves'].images.keys()]
	if args['fit_prior'] is not None:
		if args['fit_prior']==True:
			args['fit_prior']=args['curves']
		ref=args['fit_prior'].parallel.fitOrder[0]
		refnum=ref[-1]
	else:
		ref=args['refImage']
		refnum=ref[-1]
	nimage=len(imnums)
	snParams=[['dt_%s'%i,'mu_%s'%i] for i in imnums if i!=refnum]
	all_vparam_names=np.append(args['params'],
							   snParams).flatten()
	
	ims=list(args['curves'].images.keys())

	
	for param in all_vparam_names:
		if param not in args['bounds'].keys():
			if param.startswith('dt_'):
				if args['fit_prior'] is not None:
					im=[x for x in ims if x[-1]==param[-1]][0]
					args['bounds'][param]=np.array([-1,1])*3*np.sqrt(args['fit_prior'].parallel.time_delay_errors[im]**2+ \
													args['fit_prior'].parallel.time_delay_errors[ref]**2)+\
							(args['fit_prior'].parallel.time_delays[im]-args['fit_prior'].parallel.time_delays[ref])
					
				else:
					args['bounds'][param]=np.array(args['bounds']['td'])#+time_delays[im]
			elif param.startswith('mu_'):
				if args['fit_prior'] is not None:

					im=[x for x in ims if x[-1]==param[-1]][0]
					
					args['bounds'][param]=np.array([-1,1])*3*(args['fit_prior'].parallel.magnifications[im]/args['fit_prior'].parallel.magnifications[ref])*\
							np.sqrt((args['fit_prior'].parallel.magnification_errors[im]/args['fit_prior'].parallel.magnifications[im])**2+ \
									(args['fit_prior'].parallel.magnification_errors[ref]/args['fit_prior'].parallel.magnifications[ref])**2)\
							+(args['fit_prior'].parallel.magnifications[im]/args['fit_prior'].parallel.magnifications[ref])
					

				else:
					args['bounds'][param]=np.array(args['bounds']['mu'])#*magnifications[im]
			elif args['fit_prior'] is not None:
				par_ref=args['fit_prior'].parallel.fitOrder[0]
				if param not in args['fit_prior'].images[par_ref].param_quantiles.keys():
					continue
				args['bounds'][param]=3*np.array([args['fit_prior'].images[par_ref].param_quantiles[param][0]- \
												args['fit_prior'].images[par_ref].param_quantiles[param][1],
												args['fit_prior'].images[par_ref].param_quantiles[param][2]- \
												args['fit_prior'].images[par_ref].param_quantiles[param][1]])+ \
									  args['fit_prior'].images[par_ref].param_quantiles[param][1]

		elif args['fit_prior'] is not None:
			par_ref=args['fit_prior'].parallel.fitOrder[0]
			if param not in args['fit_prior'].images[par_ref].param_quantiles.keys():
				continue
			args['bounds'][param]=3*np.array([args['fit_prior'].images[par_ref].param_quantiles[param][0]- \
											  args['fit_prior'].images[par_ref].param_quantiles[param][1],
											  args['fit_prior'].images[par_ref].param_quantiles[param][2]- \
											  args['fit_prior'].images[par_ref].param_quantiles[param][1]])+ \
								  args['fit_prior'].images[par_ref].param_quantiles[param][1]
	
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

	if args['fit_prior'] is not None and args['fit_prior'].images[args['fit_prior'].parallel.fitOrder[0]].fits.model._source.name not in args['models']:
		print('Wanted to use a fit prior but do not have the same model as an option.')
		raise RuntimeError
	elif args['fit_prior'] is not None:
		args['models']=args['fit_prior'].images[args['fit_prior'].parallel.fitOrder[0]].fits.model._source.name


	if args['max_n_bands'] is not None:
		best_bands=band_SNR[ref][:min(len(band_SNR[ref]),args['max_n_bands'])]
		temp_bands=[]
		for b in best_bands:
			temp_bands=np.append(temp_bands,np.where(args['curves'].images[ref].table['band']==b)[0])
		inds=temp_bands.astype(int)
	else:
		best_bands=args['bands']
		inds=np.arange(0,len(args['curves'].images[ref].table),1).astype(int)

	if 'ignore_models' in args['set_from_simMeta'].keys():
		to_ignore=args['curves'].images[ref].simMeta[args['set_from_simMeta']['ignore_models']]
		if isinstance(to_ignore,str):
			to_ignore=[to_ignore]
		args['models']=[x for x in np.array(args['models']).flatten() if x not in to_ignore]
	all_fit_dict={}
	if not args['curves'].quality_check(min_n_bands=args['min_n_bands'],
							min_n_points_per_band=args['min_points_per_band'],clip=False,method='parallel'):
		return
	if args['fast_model_selection'] and len(np.array(args['models']).flatten())>1:
		for b in args['force_positive_param']:
			if b in args['bounds'].keys():
				args['bounds'][b]=np.array([max([args['bounds'][b][0],0]),max([args['bounds'][b][1],0])])
			else:
				args['bounds'][b]=np.array([0,np.inf])
		
		minchisq=np.inf
		init_inds=deepcopy(inds)
		for mod in np.array(args['models']).flatten():
			inds=deepcopy(init_inds)
			if isinstance(mod,str):
				if mod.upper() in ['BAZIN','BAZINSOURCE']:
					mod='BAZINSOURCE'
					if len(np.unique(args['curves'].images[ref].table['band']))>1 and args['color_curve'] is None:
						best_band=band_SNR[args['fitOrder'][0]][0]
						inds=np.where(args['curves'].images[ref].table['band']==best_band)[0]
						

					source=BazinSource(data=args['curves'].images[ref].table[inds],colorCurve=args['color_curve'])
				else:
					source=sncosmo.get_source(mod)
				tempMod = sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
			else:
				tempMod=copy(mod)
			
			tempMod.set(**{k:args['constants'][k] for k in args['constants'].keys() if k in tempMod.param_names})
			tempMod.set(**{k:args['curves'].images[args['refImage']].simMeta[args['set_from_simMeta'][k]] for k in args['set_from_simMeta'].keys() if k in tempMod.param_names})
			if not np.all([tempMod.bandoverlap(x) for x in best_bands]):
				if args['verbose']:
					print('Skipping %s because it does not cover the bands...')
				continue
			if mod=='BAZINSOURCE':
				tempMod.set(z=0)
			try:
				res,fit=sncosmo.fit_lc(args['curves'].images[ref].table[inds],tempMod,[x for x in args['params'] if x in tempMod.param_names],
										bounds={b:args['bounds'][b] for b in args['bounds'] if b not in ['t0',tempMod.param_names[2]]},
										minsnr=args.get('minsnr',0))
			except:
				if args['verbose']:
					print('Issue with %s, skipping...'%mod)
				continue
			tempchisq=res.chisq/(len(inds)+len([x for x in args['params'] if x in tempMod.param_names])-1)
			if tempchisq<minchisq:
				minchisq=tempchisq
				bestres=copy(res)
				bestfit=copy(fit)
				bestmodname=copy(mod)
			all_fit_dict[mod]=[copy(fit),copy(res)]
		try:
			args['models']=[bestmodname]
		except:
			print('Every model had an error.')
			sys.exit(1)


	for mod in np.array(args['models']).flatten():



		if isinstance(mod,str):
			if mod.upper() in ['BAZIN','BAZINSOURCE']:
				source=BazinSource(data=args['curves'].images[args['fitOrder'][0]].table)
			else:
				source=sncosmo.get_source(mod)
			
			tempMod = sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
		else:
			tempMod=copy(mod)
		tempMod.set(**{k:args['constants'][k] for k in args['constants'].keys() if k in tempMod.param_names})
		tempMod.set(**{k:args['curves'].images[args['refImage']].simMeta[args['set_from_simMeta'][k]] for k in args['set_from_simMeta'].keys() if k in tempMod.param_names})
		

		if args['fit_prior'] is not None:

			par_ref=args['fit_prior'].parallel.fitOrder[0]
			if mod!=args['fit_prior'].images[par_ref].fits.model._source.name:
				continue

			temp_delays={k:args['fit_prior'].parallel.time_delays[k]-args['fit_prior'].parallel.time_delays[par_ref]\
						 for k in args['fit_prior'].parallel.fitOrder}
			temp_mags={k:args['fit_prior'].parallel.magnifications[k]/args['fit_prior'].parallel.magnifications[par_ref] \
						 for k in args['fit_prior'].parallel.fitOrder}
			args['curves'].combine_curves(time_delays={im:0 for im in args['curves'].images.keys()},
				magnifications={im:1 for im in args['curves'].images.keys()},minsnr=args.get('minsnr',0))
			args['curves'].series.meta['reft0']=args['fit_prior'].images[par_ref].fits.model.get('t0')
			args['curves'].series.meta['refamp']=args['fit_prior'].images[par_ref].fits.model.get(tempMod.param_names[2])
			args['curves'].series.meta['td']=temp_delays
			args['curves'].series.meta['mu']=temp_mags
		else:
			par_ref=args['refImage']
			im_name=args['refImage'][:-1]
			if args['trial_fit']:
				for b in args['force_positive_param']:
					if b in args['bounds'].keys():
						args['bounds'][b]=np.array([max([args['bounds'][b][0],0]),max([args['bounds'][b][1],0])])
					else:
						args['bounds'][b]=np.array([0,np.inf])

				nbands=args['max_n_bands'] if args['max_n_bands'] is not None else 2
				best_bands=band_SNR[args['refImage']][:min(len(band_SNR[args['refImage']]),nbands)]
				temp_delays={}
				temp_mags={}
				fit_order=np.flip(args['fitOrder']) if args['fitOrder'] is not None else \
										[x for x in args['curves'].images.keys() if x != args['refImage']]+[args['refImage']]
				for im in fit_order:
					temp_bands=[]
					for b in best_bands:
						temp_bands=np.append(temp_bands,np.where(args['curves'].images[im].table['band']==b)[0])
					temp_inds=temp_bands.astype(int)
					
					res,fit=sncosmo.fit_lc(deepcopy(args['curves'].images[im].table[temp_inds]),tempMod,[x for x in args['params'] if x in tempMod.param_names],
											bounds={b:args['bounds'][b] for b in args['bounds'].keys() if b not in ['t0',tempMod.param_names[2]]},
											minsnr=args.get('minsnr',0))
					temp_delays[im]=fit.get('t0')

					temp_mags[im]=fit.parameters[2]
				args['curves'].series.meta['reft0']=temp_delays[args['refImage']]
				args['curves'].series.meta['refamp']=temp_mags[args['refImage']]

				
				temp_delays={im:temp_delays[im]-temp_delays[args['refImage']] for im in temp_delays.keys()}
				temp_mags={im:temp_mags[im]/temp_mags[args['refImage']] for im in temp_mags}
				
				for b in args['bounds']:
					if b in list(res.errors.keys()):
						if b not in ['t0',tempMod.param_names[2]]:
							args['bounds'][b]=(np.array(args['bounds'][b])-np.median(args['bounds'][b]))/2+fit.get(b)
						elif b=='t0':
							args['bounds'][b]=(np.array(args['bounds'][b])-np.median(args['bounds'][b]))/2+args['curves'].series.meta['reft0']
						else:
							args['bounds'][b]=(np.array(args['bounds'][b])-np.median(args['bounds'][b]))/2+args['curves'].series.meta['refamp']

					elif b.startswith('dt_'):
						args['bounds'][b]=np.array(args['bounds']['td'])/2+temp_delays[im_name+b[-1]]
					elif b.startswith('mu_'):
						args['bounds'][b]=(np.array(args['bounds']['mu'])*temp_mags[im_name+b[-1]]+temp_mags[im_name+b[-1]])/2
				if tempMod.param_names[2] not in args['bounds'].keys():
					args['bounds'][tempMod.param_names[2]]=(np.array(args['bounds']['mu'])*args['curves'].series.meta['refamp']+args['curves'].series.meta['refamp'])/2
				if 't0' not in args['bounds'].keys():
					args['bounds']['t0']=np.array(args['bounds']['td'])/2+args['curves'].series.meta['reft0']

				args['curves'].combine_curves(time_delays={im:0 for im in args['curves'].images.keys()},
						magnifications={im:1 for im in args['curves'].images.keys()},minsnr=args.get('minsnr',0))
				args['curves'].series.meta['td']=temp_delays
				args['curves'].series.meta['mu']=temp_mags
			else:
				args['curves'].combine_curves(referenceImage=args['refImage'],static=True,model=tempMod,minsnr=args.get('minsnr',0))
				for b in args['bounds']:
					if b.startswith('dt_'):
						args['bounds'][b]=np.array(args['bounds']['td'])+args['curves'].series.meta['td'][im_name+b[-1]]
					elif b.startswith('mu_'):
						args['bounds'][b]=np.array(args['bounds']['mu'])*args['curves'].series.meta['mu'][im_name+b[-1]]
					elif b=='t0':
						args['bounds'][b]=np.array(args['bounds'][b])+args['curves'].series.meta['reft0']

				if tempMod.param_names[2] not in args['bounds'].keys():
					args['bounds'][tempMod.param_names[2]]=np.array(args['bounds']['mu'])*args['curves'].series.meta['refamp']
				if 't0' not in args['bounds'].keys():
					args['bounds']['t0']=np.array(args['bounds']['td'])+args['curves'].series.meta['reft0']


		

		

		if args['cut_time'] is not None:
			
			for im in args['curves'].images.keys():
				args['curves'].series.table=args['curves'].series.table[np.where(np.logical_or(args['curves'].series.table['image']!=im,
																	np.logical_and(args['curves'].series.table['time']>= \
																	args['cut_time'][0]*(1+tempMod.get('z'))+args['curves'].series.meta['reft0']+\
																	args['curves'].series.meta['td'][im],
																	args['curves'].series.table['time']<= \
																	args['cut_time'][1]*(1+tempMod.get('z'))+args['curves'].series.meta['reft0']+\
																	args['curves'].series.meta['td'][im])))[0]]
		for b in args['force_positive_param']:
			if b in args['bounds'].keys():
				args['bounds'][b]=np.array([max([args['bounds'][b][0],0]),max([args['bounds'][b][1],0])])
			else:
				args['bounds'][b]=np.array([0,np.inf])
		for b in [x for x in np.unique(args['curves'].series.table['band']) if x not in args['curves'].series.bands]:
			args['curves'].series.table=args['curves'].series.table[args['curves'].series.table['band']!=b]
		if not args['curves'].quality_check(min_n_bands=args['min_n_bands'],
							min_n_points_per_band=args['min_points_per_band'],clip=args['clip_data'],method='series'):
			return
		params,res,model=nest_series_lc(args['curves'].series.table,tempMod,nimage,bounds=args['bounds'],
									  vparam_names=[x for x in all_vparam_names if x in tempMod.param_names or x in np.array(snParams).flatten()],ref=par_ref,
									  minsnr=args.get('minsnr',5.),priors=args.get('priors',None),ppfs=args.get('ppfs',None),
									  method=args.get('nest_method','single'),maxcall=args.get('maxcall',None),
									  modelcov=args.get('modelcov',None),rstate=args.get('rstate',None),
									  maxiter=args.get('maxiter',None),npoints=args.get('npoints',100))
		if finallogz<res.logz:
			finallogz=res.logz
			final_param_quantiles,finalres,finalmodel=params,res,model
			time_delays=args['curves'].series.meta['td']
			magnifications=args['curves'].series.meta['mu']

	args['curves'].series.param_quantiles={d:final_param_quantiles[finalres.vparam_names.index(d)] \
											  for d in finalres.vparam_names}
	args['curves'].series.time_delays=dict([])
	args['curves'].series.magnifications=dict([])
	args['curves'].series.magnification_errors=dict([])
	args['curves'].series.time_delay_errors=dict([])

	args['curves'].series.t_peaks=dict([])
	args['curves'].series.a_peaks=dict([])

	
	if par_ref==args['refImage']:
		args['curves'].series.time_delays[par_ref]=0
		args['curves'].series.time_delay_errors[par_ref]=np.array([0,0])
		args['curves'].series.magnifications[par_ref]=1
		args['curves'].series.magnification_errors[par_ref]=np.array([0,0])
		args['curves'].series.t_peaks[par_ref]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('t0')],.5,finalres.weights)
		args['curves'].series.a_peaks[par_ref]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index(finalmodel.param_names[2])],
																			.5,finalres.weights)
		for k in args['curves'].images.keys():
			if k==par_ref:
				continue
			else:
				args['curves'].series.t_peaks[k]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('dt_'+k[-1])]+\
												finalres.samples[:,finalres.vparam_names.index('t0')],
												.5,finalres.weights)
				args['curves'].series.a_peaks[k]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('mu_'+k[-1])]*\
												finalres.samples[:,finalres.vparam_names.index(finalmodel.param_names[2])],
												.5,finalres.weights)
				dt_quant=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('dt_'+k[-1])],[.16,.5,.84],finalres.weights)
				mu_quant=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('mu_'+k[-1])],[.16,.5,.84],finalres.weights)
				args['curves'].series.time_delays[k]=dt_quant[1]
				args['curves'].series.time_delay_errors[k]=np.array([dt_quant[0]-dt_quant[1],dt_quant[2]-dt_quant[1]])
				args['curves'].series.magnifications[k]=mu_quant[1]
				args['curves'].series.magnification_errors[k]=np.array([mu_quant[0]-mu_quant[1],mu_quant[2]-mu_quant[1]])
	else:
		args['curves'].series.time_delays[args['refImage']]=0
		args['curves'].series.time_delay_errors[args['refImage']]=np.array([0,0])
		args['curves'].series.magnifications[args['refImage']]=1
		args['curves'].series.magnification_errors[args['refImage']]=np.array([0,0])
		trefSamples=finalres.samples[:,finalres.vparam_names.index('dt_'+args['refImage'][-1])]
		arefSamples=finalres.samples[:,finalres.vparam_names.index('mu_'+args['refImage'][-1])]
		args['curves'].series.t_peaks[args['refImage']]=weighted_quantile(trefSamples+finalres.samples[:,finalres.vparam_names.index('t0')],.5,finalres.weights)
		args['curves'].series.a_peaks[args['refImage']]=weighted_quantile(arefSamples*finalres.samples[:,finalres.vparam_names.index(finalmodel.param_names[2])],
																			.5,finalres.weights)
		for k in args['curves'].images.keys():
			if k==args['refImage']:
				continue
			elif k==par_ref:
				args['curves'].series.t_peaks[k]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('t0')],.5,finalres.weights)
				args['curves'].series.a_peaks[k]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index(finalmodel.param_names[2])],
																			.5,finalres.weights)
				dt_quant=weighted_quantile(-1*trefSamples,[.16,.5,.84],finalres.weights)
				mu_quant=weighted_quantile(1./arefSamples,[.16,.5,.84],finalres.weights)
				args['curves'].series.time_delays[k]=dt_quant[1]
				args['curves'].series.time_delay_errors[k]=np.array([dt_quant[0]-dt_quant[1],dt_quant[2]-dt_quant[1]])
				args['curves'].series.magnifications[k]=mu_quant[1]
				args['curves'].series.magnification_errors[k]=np.array([mu_quant[0]-mu_quant[1],mu_quant[2]-mu_quant[1]])
			else:
				args['curves'].series.t_peaks[k]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('dt_'+k[-1])]+\
												finalres.samples[:,finalres.vparam_names.index('t0')],
												.5,finalres.weights)
				args['curves'].series.a_peaks[k]=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('mu_'+k[-1])]*\
												finalres.samples[:,finalres.vparam_names.index(finalmodel.param_names[2])],
												.5,finalres.weights)
				dt_quant=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('dt_'+k[-1])]-trefSamples,[.16,.5,.84],finalres.weights)
				mu_quant=weighted_quantile(finalres.samples[:,finalres.vparam_names.index('mu_'+k[-1])]/arefSamples,[.16,.5,.84],finalres.weights)
				args['curves'].series.time_delays[k]=dt_quant[1]
				args['curves'].series.time_delay_errors[k]=np.array([dt_quant[0]-dt_quant[1],dt_quant[2]-dt_quant[1]])
				args['curves'].series.magnifications[k]=mu_quant[1]
				args['curves'].series.magnification_errors[k]=np.array([mu_quant[0]-mu_quant[1],mu_quant[2]-mu_quant[1]])
	
	args['curves'].combine_curves(time_delays=args['curves'].series.time_delays,magnifications=args['curves'].series.magnifications,referenceImage=args['refImage'])
	args['curves'].series.meta['td']=time_delays
	args['curves'].series.meta['mu']=magnifications

	finalmodel.set(t0=args['curves'].series.t_peaks[args['refImage']])

	finalmodel.parameters[2]=args['curves'].series.a_peaks[args['refImage']]

	args['curves'].series.refImage=args['refImage']
	args['curves'].series.priorImage=par_ref
	args['curves'].series.fits=newDict()
	args['curves'].series.fits['model']=finalmodel
	args['curves'].series.fits['res']=finalres


	if args['microlensing'] is not None:
		tempTable=deepcopy(args['curves'].series.table)
		micro,sigma,x_pred,y_pred,samples=fit_micro(args['curves'].series.fits.model,tempTable,
													tempTable['zpsys'][0],args['nMicroSamples'],
													micro_type=args['microlensing'],kernel=args['kernel'])

		temp_vparam_names=args['curves'].series.fits.res.vparam_names+[finalmodel.param_names[2]]+['t0']
		for im in args['curves'].images.keys():
			try:
				temp_vparam_names.remove('dt_'+str(im[-1]))
				temp_vparam_names.remove('mu_'+str(im[-1]))
			except:
				pass
		temp_bounds={p:args['curves'].series.param_quantiles[p][[0,2]] \
					 for p in args['curves'].series.fits.res.vparam_names}

		temp_bounds['t0']=args['bounds']['td']+args['curves'].series.t_peaks[args['refImage']]
		temp_bounds={b:temp_bounds[b] for b in temp_bounds.keys() if b!=args['curves'].series.fits.model.param_names[2]}
		try:
			t0s=pyParz.foreach(samples.T,_micro_uncertainty,
						   [args['curves'].series.fits.model,np.array(tempTable),tempTable.colnames,
							x_pred,temp_vparam_names,
							temp_bounds,None,args.get('minsnr',0),args.get('maxcall',None)])
		except:
				if args['verbose']:
					print('Issue with microlensing identification, skipping...')
				return args['curves']
		else:
			return args['curves']
		mu,sigma=scipy.stats.norm.fit(t0s)

		args['curves'].series.param_quantiles['micro']=np.sqrt((args['curves'].series.fits.model.get('t0')-mu)**2 \
																  +9*sigma**2)

	return args['curves']




def nest_series_lc(data,model,nimage,vparam_names,bounds,ref='image_1',
					 minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
					 maxiter=None, maxcall=None, modelcov=False, rstate=None,
					 verbose=False, warn=True,**kwargs):
	
	####Taken from SNCosmo nest_lc
	# experimental parameters
	tied = kwargs.get("tied", None)

	vparam_names=list(vparam_names)
	if ppfs is None:
		ppfs = {}
	if tied is None:
		tied = {}

	

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

	model_param_names=[x for x in vparam_names[:len(vparam_names)-(nimage-1)*2]]
	model_idx = np.array([vparam_names.index(name) for name in model_param_names])
	td_params=[x for x in vparam_names[len(vparam_names)-nimage*2:] if x.startswith('dt')]
	td_idx=np.array([vparam_names.index(name) for name in td_params])
	amp_params=[x for x in vparam_names[len(vparam_names)-nimage*2:] if x.startswith('mu')]
	amp_idx=np.array([vparam_names.index(name) for name in amp_params])
	
	#mindat=model.mintime()
	#maxdat=model.maxtime()
	#data=data[np.where(np.logical_and(data['time']>=mindat,data['time']<=maxdat))]
	im_indices=[np.where(data['image']==i)[0] for i in np.unique(data['image']) if i !=ref]
	cov = np.diag(data['fluxerr']**2)
	def chisq_likelihood(parameters):
		model.set(**{model_param_names[k]:parameters[model_idx[k]] for k in range(len(model_idx))})
		all_data=deepcopy(data)
		for i in range(len(im_indices)):
			all_data['time'][im_indices[i]]-=parameters[td_idx[i]]
			all_data['flux'][im_indices[i]]/=parameters[amp_idx[i]]
		all_data.sort('time')
		model_observations = model.bandflux(all_data['band'],all_data['time'],
											zp=all_data['zp'],zpsys=all_data['zpsys'])
		if modelcov:
			
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



	return params,res,model


def getBandSNR(curves,bands,min_points_per_band):
	final_bands=[]
	band_dict={im:[] for im in curves.images.keys()}
	for band in list(bands):
		to_add=True
		for im in curves.images.keys():
			if len(np.where(curves.images[im].table['band']==band)[0])<min_points_per_band:
				to_add=False
			else:
				band_dict[im].append(band)
		if to_add:
			final_bands.append(band)
	
	all_SNR=[]
	band_SNR={im:[] for im in curves.images.keys()}
	for d in curves.images.keys():
		for band in final_bands:	
			inds=np.where(curves.images[d].table['band']==band)[0]
			if len(inds)==0:
				band_SNR[d].append(0)
			else:
				band_SNR[d].append(np.sum(curves.images[d].table['flux'][inds]/curves.images[d].table['fluxerr'][inds])*\
					 np.sqrt(len(inds)))
		
	band_SNR={k:np.array(final_bands)[np.flip(np.argsort(band_SNR[k]))] for k in band_SNR.keys()}
	return(np.array(final_bands),band_SNR,band_dict)

def _fitparallel(all_args):
	if isinstance(all_args,(list,tuple,np.ndarray)):
		curves,args=all_args
		if isinstance(args,list):
			args=args[0]
		if isinstance(curves,list):
			curves,single_par_vars=curves
			for key in single_par_vars:
				args[key]=single_par_vars[key]

		if isinstance(curves,str):
			args['curves']=pickle.load(open(curves,'rb'))
		else:
			args['curves']=curves
			if args['verbose']:
				print('Fitting MISN number %i...'%curves.nsn)
	else:
		args=all_args
	

	for p in args['curves'].constants.keys():
		if p not in args['constants'].keys():
			args['constants'][p]=args['curves'].constants[p]

	if 't0' in args['bounds']:
		t0Bounds=copy(args['bounds']['t0'])

	if args['clip_data']:
		for im in args['curves'].images.keys():
			args['curves'].clip_data(im=im,minsnr=args.get('minsnr',0),max_cadence=args['max_cadence'])

	args['bands'],band_SNR,band_dict=getBandSNR(args['curves'],args['bands'],args['min_points_per_band'])
	args['curves'].bands=args['bands']
	if len(args['bands'])==0:
		return(None)
	
	for d in args['curves'].images.keys():
		for b in [x for x in np.unique(args['curves'].images[d].table['band']) if x not in band_dict[d]]:
			args['curves'].images[d].table=args['curves'].images[d].table[args['curves'].images[d].table['band']!=b]


	if 'amplitude' in args['bounds']:
		args['guess_amplitude']=False


	if args['fitOrder'] is None:
		all_SNR=[np.sum(args['curves'].images[d].table['flux']/args['curves'].images[d].table['fluxerr']) \
					for d in np.sort(list(args['curves'].images.keys()))]
		sorted=np.flip(np.argsort(all_SNR))
		args['fitOrder']=np.sort(list(args['curves'].images.keys()))[sorted]

	args['curves'].parallel.fitOrder=args['fitOrder']

	if args['t0_guess'] is not None:
		if 't0' in args['bounds']:
			args['bounds']['t0']=(t0Bounds[0]+args['t0_guess'][args['fitOrder'][0]],t0Bounds[1]+args['t0_guess'][args['fitOrder'][0]])
			guess_t0=args['t0_guess']
		else:
			print('If you supply a t0 guess, you must also supply bounds.')
			sys.exit(1)

	if args['max_n_bands'] is not None:
		best_bands=band_SNR[args['fitOrder'][0]][:min(len(band_SNR[args['fitOrder'][0]]),args['max_n_bands'])]
		temp_bands=[]
		for b in best_bands:
			temp_bands=np.append(temp_bands,np.where(args['curves'].images[args['fitOrder'][0]].table['band']==b)[0])
		inds=temp_bands.astype(int)
	else:
		best_bands=args['bands']
		inds=np.arange(0,len(args['curves'].images[args['fitOrder'][0]].table),1).astype(int)
	initial_bounds=deepcopy(args['bounds'])
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

	if 'ignore_models' in args['set_from_simMeta'].keys():
		to_ignore=args['curves'].images[args['fitOrder'][0]].simMeta[args['set_from_simMeta']['ignore_models']]
		if isinstance(to_ignore,str):
			to_ignore=[to_ignore]
		args['models']=[x for x in np.array(args['models']).flatten() if x not in to_ignore]
	if not args['curves'].quality_check(min_n_bands=args['min_n_bands'],
								min_n_points_per_band=args['min_points_per_band'],clip=args['clip_data']):
		return
	all_fit_dict={}
	if args['fast_model_selection'] and len(np.array(args['models']).flatten())>1:
		for b in args['force_positive_param']:
			if b in args['bounds'].keys():
				args['bounds'][b]=np.array([max([args['bounds'][b][0],0]),max([args['bounds'][b][1],0])])
			else:
				args['bounds'][b]=np.array([0,np.inf])
		
		minchisq=np.inf
		init_inds=deepcopy(inds)
		for mod in np.array(args['models']).flatten():
			inds=deepcopy(init_inds)
			if isinstance(mod,str):
				if mod.upper() in ['BAZIN','BAZINSOURCE']:
					mod='BAZINSOURCE'
					if len(np.unique(args['curves'].images[args['fitOrder'][0]].table['band']))>1 and args['color_curve'] is None:
						best_band=band_SNR[args['fitOrder'][0]][0]
						inds=np.where(args['curves'].images[args['fitOrder'][0]].table['band']==best_band)[0]
						

					source=BazinSource(data=args['curves'].images[args['fitOrder'][0]].table[inds],colorCurve=args['color_curve'])
				else:
					source=sncosmo.get_source(mod)
				tempMod = sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
			else:
				tempMod=copy(mod)
			
			tempMod.set(**{k:args['constants'][k] for k in args['constants'].keys() if k in tempMod.param_names})
			tempMod.set(**{k:args['curves'].images[args['refImage']].simMeta[args['set_from_simMeta'][k]] for k in args['set_from_simMeta'].keys() if k in tempMod.param_names})
			if not np.all([tempMod.bandoverlap(x) for x in best_bands]):
				if args['verbose']:
					print('Skipping %s because it does not cover the bands...'%mod)
				continue
			if mod=='BAZINSOURCE':
				tempMod.set(z=0)
			try:
				res,fit=sncosmo.fit_lc(args['curves'].images[args['fitOrder'][0]].table[inds],tempMod,[x for x in args['params'] if x in tempMod.param_names],
										bounds={b:args['bounds'][b] for b in args['bounds'] if b not in ['t0',tempMod.param_names[2]]},
										minsnr=args.get('minsnr',0))
			except:
				if args['verbose']:
					print('Issue with %s, skipping...'%mod)
				continue
			tempchisq=res.chisq/(len(inds)+len([x for x in args['params'] if x in tempMod.param_names])-1)
			if tempchisq<minchisq:
				minchisq=tempchisq
				bestres=copy(res)
				bestfit=copy(fit)
				bestmodname=copy(mod)
			all_fit_dict[mod]=[copy(fit),copy(res)]
		try:
			args['models']=[bestmodname]
		except:
			print('Every model had an error.')
			return None
	for mod in np.array(args['models']).flatten():
		if isinstance(mod,str):
			if mod.upper() in ['BAZIN','BAZINSOURCE']:
				mod='BAZINSOURCE'
				if len(np.unique(args['curves'].images[args['fitOrder'][0]].table['band']))>1 and args['color_curve'] is None:
					best_band=band_SNR[args['fitOrder'][0]][0]
					inds=np.where(args['curves'].images[args['fitOrder'][0]].table['band']==best_band)[0]
						

				source=BazinSource(data=args['curves'].images[args['fitOrder'][0]].table[inds],colorCurve=args['color_curve'])
			else:
				source=sncosmo.get_source(mod)
			tempMod = sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
		else:
			tempMod=copy(mod)
		
		tempMod.set(**{k:args['constants'][k] for k in args['constants'].keys() if k in tempMod.param_names})
		if args['set_from_simMeta'] is not None:
			tempMod.set(**{k:args['curves'].images[args['refImage']].simMeta[args['set_from_simMeta'][k]] for k in args['set_from_simMeta'].keys() if k in tempMod.param_names})
		if mod=='BAZINSOURCE':
			tempMod.set(z=0)
		guess_t0,guess_amp=sncosmo.fitting.guess_t0_and_amplitude( \
			sncosmo.photdata.photometric_data(args['curves'].images[args['fitOrder'][0]].table[inds]),
			tempMod,args.get('minsnr',5.))
		if 't0' in args['bounds'] and args['t0_guess'] is None:

			args['bounds']['t0']=np.array(initial_bounds['t0'])+guess_t0
		
		if args['trial_fit'] and args['t0_guess'] is None:
			for b in args['force_positive_param']:
				if b in args['bounds'].keys():
					args['bounds'][b]=np.array([max([args['bounds'][b][0],0]),max([args['bounds'][b][1],0])])
				else:
					args['bounds'][b]=np.array([0,np.inf])
			if args['max_n_bands'] is None:
				best_bands=band_SNR[args['fitOrder'][0]][:min(len(band_SNR[args['fitOrder'][0]]),2)]
				temp_bands=[]
				for b in best_bands:
					temp_bands=np.append(temp_bands,np.where(args['curves'].images[args['fitOrder'][0]].table['band']==b)[0])
				temp_inds=temp_bands.astype(int)
			res,fit=sncosmo.fit_lc(args['curves'].images[args['fitOrder'][0]].table[temp_inds],tempMod,[x for x in args['params'] if x in tempMod.param_names],
									bounds={b:args['bounds'][b]+(args['bounds'][b]-np.median(args['bounds'][b]))*2 if b=='t0' else args['bounds'][b] for b in args['bounds'] if b!= tempMod.param_names[2]},
									minsnr=args.get('minsnr',0))

			args['bounds']={b:(np.array(args['bounds'][b])-np.median(args['bounds'][b]))/2+fit.get(b) if b in res.param_names else args['bounds'][b] for b in args['bounds'].keys()}
			if tempMod.param_names[2] not in args['bounds'].keys():
				args['bounds'][tempMod.param_names[2]]=np.array([.1,10])*fit.parameters[2]
			
			guess_t0=fit.get('t0')
		elif args['guess_amplitude']:
			if tempMod.param_names[2] in args['bounds']:
				args['bounds'][tempMod.param_names[2]]=np.array(args['bounds'][tempMod.param_names[2]])*\
					guess_amp
			else:
				args['bounds'][tempMod.param_names[2]]=[.1*guess_amp,10*guess_amp]

		if args['clip_data']:
			args['curves'].images[args['fitOrder'][0]].table=args['curves'].images[args['fitOrder'][0]].table[inds]
			if args['cut_time'] is not None:
				args['curves'].clip_data(im=args['fitOrder'][0],minsnr=args.get('minsnr',0),mintime=args['cut_time'][0]*(1+tempMod.get('z')),
						maxtime=args['cut_time'][1]*(1+tempMod.get('z')),peak=guess_t0)
			else:
				args['curves'].clip_data(im=args['fitOrder'][0],minsnr=args.get('minsnr',0))
			fit_table=args['curves'].images[args['fitOrder'][0]].table
		elif args['cut_time'] is not None:
			fit_table=deepcopy(args['curves'].images[args['fitOrder'][0]].table)
			fit_table=fit_table[inds]
			fit_table=fit_table[fit_table['time']>=guess_t0+(args['cut_time'][0]*(1+tempMod.get('z')))]
			fit_table=fit_table[fit_table['time']<=guess_t0+(args['cut_time'][1]*(1+tempMod.get('z')))]
			fit_table=fit_table[fit_table['flux']/fit_table['fluxerr']>=args.get('minsnr',0)]
			
		else:
			fit_table=deepcopy(args['curves'].images[args['fitOrder'][0]].table)
			fit_table=fit_table[inds]
		for b in args['force_positive_param']:
			if b in args['bounds'].keys():
				args['bounds'][b]=np.array([max([args['bounds'][b][0],0]),max([args['bounds'][b][1],0])])
			else:
				args['bounds'][b]=np.array([0,np.inf])
		res,fit=sncosmo.nest_lc(fit_table,tempMod,[x for x in args['params'] if x in tempMod.param_names],
								bounds=args['bounds'],
							  priors=args.get('priors',None), ppfs=args.get('ppfs',None),
								minsnr=args.get('minsnr',5.0), method=args.get('nest_method','single'),
							  maxcall=args.get('maxcall',None), modelcov=args.get('modelcov',False),
							  rstate=args.get('rstate',None),guess_amplitude_bound=False,
							  zpsys=args['curves'].images[args['fitOrder'][0]].zpsys,
							  maxiter=args.get('maxiter',None),npoints=args.get('npoints',100))
		all_fit_dict[mod]=[copy(fit),copy(res)]
		
		if finallogz<res.logz:
			first_res=[args['fitOrder'][0],copy(fit),copy(res)]
			finallogz=res.logz

	first_params=[weighted_quantile(first_res[2].samples[:,i],[.16,.5,.84],first_res[2].weights)\
				  for i in range(len(first_res[2].vparam_names))]

	args['curves'].images[args['fitOrder'][0]].fits=newDict()
	args['curves'].images[args['fitOrder'][0]].fits['model']=first_res[1]
	args['curves'].images[args['fitOrder'][0]].fits['res']=first_res[2]




	t0ind=first_res[2].vparam_names.index('t0')
	ampind=first_res[2].vparam_names.index(first_res[1].param_names[2])

	args['curves'].images[args['fitOrder'][0]].param_quantiles={k:first_params[first_res[2].vparam_names.index(k)] for\
																 k in first_res[2].vparam_names}
	for i in range(len(first_res[2].vparam_names)):
		if first_res[2].vparam_names[i]==first_res[1].param_names[2] or first_res[2].vparam_names[i]=='t0':
			continue
		initial_bounds[first_res[2].vparam_names[i]]=3*np.array([first_params[i][0],first_params[i][2]])-2*first_params[i][1]
	for d in args['fitOrder'][1:]:
		if args['max_n_bands'] is not None:
			best_bands=band_SNR[d][:min(len(band_SNR[d]),args['max_n_bands'])]
			temp_bands=[]
			for b in best_bands:
				temp_bands=np.append(temp_bands,np.where(args['curves'].images[d].table['band']==b)[0])
			inds=temp_bands.astype(int)
		else:
			inds=np.arange(0,len(args['curves'].images[d].table),1).astype(int)
		args['curves'].images[d].fits=newDict()
		initial_bounds['t0']=deepcopy(t0Bounds)
		
		if args['t0_guess'] is not None:
			if 't0' in args['bounds']:
				initial_bounds['t0']=(t0Bounds[0]+args['t0_guess'][d],t0Bounds[1]+args['t0_guess'][d])
		else:
			best_bands=band_SNR[d][:min(len(band_SNR[d]),2)]
			temp_bands=[]
			for b in best_bands:
				temp_bands=np.append(temp_bands,np.where(args['curves'].images[d].table['band']==b)[0])
			inds=temp_bands.astype(int)

		if mod=='BAZINSOURCE':
			minds=np.where(args['curves'].images[d].table['band']==best_band)[0]
			inds=None
		else:
			minds=np.arange(0,len(args['curves'].images[d].table),1).astype(int)
		if args['clip_data']:
			fit_table=args['curves'].images[d].table[minds]
		else:
			fit_table=deepcopy(args['curves'].images[d].table)
			fit_table=fit_table[minds]
		
		par_output=nest_parallel_lc(fit_table,first_res[1],first_res[2],initial_bounds,min_n_bands=args['min_n_bands'],
						min_n_points_per_band=args['min_points_per_band'],
						guess_amplitude_bound=True,priors=args.get('priors',None), ppfs=args.get('None'),
						method=args.get('nest_method','single'),cut_time=args['cut_time'],snr_band_inds=inds,
						maxcall=args.get('maxcall',None), modelcov=args.get('modelcov',False),
						rstate=args.get('rstate',None),minsnr=args.get('minsnr',5),
						maxiter=args.get('maxiter',None),npoints=args.get('npoints',1000))
		if par_output is None:
			return
		params,args['curves'].images[d].fits['model'],args['curves'].images[d].fits['res']=par_output

	sample_dict={args['fitOrder'][0]:[first_res[2].samples[:,t0ind],first_res[2].samples[:,ampind]]}
	for k in args['fitOrder'][1:]:
		sample_dict[k]=[args['curves'].images[k].fits['res'].samples[:,t0ind],
						args['curves'].images[k].fits['res'].samples[:,ampind]]
		args['curves'].images[k].param_quantiles={d:params[args['curves'].images[k].fits['res'].vparam_names.index(d)] \
												  for d in args['curves'].images[k].fits['res'].vparam_names}
	trefSamples,arefSamples=sample_dict[args['refImage']]
	refWeights=args['curves'].images[args['refImage']].fits['res'].weights

	args['curves'].parallel.time_delays={args['refImage']:0}
	args['curves'].parallel.magnifications={args['refImage']:1}
	args['curves'].parallel.time_delay_errors={args['refImage']:np.array([0,0])}
	args['curves'].parallel.magnification_errors={args['refImage']:np.array([0,0])}
	for k in args['curves'].images.keys():
		if k==args['refImage']:
			continue
		else:
			ttempSamples,atempSamples=sample_dict[k]

			if len(ttempSamples)>len(trefSamples):
				inds=np.flip(np.argsort(args['curves'].images[k].fits['res'].weights)[len(ttempSamples)-len(trefSamples):])
				inds_ref=np.flip(np.argsort(refWeights))
			else:
				inds_ref=np.flip(np.argsort(refWeights)[len(trefSamples)-len(ttempSamples):])
				inds=np.flip(np.argsort(args['curves'].images[k].fits['res'].weights))

			t_quant=weighted_quantile(ttempSamples[inds]-trefSamples[inds_ref],[.16,.5,.84],refWeights[inds_ref]* \
									  args['curves'].images[k].fits['res'].weights[inds])
			a_quant=weighted_quantile(atempSamples[inds]/arefSamples[inds_ref],[.16,.5,.84],refWeights[inds_ref]* \
									  args['curves'].images[k].fits['res'].weights[inds])
			args['curves'].parallel.time_delays[k]=t_quant[1]
			args['curves'].parallel.magnifications[k]=a_quant[1]
			args['curves'].parallel.time_delay_errors[k]=np.array([t_quant[0]-t_quant[1],t_quant[2]-t_quant[1]])
			args['curves'].parallel.magnification_errors[k]= \
				np.array([a_quant[0]-a_quant[1],a_quant[2]-a_quant[1]])
		if args['clip_data']:
			if args['cut_time'] is not None:
				args['curves'].clip_data(im=k,minsnr=args.get('minsnr',0),mintime=args['cut_time'][0]*(1+args['curves'].images[k].fits.model.get('z')),
					maxtime=args['cut_time'][1]*(1+args['curves'].images[k].fits.model.get('z')),
						peak=args['curves'].images[k].fits.model.get('t0'))
			else:
				args['curves'].clip_data(im=k,minsnr=args.get('minsnr',0))

	if args['microlensing'] is not None:
		for k in args['curves'].images.keys():
			tempTable=deepcopy(args['curves'].images[k].table)
			micro,sigma,x_pred,y_pred,samples=fit_micro(args['curves'].images[k].fits.model,tempTable,
														args['curves'].images[k].zpsys,args['nMicroSamples'],
														micro_type=args['microlensing'],kernel=args['kernel'])


			try:

					
				t0s=pyParz.foreach(samples.T,_micro_uncertainty,
							   [args['curves'].images[k].fits.model,np.array(tempTable),tempTable.colnames,
								x_pred,args['curves'].images[k].fits.res.vparam_names,
								{p:args['curves'].images[k].param_quantiles[p][[0,2]]\
								 for p in args['curves'].images[k].fits.res.vparam_names if p != \
								 args['curves'].images[k].fits.model.param_names[2]},None,
								args.get('minsnr',0),args.get('maxcall',None)],numThreads=args['npar_cores'])
			except:
				if args['verbose']:
					print('Issue with microlensing identification, skipping...')
				return args['curves']
			mu,sigma=scipy.stats.norm.fit(t0s)
			args['curves'].images[k].param_quantiles['micro']=np.sqrt((args['curves'].images[k].fits.model.get('t0')-mu)**2\
																	  +9*sigma**2)


	return args['curves']

def nest_parallel_lc(data,model,prev_res,bounds,guess_amplitude_bound=False,cut_time=None,snr_band_inds=None,
				   min_n_bands=1,min_n_points_per_band=3,
				   minsnr=5., priors=None, ppfs=None, npoints=100, method='single',
				   maxiter=None, maxcall=None, modelcov=False, rstate=None,
				   verbose=False, warn=True,**kwargs):
	
	####Taken from SNCosmo nest_lc
	# experimental parameters
	tied = kwargs.get("tied", None)

	vparam_names=list(prev_res.vparam_names)
	if ppfs is None:
		ppfs = {}
	if tied is None:
		tied = {}

	model=copy(model)
	if guess_amplitude_bound:

		if snr_band_inds is None:
			snr_band_inds=np.arange(0,len(data),1).astype(int)

		guess_t0,guess_amp=sncosmo.fitting.guess_t0_and_amplitude(sncosmo.photdata.photometric_data(data[snr_band_inds]),
																  model,minsnr)

		model.set(t0=guess_t0)
		model.parameters[2]=guess_amp

		bounds[model.param_names[2]]=(0,10*guess_amp)
		bounds['t0']=np.array(bounds['t0'])+guess_t0

	if cut_time is not None and guess_amplitude_bound:
		data=data[data['time']>=cut_time[0]*(1+model.get('z'))+guess_t0]
		data=data[data['time']<=cut_time[1]*(1+model.get('z'))+guess_t0]

	data,quality=check_table_quality(data,min_n_bands=min_n_bands,min_n_points_per_band=min_n_points_per_band,clip=True)
	if not quality:
		return
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
		if p in __thetaL__:
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



def _micro_uncertainty(args):
	sample,other=args
	nest_fit,data,colnames,x_pred,vparam_names,bounds,priors,minsnr,maxcall=other
	data=Table(data,names=colnames)

	temp_nest_mod=deepcopy(nest_fit)
	tempMicro=AchromaticMicrolensing(x_pred/(1+nest_fit.get('z')),sample,magformat='multiply')
	temp_nest_mod.add_effect(tempMicro,'microlensing','rest')
	tempRes,tempMod=nest_lc(data,temp_nest_mod,vparam_names=vparam_names,bounds=bounds,minsnr=minsnr,maxcall=maxcall,
							guess_amplitude_bound=True,maxiter=None,npoints=200,priors=priors)

	return float(tempMod.get('t0'))

def fit_micro(fit,dat,zpsys,nsamples,micro_type='achromatic',kernel='RBF'):
	t0=fit.get('t0')
	fit.set(t0=t0)
	data=deepcopy(dat)
	data['time']-=t0

	data=data[data['time']<=40.]
	data=data[data['time']>=-15.]
	if len(data)==0:
		data=deepcopy(dat)

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
			plt.close()


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




		res,mod=nest_lc(args['curve'].table,mod,vparam_names=args['params'],bounds=args['bounds'],
						guess_amplitude_bound=guess_amp_bound,maxiter=1000,npoints=200)
	return({'res':res,'model':mod})


def identify_micro_func(args):
	print('Only a development function for now!')
	return args['bands'],args['bands']
	if len(args['bands'])<=2:
		return args['bands'],args['bands']
	res_dict={}
	original_args=copy(args)
	combos=[]
	for r in range(len(args['bands'])-1):
		temp=[x for x in itertools.combinations(original_args['bands'],r)]
		for t in temp:
			combos.append(t)
	if 'td' not in args['bounds'].keys():
		args['bounds']['td']=args['bounds']['t0']
	for bands in itertools.combinations(args['bands'],2):
		good=True
		for b in bands:
			if not np.all([len(np.where(original_args['curves'].images[im].table['band']==b)[0])>=3 for im in original_args['curves'].images.keys()]):
				good=False
		if not good:
			continue
		temp_args=copy(original_args)

		temp_args['bands']=[x for x in bands]
		temp_args['npoints']=200
		temp_args['fit_prior']=None

		fitCurves=_fitColor(temp_args)
		if np.all([np.isfinite(fitCurves.color.time_delays[x]) for x in fitCurves.images.keys()]):
			res_dict[bands[0]+'-'+bands[1]]=copy(fitCurves.color.fits.res)

	if len(list(res_dict.keys()))==0:
		print('No good fitting.',args['bands'])
		return(args['bands'],args['bands'])
	ind=res_dict[list(res_dict.keys())[0]].vparam_names.index('c')
	print([(x,weighted_quantile(res_dict[x].samples[:,ind],[.16,.5,.84],res_dict[x].weights)) for x in res_dict.keys()])
	dev_dict={}
	
	for bs in combos:
		dev_dict[','.join(list(bs))]=(np.average([weighted_quantile(res_dict[x].samples[:,ind],.5,res_dict[x].weights) \
												  for x in res_dict.keys() if np.all([b not in x for b in bs])],weights= \
													 1/np.abs([res_dict[x].logz for x in res_dict.keys() if np.all([b not in x for b in bs])])),
					np.std([weighted_quantile(res_dict[x].samples[:,ind],.5,res_dict[x].weights) \
			  for x in res_dict.keys() if np.all([b not in x for b in bs])]))
	print(dev_dict)

	to_remove=None
	
	best_std=dev_dict[''][1]/np.sqrt(len(args['bands']))
	
	if len(args['bands'])>3:
		for bands in dev_dict.keys():
			if dev_dict[bands][1]!=0:#len(args['bands'])-len(bands.split(','))==2:
				print(bands,dev_dict[bands])
				if dev_dict[bands][1]/np.sqrt(len(args['bands'])-len(bands.split(',')))<best_std:
					to_remove=bands.split(',')
					best_std=dev_dict[bands][1]/np.sqrt(len(args['bands'])-len(to_remove))
	print(to_remove,best_std)
	sys.exit()
	final_color_bands=None
	best_logz=-np.inf
	best_logzerr=0

	for bands in res_dict.keys():
		logz,logzerr=calc_ev(res_dict[bands],args['npoints'])
		if logz>best_logz:
			final_color_bands=bands
			best_logz=logz
			best_logzerr=logzerr
	print(bands, best_logz,best_logzerr)
	final_all_bands=[]
	for bands in res_dict.keys():
		logz,logzerr=calc_ev(res_dict[bands],args['npoints'])
		print(bands,logz,logzerr)
		if logz+3*logzerr>=best_logz-3*best_logzerr:
			final_all_bands=np.append(final_all_bands,bands.split('-'))


	print(np.unique(final_all_bands),np.array(final_color_bands.split('-')))
	sys.exit()
	return(np.unique(final_all_bands),np.array(final_color_bands.split('-')))

	#else:
	#	print([[x for x in args['bands'] if x not in to_remove]]*2)
	#	sys.exit()
	#	return [[x for x in args['bands'] if x not in to_remove]]*2

	# else:
	# 	best_bands=None
	# 	best_logz=-np.inf
	# 	for bands in res_dict.keys():
	#
	# 		if res_dict[bands].logz>best_logz:
	# 			best_bands=bands
	# 			best_logz=res_dict[bands].logz
	#
	# 	return [best_bands.split('-')]*2

def calc_ev(res,nlive):
	logZnestle = res.logz                         # value of logZ
	infogainnestle = res.h                        # value of the information gain in nats
	if not np.isfinite(infogainnestle):
		infogainnestle=.1*logZnestle
	logZerrnestle = np.sqrt(infogainnestle)#/nlive) # estimate of the statistcal uncertainty on logZ

	return logZnestle, logZerrnestle