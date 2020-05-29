from matplotlib import pyplot as plt
import numpy as np
import scipy
from scipy.integrate import quad
from matplotlib import ticker, rcParams

from copy import deepcopy,copy
import sncosmo
import nestle
import corner

from .util import weighted_quantile,plot

def EA(z1, z2, cosmo):
	"""The integral of the inverse of the normalized 
	Hubble parameter, from z1 to z2:  int_z1^z2{1/E(z')dz}
	Eq 8 of Coe & Moustakas 2009.
	Assuming a flat universe."""
	if 'Ode0' not in cosmo.keys():
		cosmo['Ode0'] = 1-cosmo['Om0'] # assuming flat universe
	#if w0>-0.3 or w0<-2:
	#    result = -np.inf
	#elif np.abs(wa)>5:
	#    result = -np.inf
	#else:
	if 'w' in cosmo.keys():
		Ez2 = lambda z: ( cosmo['Om0']*(1+z)**3 + 
					 cosmo['Ode0']*(1+z)**(3*(1+cosmo['w'])))
	else:
		Ez2 = lambda z: ( cosmo['Om0']*(1+z)**3 + 
					 cosmo['Ode0']*(1+z)**(3*(1+cosmo['w0']+cosmo['wa'])) * 
					 np.exp(-3*cosmo['wa']*z/(1+z)) )
	Ezinv = lambda z: 1 / np.sqrt(Ez2(z))
	# if Ez2(z1)>0 and Ez2(z2)>0:
	allz=np.append(z1,z2)
	sort=np.argsort(allz)
	n1=len(z1)
	result=np.zeros(n1*2)

	last=0
	for i in range(len(sort)):
		if i>0:
			temp,_=quad(Ezinv,allz[sort[i-1]],allz[sort[i]])
		else:
			temp,_=quad(Ezinv,0,allz[sort[i]])

		result[sort[i]]=last+temp
		last+=temp
		
	
	#result, precision = quad(Ezinv, z1, z2)
	#else:
	#    result = -np.inf
	return(result[:n1],result[n1:])


def Eratio(zl, zs, cosmo):
	"""The time delay expansion function ratio 
	(script E, eq 12 in Coe & Moustakas 2009).
	Assuming a flat universe."""

	EL,ES = EA(zl,zs, cosmo)
	#ES = EA(0,zs, Om0, w0, wa,iszero=True)
	ELS = ES-EL#EA(zl, zs, Om0, w0, wa)
	
	#if EL<=0 or ES<=0 or ELS<=0:
	#    Erat = -np.inf
	#else:
	Erat = EL * ES / ELS
	return(Erat/cosmo['h'])

def Evariance(Erat, dTc):
	"""The variance on the E ratio due to the 
	time delay and lens potential uncertainty. 
	dTc : percent uncertainty on time delay distance ratio 
		(combining time delay measurement + lens modeling error)  
		(This is a percentage:  enter "5" for 5% precision)
	"""
	return( Erat**2 * (dTc/100.)**2 )


def loglikelihoodE(testVars,testPoint, zl, zs, dTc=0.1, set_cosmo={},
				   Om0=0.3, w0true=-1., watrue=0.,htrue=.7): 
	"""log10(likelihood) for a time delay distance 
	ratio constraint, comparing 
	a given test point in w0,wa space to the true (w0,wa) position. 
	
	w0wa: [w0,wa], the test point
	zl : lens redshift
	zs : source redshift
	dTc : percent precision on cosmological distance ratio, combining
	   uncertainty from the time delay measurement + lens modeling
	Om0, w0true, watrue : the true cosmology
	"""
	#Note: the Eratio function assumes a flat universe
	cosmo={testVars[i]:testPoint[i] for i in range(len(testVars))}
	if 'w' not in testVars:
		if 'w0' not in cosmo.keys():
			cosmo['w0']=w0true
		if 'wa' not in cosmo.keys():
			cosmo['wa']=watrue
		true_cosmo={'Om0':Om0,'w0':w0true,'wa':watrue,'h':htrue}
	else:
		true_cosmo={'Om0':Om0,'w':w0true,'h':htrue}
	for k in set_cosmo.keys():
		if k not in cosmo.keys():
			cosmo[k]=set_cosmo[k]
	if 'Om0' not in cosmo.keys():
		if 'Ode0' in cosmo.keys():
			cosmo['Om0']=1-cosmo['Ode0']
			cosmo.pop('Ode0')
		else:
			cosmo['Om0']=Om0	
	if 'h' not in cosmo.keys():
		cosmo['h']=htrue
	
	Eratio_true = Eratio(zl, zs,true_cosmo )
	Eratio_test = Eratio(zl, zs, cosmo)

	#if np.isfinite(Eratio_w0wa):
	loglike = -0.5 * ( (Eratio_true-Eratio_test)**2 / 
					  Evariance(Eratio_test, dTc))
	#else:
	#    loglike = -np.inf
	return(loglike)


def rescale_likelihood( a ):
	"""
	Rescale the likelihood array using 
	a Sorted Cumulative Sum function (e.g., for making contours)
	Construct an array "sumabove" such that the cell 
	at index i in sumabove is equal to the sum of all 
	cells from the input array "a" that have a
	cell value higher than a[i]
	"""
	# Collapse the array into 1 dimension
	sumabove = deepcopy(a).ravel()

	# Sort the raveled array by descending cell value
	iravelsorted = sumabove.argsort( axis=0 )[::-1]

	# Reassign each cell to be the cumulative sum of all
	# input array cells with a higher value :
	sumabove[iravelsorted] = sumabove[iravelsorted].cumsum()

	# Normalize by the max value in the array
	sumabove /= sumabove.max()

	# Now unravel back into shape of original array and return
	return( sumabove.reshape( a.shape ) )



class Survey(object):
	"""
	A Survey class that enables cosmology tests with assumptions of
	redshift distributions and time delay/lens model precision.

	Parameters
	----------
	N: int
		The number of discovered glSN (overruled by zl/zs below)
	dTL: float
		The percent precision on each lens model.
	dTT: float
		The percent precision on each time delay measurement
	zl: float or list
		Redshift(s) of the lens(es). If a float, assumes you want N identical SN.
	zs: float or list
		Redshift(s) of the source(s). If a float, assumes you want N identical SN.
	name: str
		Name of your survey

	"""
	def __init__(self, N=10,dTL=5, dTT=5, zl=0.3, zs=0.8, name='mySurvey', **kwargs):
		if not isinstance(zl,(list,tuple,np.ndarray)):
			zl=[zl]
		if not isinstance(zs,(list,tuple,np.ndarray)):
			zs=[zs]
		assert(len(zl)==len(zs))

		if len(zl)>1:
			self.N = len(zl)
		else:
			self.N = N
		self.dTL = dTL
		self.dTT = dTT
		self.zl = zl
		self.zs = zs
		
		self.grid_likelihood = None
		self.nestle_result = None
		self.w0 = -1
		self.wa = 0
		self.Om0 = 0.3
		self.w = -1
		self.h = .7
		self.cosmo_truths={'h':self.h,'w':self.w,'Ode0':1-self.Om0,
					'w0':self.w0,'wa':self.wa,'Om0':self.Om0}
		self.name=name

	@property
	def dTc(self):
		return np.sqrt(self.dTL**2 + self.dTT**2) / np.sqrt(self.N)   

	def survey_grid(self,vparam_names,bounds,npoints=100,grad_param=None,
			grad_param_bounds=None,ngrad=10,**kwargs):
		if len(vparam_names)!=2:
			print('For grid mode, must provide exactly 2 parameters.')
			return
		for p in vparam_names:
			if p not in bounds.keys():
				print('Must provide bounds for every parameter.')
				return

		param1_list = np.linspace(bounds[vparam_names[0]][0], bounds[vparam_names[0]][1], npoints)
		param2_list = np.linspace(bounds[vparam_names[1]][0], bounds[vparam_names[1]][1], npoints)
		self.param1_list=param1_list
		self.param2_list=param2_list

		p1grid, p2grid = np.meshgrid(param1_list, param2_list) 

		if grad_param is None:
			loglE = np.array([[loglikelihoodE(vparam_names,[param1_list[i],param2_list[j]], 
						self.zl, self.zs,  dTc=kwargs.get('dTc',self.dTc), 
					    Om0=self.cosmo_truths['Om0'], w0true=self.cosmo_truths['w0'], 
					    watrue=self.cosmo_truths['wa'],htrue=self.cosmo_truths['h']) for i in range(len(param1_list))] for j in range(len(param2_list))])[:,:,0]
							   
			likelihood = np.exp(loglE)
			like_rescaled = rescale_likelihood(likelihood)
			self.grid_likelihood = like_rescaled
			
		else:

			all_res={}
			grad_param_range=np.linspace(grad_param_bounds[0],grad_param_bounds[1],ngrad)
			for g in grad_param_range:
				loglE = np.array([[loglikelihoodE(vparam_names,[param1_list[i],param2_list[j]], 
						self.zl, self.zs,  dTc=kwargs.get('dTc',self.dTc), set_cosmo={grad_param:g},
					    Om0=self.cosmo_truths['Om0'], w0true=self.cosmo_truths['w0'], 
					    watrue=self.cosmo_truths['wa'],htrue=self.cosmo_truths['h']) for i in range(len(param1_list))] for j in range(len(param2_list))])[:,:,0]
							   
				likelihood = np.exp(loglE)
				like_rescaled = rescale_likelihood(likelihood)
				all_res[g] = like_rescaled
				
			self.grid_grad_res=all_res
			self.grid_grad_param=grad_param
			self.grid_grad_param_range=grad_param_range
			self.grid_likelihood = np.median([all_res[g] for g in grad_param_range])
		self.grid_samples = [p1grid,p2grid]
		self.grid_vparam_names=vparam_names

	def survey_nestle(self,vparam_names,bounds,constants={},npoints=100, method='single',
			maxiter=None, maxcall=None, modelcov=False, rstate=None,
			verbose=False, warn=True, **kwargs):
		
		"""
		MCMC-like measurement of precision on cosmological parameters
		based on your survey ensemble of glSN and precsions.
		"""
		ppfs={}
		if bounds is not None:
			for key, val in bounds.items():
				a, b = val
				f = sncosmo.utils.Interp1D(0., 1., np.array([a, b]))
				ppfs[key] = f

		
		
		iparam_names = [key for key in vparam_names if key in ppfs]
		
		ppflist = [ppfs[key] for key in iparam_names]
		npdim = len(iparam_names)  # length of u
		ndim = len(vparam_names)  # length of v

		def prior_transform(u):
			d = {}
			for i in range(npdim):
				d[iparam_names[i]] = ppflist[i](u[i])
			v = np.empty(ndim, dtype=np.float)
			for i in range(ndim):
				key = vparam_names[i]
				v[i] = d[key]
			return v

		
		def likelihood(parameters):
			newParams={}
			newParams={iparam_names[i]:parameters[i] for i in range(len(parameters))}
			for param in constants.keys():
				newParams[param]=constants[param]

			return(np.sum(loglikelihoodE(vparam_names,[newParams[p] for p in iparam_names], 
					self.zl, self.zs, dTc=kwargs.get('dTc',self.dTc), 
					   Om0=self.cosmo_truths['Om0'], w0true=self.cosmo_truths['w0'], 
					   watrue=self.cosmo_truths['wa'],htrue=self.cosmo_truths['h'])))
		res = nestle.sample(likelihood, prior_transform, ndim, npdim=npdim,
						npoints=npoints, method=method, maxiter=maxiter,
						maxcall=maxcall, rstate=rstate,#decline_factor=.5,
						callback=(nestle.print_progress if verbose else None))

		self.nestle_result=res
		self.nestle_cosmology_fit={p:weighted_quantile(res.samples[:,iparam_names.index(p)],[.16,.5,.84],res.weights) for p in iparam_names}
		self.nestle_vparam_names=vparam_names
	def plot_survey_gradient(self,math_labels=None):
		if self.grid_likelihood is None:
			print('Cannot plot survey gradient without running grid first.')
			return	
		
		if math_labels is None:
			math_labels=self.grid_vparam_names

		ax=plot('scatter',[],[],x_lab=math_labels[0],y_lab=math_labels[1])
		ind=-1
		contours=[]
		for g in self.grid_grad_param_range:

			contours.append(ax.contour(self.grid_samples[0],self.grid_samples[1],self.grid_grad_res[g],levels=[.4]))
			ind+=1
			strs=['%.2f'%g]
			fmt={}
			for l, s in zip(contours[ind].levels, strs):
			    fmt[l] = s

			ax.clabel(contours[ind],contours[ind].levels,inline=True,fmt=fmt,fontsize=10)

			
		plt.show()

		
			
			
	def plot_survey_contour(self,params=None,math_labels=None):
		if self.nestle_result is None and self.grid_likelihood is None:
			print('Cannot plot nestle result without running nestle or grid first.')
			return

		if self.nestle_result is not None:
			if params is None:
				params=self.nestle_vparam_names
			if math_labels is None:
				math_labels=params
			param_inds=[]
			for p in params:
				param_inds.append(self.nestle_vparam_names.index(p))
			fig=corner.corner(self.nestle_result.samples[:,param_inds],weights=self.nestle_result.weights,
				labels=math_labels,
				truths=[self.cosmo_truths[p] for p in params],
				no_fill_contours=True,fill_contours=False,levels=[.68])
		else:
			fig = plt.figure(figsize=[6,4])
			ax = fig.add_subplot(1,1,1)
			if params is None:
				params=self.grid_vparam_names
			if math_labels is None:
				math_labels=params

			ax.contourf(self.grid_samples[0], self.grid_samples[1], self.grid_likelihood, colors='#1f77b4',
					levels=[0,0.68], alpha=.9, zorder=10)
			ax.contourf(self.grid_samples[0], self.grid_samples[1], self.grid_likelihood, colors='#1f77b4',
					levels=[0,0.95], alpha=.3, zorder=10)

			ax.plot(self.cosmo_truths[params[0]],self.cosmo_truths[params[1]],marker='o', mec='k', mfc='w', ms=6)
			ax.axhline(self.cosmo_truths[params[1]],ls='--', color='0.6', lw=0.8)
			ax.axvline(self.cosmo_truths[params[0]],ls='--', color='0.6', lw=0.8)
			ax.set_xlabel(math_labels[0], fontsize=20)
			ax.set_ylabel(math_labels[1], fontsize=20)
			plt.setp(ax.xaxis.get_ticklabels(), fontsize='large')
			plt.setp(ax.yaxis.get_ticklabels(), fontsize='large')
			plt.tight_layout()



