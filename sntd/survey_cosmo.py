from matplotlib import pyplot as plt
import numpy as np
import scipy
from scipy.integrate import quad
from matplotlib import ticker, rcParams
from matplotlib.lines import Line2D
from scipy.misc import derivative as deriv
import pandas

from copy import deepcopy,copy
import sncosmo
import nestle
import corner


from .coetools import *
from .util import *



def EA1(z1, z2, cosmo):
	"""The integral of the inverse of the normalized 
	Hubble parameter, from z1 to z2:  int_z1^z2{1/E(z')dz}
	Eq 8 of Coe & Moustakas 2009 (for non-flat)"""
	
	if 'w' in cosmo.keys():
		Ez2 = lambda z: ( cosmo['Om0']*(1+z)**3 + cosmo['Ok']*(1+z)**2+
					 cosmo['Ode0']*(1+z)**(3*(1+cosmo['w'])))
	else:
		Ez2 = lambda z: ( cosmo['Om0']*(1+z)**3 + cosmo['Ok']*(1+z)**2 + 
					 cosmo['Ode0']*(1+z)**(3*(1+cosmo['w0']+cosmo['wa'])) * 
					 np.exp(-3*cosmo['wa']*z/(1+z)) )
	Ezinv = lambda z: 1 / np.sqrt(Ez2(z))
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
		

	return(result[:n1],result[n1:])

def EA2(z1, z2, cosmo):
	"""The integral of the inverse of the normalized 
	Hubble parameter, from z1 to z2:  int_z1^z2{1/E(z')dz}
	Eq 8 of Coe & Moustakas 2009. Assumes flat universe. """

	if 'w' in cosmo.keys():
		Ez2 = lambda z: ( cosmo['Om0']*(1+z)**3 + cosmo['Ok']*(1+z)**2+
					 cosmo['Ode0']*(1+z)**(3*(1+cosmo['w'])))
	else:
		Ez2 = lambda z: ( cosmo['Om0']*(1+z)**3 + cosmo['Ok']*(1+z)**2 + 
					 cosmo['Ode0']*(1+z)**(3*(1+cosmo['w0']+cosmo['wa'])) * 
					 np.exp(-3*cosmo['wa']*z/(1+z)) )
	Ezinv = lambda z: 1 / np.sqrt(Ez2(z))
	# if Ez2(z1)>0 and Ez2(z2)>0:
	return(np.array([quad(Ezinv,z1[i],z2[i])[0] for i in range(len(z1))]))



def Eratio(zl, zs, cosmo):
	"""The time delay expansion function ratio 
	(script E, eq 12 in Coe & Moustakas 2009)."""

	EL,ES = EA1(zl,zs, cosmo)
	if cosmo['Ok']<-.0001:
		ELS=EA2(zl,zs,cosmo)
		EL=np.sin(np.sqrt(np.abs(cosmo['Ok']))*EL)/np.sqrt(np.abs(cosmo['Ok']))
		ES=np.sin(np.sqrt(np.abs(cosmo['Ok']))*ES)/np.sqrt(np.abs(cosmo['Ok']))
		ELS=np.sin(np.sqrt(np.abs(cosmo['Ok']))*ELS)/np.sqrt(np.abs(cosmo['Ok']))
	elif cosmo['Ok']>.0001:
		ELS=EA2(zl,zs,cosmo)
		EL=np.sinh(np.sqrt(np.abs(cosmo['Ok']))*EL)/np.sqrt(np.abs(cosmo['Ok']))
		ES=np.sinh(np.sqrt(np.abs(cosmo['Ok']))*ES)/np.sqrt(np.abs(cosmo['Ok']))
		ELS=np.sinh(np.sqrt(np.abs(cosmo['Ok']))*ELS)/np.sqrt(np.abs(cosmo['Ok']))
	else:# np.abs(cosmo['Ok'])<.001:
		ELS = ES-EL
	
	

	Erat = EL * ES / ELS
	return(Erat/(cosmo['h']*100))

def Evariance(Erat, dTc):
	"""The variance on the E ratio due to the 
	time delay and lens potential uncertainty. 
	dTc : percent uncertainty on time delay distance ratio 
		(combining time delay measurement + lens modeling error)  
		(This is a percentage:  enter "5" for 5% precision)
	"""
	return( Erat**2 * (dTc/100.)**2 )


def loglikelihoodE(testVars,testPoint, zl, zs, dTc=2, set_cosmo={},Eratio_true=None,
				   Om0true=0.3,Oktrue=0, Ode0true=.7,w0true=-1., watrue=0.,htrue=.7,
				   return_ratios=False,P=1): 
	"""log(likelihood) for a time delay distance 
	ratio constraint, comparing 
	a given test point to the true position. 
	
	testVars: the test variables
	testPOint: the test point
	zl : lens redshift(s)
	zs : source redshift(s)
	dTc : percent precision on cosmological distance ratio, combining
	   uncertainty from the time delay measurement + lens modeling
	set_cosmo: dictionary of constants to set (if not defaults)
	Eratio_true: The true ratio if computed beforehand for speed
	Om0true: True value of Om0
	Oktrue: True value of Ok
	Ode0true: True value of Ode0
	w0true: True value of w0
	watrue: True value of wa
	htrue: True value of h
	return_ratios: if true, return the actual ratios in additon to the likelihood
	P: The probability of each source/lens redshift combo (see C&M 2009 equation 17)
	"""

	cosmo={testVars[i]:testPoint[i] for i in range(len(testVars))}
	if 'w' not in testVars:
		if 'w0' not in cosmo.keys():
			cosmo['w0']=w0true
		if 'wa' not in cosmo.keys():
			cosmo['wa']=watrue
		true_cosmo={'Om0':Om0true,'Ok':Oktrue,'Ode0':Ode0true,'w0':w0true,'wa':watrue,'h':htrue}
	else:
		true_cosmo={'Om0':Om0true,'Ok':Oktrue,'Ode0':Ode0true,'w':w0true,'h':htrue}
	for k in set_cosmo.keys():
		if k not in cosmo.keys():
			cosmo[k]=set_cosmo[k]

	if 'Om0' not in cosmo.keys():
		if 'Ode0' in cosmo.keys() and 'Ok' in cosmo.keys():
			cosmo['Om0']=1-cosmo['Ode0']-cosmo['Ok']
			#cosmo.pop('Ode0')
		else:
			cosmo['Om0']=Om0true	
	if 'Ode0' not in cosmo.keys():
		if 'Ok' in cosmo.keys():
			cosmo['Ode0']=1-cosmo['Om0']-cosmo['Ok']
		else:
			cosmo['Ode0']=Ode0true
	if 'Ok' not in cosmo.keys():
		cosmo['Ok']=1-cosmo['Om0']-cosmo['Ode0']
	if 'h' not in cosmo.keys():
		cosmo['h']=htrue


	if Eratio_true is None:
		Eratio_true = Eratio(zl, zs,true_cosmo )
	Eratio_test = Eratio(zl, zs, cosmo)

	
	#if np.isfinite(Eratio_w0wa):
	loglike = -0.5 * np.sum( P*(Eratio_true-Eratio_test)**2 / 
					  Evariance(Eratio_test, dTc))
	if return_ratios:
		return(Eratio_true,Eratio_test,loglike)
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

def deriv_like(p1,p2,names,zl,zs,dTc,Om0,Ok,Ode0,w0true,watrue,htrue,P):
	"""
	helper to calculate fisher matrix
	"""
	if p2 is not None:
		return -2*loglikelihoodE(names,[p1,p2],zl,zs,dTc=dTc,Om0true=Om0,Oktrue=Ok,Ode0true=Ode0,
					w0true=w0true,watrue=watrue,htrue=htrue,P=P)
	else:
		return -2*loglikelihoodE(names,[p1],zl,zs,dTc=dTc,Om0true=Om0,Oktrue=Ok,
					Ode0true=Ode0,w0true=w0true,watrue=watrue,htrue=htrue,P=P)


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
	P: float or list
		The probability of each source/lens redshift combo, defaults to calculating
		assuming gaussian redshift distributions
	name: str
		Name of your survey

	"""
	def __init__(self, N=10,dTL=5, dTT=5, zl=0.3, zs=0.8,P=None, name='mySurvey', **kwargs):
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
		if P is not None:
			self.P=P if isinstance(P,(list,tuple,np.ndarray)) else [P]*len(zl)
		else:
			zl_p=scipy.stats.norm(np.mean(zl),np.std(zl))
			zs_p=scipy.stats.norm(np.mean(zs),np.std(zs))			

			
			self.P = np.array([scipy.integrate.simps(zl_p.pdf([z-.01,z+.01]),[z-.01,z+.01]) for z in zl])*\
						np.array([scipy.integrate.simps(zs_p.pdf([z-.01,z+.01]),[z-.01,z+.01]) for z in zs])
		
		self.grid_likelihood = None
		self.nestle_result = None
		self.w0 = -1
		self.wa = 0
		self.Om0 = 0.3
		self.Ok = 0
		self.Ode0 = 0.7
		self.w = -1
		self.h = .7
		self.cosmo_truths={'h':self.h,'w':self.w,'Ode0':self.Ode0,'Ok':self.Ok,
					'w0':self.w0,'wa':self.wa,'Om0':self.Om0}
		self.name=name

	@property
	def dTc(self):
		"""
		Returns the combined lens/delay uncertainty assuming statistical uncertainties
		"""
		return np.sqrt(self.dTL**2 + self.dTT**2) / np.sqrt(self.N)   

	def survey_grid(self,vparam_names,bounds,npoints=100,grad_param=None,constants={},
			grad_param_bounds=None,ngrad=10,**kwargs):
		"""Calculate cosmological contours by varying 2 parameters in a grid.

		Parameters
		-------
		vparam_names: list
			The names of parameters to vary
		bounds: dict
			Dictionary with param names as keys and list/tuple/array of bounds as values
		npoints: int
			The number of sample points
		grad_param: str
			Parameter to assume we're to have measured wrong (see C&M 2009 Figure 4)
		constants: dict
			Constants that are not the defaults
		grad_param_bounds: dict
			Bounds for grad_param, same format as bounds
		ngrad: int
			Number of grid points to vary grad_param

		Returns
		-------
		Adds to class attribute "grid" (a dictionary), with a comma-separated list of 
		the vparam_names as the key and grid values as the value. 
		"""

		if len(vparam_names)!=2:
			print('For grid mode, must provide exactly 2 parameters.')
			return
		for p in vparam_names:
			if p not in bounds.keys():
				print('Must provide bounds for every parameter.')
				return

		param1_list = np.linspace(bounds[vparam_names[0]][0], bounds[vparam_names[0]][1], npoints)
		param2_list = np.linspace(bounds[vparam_names[1]][0], bounds[vparam_names[1]][1], npoints)

		

		p1grid, p2grid = np.meshgrid(param1_list, param2_list) 

		true_ratio=Eratio(self.zl, self.zs,self.cosmo_truths)
		
		if grad_param is None:
			loglE = np.array([[loglikelihoodE(vparam_names,[param1_list[i],param2_list[j]], 
						self.zl, self.zs,  dTc=kwargs.get('dTc',self.dTc), set_cosmo=constants,
						Om0true=self.cosmo_truths['Om0'],Oktrue=self.cosmo_truths['Ok'],Eratio_true=true_ratio,
						Ode0true=self.cosmo_truths['Ode0'], w0true=self.cosmo_truths['w0'], 
						watrue=self.cosmo_truths['wa'],htrue=self.cosmo_truths['h']) for i in range(len(param1_list))] for j in range(len(param2_list))])[:,:]
			likelihood = np.exp(loglE)
			likelihood[np.isnan(likelihood)]=0

			like_rescaled = rescale_likelihood(likelihood)
			if self.grid_likelihood is None:
				self.grid_likelihood = {}
				self.grid_samples = {}
				self.param1_list = {}
				self.param2_list = {}

			self.grid_likelihood[','.join(vparam_names)] = like_rescaled
			
		else:

			all_res={}
			grad_param_range=np.linspace(grad_param_bounds[0],grad_param_bounds[1],ngrad)
			for g in grad_param_range:
				loglE = np.array([[loglikelihoodE(vparam_names,[param1_list[i],param2_list[j]], 
						self.zl, self.zs,  dTc=kwargs.get('dTc',self.dTc), set_cosmo={grad_param:g},
						Om0true=self.cosmo_truths['Om0'],Oktrue=self.cosmo_truths['Ok'],
						Ode0true=self.cosmo_truths['Ode0'], w0true=self.cosmo_truths['w0'], 
						watrue=self.cosmo_truths['wa'],htrue=self.cosmo_truths['h']) for i in range(len(param1_list))] for j in range(len(param2_list))])[:,:,0]
							   
				likelihood = np.exp(loglE)
				likelihood[np.isnan(likelihood)]=0
				like_rescaled = rescale_likelihood(likelihood)
				all_res[g] = like_rescaled
				
			self.grid_grad_res=all_res
			self.grid_grad_param=grad_param
			self.grid_grad_param_range=grad_param_range
			if self.grid_likelihood is None:
				self.param1_list = {}
				self.param2_list = {}
				self.grid_samples = {}
				self.grid_likelihood = {}
				self.grid_grad_res={}
				self.grid_grad_param={}
				self.grid_grad_param_range={}

			self.grid_grad_res[','.join(vparam_names)]=all_res
			self.grid_grad_param[','.join(vparam_names)]=grad_param
			self.grid_grad_param_range[','.join(vparam_names)]=grad_param_range
			self.grid_likelihood[','.join(vparam_names)] = np.median([all_res[g] for g in grad_param_range])
		self.param1_list[','.join(vparam_names)]=param1_list
		self.param2_list[','.join(vparam_names)]=param2_list
		self.grid_samples[','.join(vparam_names)] = [p1grid,p2grid]
		

	def survey_nestle(self,vparam_names,bounds,constants={},npoints=100, **kwargs):
		
		"""Calculate cosmological contours in an MCMC-like fashion.

		Parameters
		-------
		vparam_names: list
			The names of parameters to vary
		bounds: dict
			Dictionary with param names as keys and list/tuple/array of bounds as values
		npoints: int
			The number of sample points
		grad_param: str
			Parameter to assume we're to have measured wrong (see C&M 2009 Figure 4)
		constants: dict
			Constants that are not the defaults
		grad_param_bounds: dict
			Bounds for grad_param, same format as bounds
		ngrad: int
			Number of grid points to vary grad_param

		Returns
		-------
		Adds to class attribute "grid" (a dictionary), with a comma-separated list of 
		the vparam_names as the key and grid values as the value. 
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

			return(loglikelihoodE(vparam_names,[newParams[p] for p in iparam_names], 
					self.zl, self.zs, dTc=kwargs.get('dTc',self.dTc), 
					   Om0true=self.cosmo_truths['Om0'],Oktrue=self.cosmo_truths['Ok'],
						Ode0true=self.cosmo_truths['Ode0'], w0true=self.cosmo_truths['w0'], 
					   watrue=self.cosmo_truths['wa'],htrue=self.cosmo_truths['h']))
		res = nestle.sample(likelihood, prior_transform, ndim, npdim=npdim,
						npoints=npoints, method=kwargs.get('method'), maxiter=kwargs.get('maxiter',None),
						maxcall=kwargs.get('maxcall',None), rstate=kwargs.get('rstate',None),
						callback=(nestle.print_progress if kwargs.get('verbose',False) else None),**kwargs)
		if self.nestle_result is None:
			self.nestle_result = {}
			self.nestle_cosmology_fit = {}

		self.nestle_result[','.join(vparam_names)] = res
		self.nestle_cosmology_fit[','.join(vparam_names)] = {p:weighted_quantile(res.samples[:,iparam_names.index(p)],[.16,.5,.84],res.weights) for p in iparam_names}
		
	

	def survey_fisher(self,params,dx=1e-6):
		def take_deriv(p2,p2name,p1,p1name):
			return deriv(deriv_like,p1,dx=dx,
					args=(p2,[p1name,p2name], self.zl, self.zs,self.dTc,
				   self.cosmo_truths['Om0'],self.cosmo_truths['Ok'],
				   self.cosmo_truths['Ode0'], self.cosmo_truths['w0'], self.cosmo_truths['wa'],
								self.cosmo_truths['h'],self.P))
		
		fisher_matrix=np.zeros((len(params),len(params)))
		for i in range(len(params)):
			for j in range(len(params)):
				if i==j:
					fisher_matrix[i][j]=.5*np.abs(deriv(deriv_like,self.cosmo_truths[params[i]],dx=dx,
						args=(None,[params[i]], self.zl, self.zs, self.dTc,
					   self.cosmo_truths['Om0'],self.cosmo_truths['Ok'],self.cosmo_truths['Ode0'],
					   self.cosmo_truths['w0'], self.cosmo_truths['wa'],
									self.cosmo_truths['h'],self.P),n=2))
				else:
					fisher_matrix[i][j]=.5*deriv(take_deriv,self.cosmo_truths[params[i]],dx=dx,
						args=(params[i],self.cosmo_truths[params[j]],params[j]))
					
		
		#fisher_matrix=np.array([[49824.9224,-1829.7018,-4434.2995,4546.8899,122.5319],
		#						[-1829.7018, 88.3760, 200.9795, -189.2658, -8.4386],
		#						[-4434.2995, 200.9795, 463.5732, -445.5690, -17.9694],
		#						[4546.8899, -189.2658, -445.5690, 441.9725, 15.2981],
		#						[122.5319, -8.4386, -17.9694, 15.2981, 1.0394]])
		self.fisher_matrix=Fisher(data=fisher_matrix,params=params,name=self.name,cosmo_truths=self.cosmo_truths)
		print(self.fisher_matrix.pretty_fish)

	def plot_survey_gradient(self,params, math_labels=None):
		if self.grid_likelihood is None or (','.join(params) not in self.grid_likelihood.keys() and ','.join([params[1],params[0]]) not in self.grid_likelihood.keys()):
			print('Cannot plot survey gradient without running grid first.')
			return	
		
		if ','.join(params) not in self.grid_likelihood.keys():
			params=[params[1],params[0]]
		if math_labels is None:
			math_labels=self.grid_vparam_names

		ax=plot('scatter',[],[],x_lab=math_labels[0],y_lab=math_labels[1])
		ind=-1
		contours=[]
		for g in self.grid_grad_param_range:

			contours.append(ax.contour(self.grid_samples[','.join(params)][0],self.grid_samples[','.join(params)][1],
						self.grid_grad_res[','.join(params)][g],levels=[.4]))
			ind+=1
			strs=['%.2f'%g]
			fmt={}
			for l, s in zip(contours[ind].levels, strs):
				fmt[l] = s

			ax.clabel(contours[ind],contours[ind].levels,inline=True,fmt=fmt,fontsize=10)

			
		plt.show()

			
			
	def plot_survey_contour(self,params,math_labels=None,color='#1f77b4',filled=True,confidence=[.68,.95],
				fom=False,ax=None,alphas=[.9,.3],**kwargs):
		if self.nestle_result is None and self.grid_likelihood is None \
			or (self.nestle_result is not None and (','.join(params) not in self.nestle_result.keys() and ','.join([params[1],params[0]]) not in self.nestle_result.keys()) and\
			 (self.grid_likelihood is not None and ','.join(params) not in self.grid_likelihood.keys() and ','.join([params[1],params[0]]) not in self.grid_likelihood.keys())):
			print('Cannot plot nestle result without running nestle or grid first.')
			return

		if self.nestle_result is not None:
			if ','.join(params) not in self.nestle_result.keys():
				params=[params[1],params[0]]
			if math_labels is None:
				math_labels=params
			param_inds=[]
			param_inds=[0,1]
			fig=corner.corner(self.nestle_result[','.join(params)].samples[:,param_inds],weights=self.nestle_result[','.join(params)].weights,
				labels=math_labels,
				truths=[self.cosmo_truths[p] for p in params],
				no_fill_contours=filled,fill_contours=filled,levels=[.68])
		else:

			if ax is None:
				fig = plt.figure(figsize=[9,6])
				ax = fig.gca()#add_subplot(1,1,1)
				plt.setp(ax.xaxis.get_ticklabels(), fontsize='large')
				plt.setp(ax.yaxis.get_ticklabels(), fontsize='large')
				#plt.tight_layout()
			if ','.join(params) not in self.grid_likelihood.keys():
				params=[params[1],params[0]]
			if math_labels is None:
				math_labels=params
			if isinstance(confidence,(float,int)):
				confidence=[confidence]
			alphas=np.linspace(.9,.3,len(confidence)) if alphas is None else alphas
			if len(alphas)!=len(confidence):
				alphas=[alphas[0]]*len(confidence)
			for i in range(len(confidence)):
				if filled:
					CS=ax.contourf(self.grid_samples[','.join(params)][0], self.grid_samples[','.join(params)][1], self.grid_likelihood[','.join(params)], colors=color,
						levels=[0,np.sort(confidence)[i]], alpha=alphas[i], zorder=10,**kwargs)
					lines=Line2D([0],[0],color=color,alpha=alphas[0],linewidth=10)
				else:
					CS=ax.contour(self.grid_samples[','.join(params)][0], self.grid_samples[','.join(params)][1], self.grid_likelihood[','.join(params)], colors=color,
						levels=[0,np.sort(confidence)[i]],**kwargs)
					lines=Line2D([0],[0],color=color,alpha=alphas[0],linewidth=10,marker='.')

			if fom:
				CS954=ax.contour(self.grid_samples[','.join(params)][0], self.grid_samples[','.join(params)][1], self.grid_likelihood[','.join(params)],alpha=0,
						levels=[.954])

				contour = CS954.collections[0]
				vs = contour.get_paths()[0].vertices
				x=vs[:,0]
				y=vs[:,1]
				area=0.5*np.sum(y[:-1]*np.diff(x) - x[:-1]*np.diff(y))
				area=np.abs(area)
				line_name=self.name+': FOM=%.2f'%(np.pi/area)
			else:

				line_name=self.name
		

			ax.plot(self.cosmo_truths[params[0]],self.cosmo_truths[params[1]],marker='o', mec='k', mfc='w', ms=6)
			ax.axhline(self.cosmo_truths[params[1]],ls='--', color='0.6', lw=0.8)
			ax.axvline(self.cosmo_truths[params[0]],ls='--', color='0.6', lw=0.8)
			ax.set_xlabel(math_labels[0], fontsize=20)
			ax.set_ylabel(math_labels[1], fontsize=20)
			

			return(ax,lines,line_name)

class Fisher:
	def __init__(self, inroot='',xvar='', yvar='', fixes=[], margs=[],
				 data=[], params=[], silent=False,name='my_fisher',cosmo_truths={'h':.7,'w':0,'Ode0':.7,
					'w0':0,'wa':-1,'Om0':.3}):
		self.name=name
		self.xvar = xvar
		self.yvar = yvar
		self.fixes = fixes
		self.margs = margs
		self.data = data
		self.params = params
		self.silent = silent 
		self.nFish=1
		self.inroot=inroot
		self.cosmo_truths=cosmo_truths

		if len(self.params)>0 and len(self.data)>0:

			self.pretty_fish=pandas.DataFrame.from_dict({k:list(self.data[self.params.index(k)][:]) for k in self.params},
																	orient='index',columns=self.params)
		else:
			self.pretty_fish=None

		self.fish_list=[self]
	def load(self,fishdir=''):  # DETFast join
		txt = loadfile(self.inroot+'.fisher', dir=fishdir, silent=self.silent)
		nparam = int(txt[0].split()[1])
		self.params = []
		for ivar in range(nparam):
			param = txt[ivar+4].split()[0]
			self.params.append(param)

		self.data = loaddata(self.inroot+'.fisher+', dir=fishdir, headlines=4+nparam, silent=1)

	def ii(self):
		self.ix = None
		self.iy = None
		if self.xvar:
			self.ix = self.params.index(self.xvar)
		if self.yvar:
			self.iy = self.params.index(self.yvar)

		self.ifixes = []
		for i, param in enumerate(self.params):
			if param in self.fixes:
				self.ifixes.append(i)

		self.imargs = []
		for i, param in enumerate(self.params):
			if param in self.margs:
				self.imargs.append(i)

	def pindex(self, param):
		return self.params.index(param)

	def take(self, iparams):
		self.data = self.data.take(iparams, 0)
		self.data = self.data.take(iparams, 1)
		self.params = list(take(self.params, iparams))

	def reorder(self, params):
		"""Matrix with params in new order"""
		self.repar(params)
		iparams = list(map(self.pindex, params))
		self.take(iparams)

	def rename(self, pdict1=None):
		"""Rename parameters given a dictionary of names & nicknames"""
		pdict1 = pdict1 or pdict

		for i, param in enumerate(self.params):
			self.params[i] = pdict1.get(param, param)
		self.pretty_fish=pandas.DataFrame.from_dict({k:list(self.data[self.params.index(k)][:]) for k in self.params},
													 orient='index',columns=self.params)

	def fix(self, fixes=[]):
		"""Fix parameters constant <==> Remove them from the Fisher matrix"""
		self.fixes = fixes or self.fixes
		self.fixes = strspl(self.fixes)

		self.ii()

		iall = arange(len(self.params))
		ikeep = set(iall) - set(self.ifixes)  # Sorts result
		ikeep = list(ikeep)

		self.take(ikeep)

	def marg(self, margs=[]):
		"""Marginalize over variables: Remove them from the covariance matrix"""
		self.margs = margs or self.margs
		self.ii()

		#ikeep = invertselection(arange(len(self.params)), self.fixes)
		iall = arange(len(self.params))
		ikeep = set(iall) - set(self.imargs)  # Sorts result
		ikeep = list(ikeep)

		C = self.cov()
		C = C.take(ikeep, 0)
		C = C.take(ikeep, 1)
		self.data = inv(C)

		self.params = list(take(self.params, ikeep))

	def transform(self, params, M):
		"""Transform to new set of parameters using matrix provided"""
		self.data = matrix_multiply([transpose(M), self.data, M])
		self.params = params

	def prior(self, param, sig):
		"""Set a prior of sig on param"""
		ivar = self.pindex(param)
		self.data[ivar,ivar] = self.data[ivar,ivar] + 1 / sig**2


	def cov(self):
		"""Covariance matrix"""
		return inv(self.data)

	def dxdyp(self, xvar='', yvar=''):  # , fixes=None
		"""Return uncertainty in two parameters and their correlation"""
		self.xvar = xvar or self.xvar
		self.yvar = yvar or self.yvar
		#self.fixes = strspl(fixes or self.fixes)
		self.ii()

		C = self.cov()
		C = C.take((self.ix,self.iy),0)
		C = C.take((self.ix,self.iy),1)
		dx = np.sqrt(C[0,0])
		dy = np.sqrt(C[1,1])
		dxy = C[0,1]
		p = dxy / (dx * dy)
		self.C = C
		return dx, dy, p

	def dx(self, xvar=''): # , fixes=None
		"""Return uncertainty in parameter (if marginalizing over others)"""
		self.xvar = xvar or self.xvar
		#self.fixes = strspl(fixes or self.fixes)
		self.ii()

		#dx = 1 / sqrt(self.data[self.ix,self.ix])
		self.C = C = self.cov()
		dx = sqrt(C[self.ix,self.ix])
		return dx

	def addpar(self, param):
		npar = len(self.params)
		data = np.zeros((npar+1, npar+1))
		data[:npar,:npar] = self.data
		self.data = data
		self.params.append(param)

	def repar(self, params):
		for param in params:
			if param not in self.params:
				self.addpar(param)

	def pr(self):
		"""Print contents"""
		print(self.params)
		pint(self.data)

	def __add__(self, sel2):
		"""Add Fisher matrices"""
		pl1 = self.params
		pl2 = sel2.params
		n1 = len(pl1)
		n2 = len(pl2)
		F1 = self.data
		F2 = sel2.data
		if pl1 == pl2:
			pl = pl1
		else:
			params1 = set(pl1)
			params2 = set(pl2)
			params = params1 | params2
			pl = list(params)
		#pl = list(sort(pl))
		n = len(pl)

		# Put F1 in merged F join
		ii = []
		for p in pl1:
			i = pl.index(p)
			ii.append(i)

		FF1 = np.zeros((n,n))
		for i1, i in enumerate(ii):
			FF1[i].put(ii, F1[i1])

		# Put F2 in merged F join
		ii = []
		for p in pl2:
			i = pl.index(p)
			ii.append(i)

		FF2 = np.zeros((n,n))
		for i2, i in enumerate(ii):
			FF2[i].put(ii, F2[i2])
		# Add
		new = Fisher()
		new.data = FF1 + FF2
		new.params = pl
		new.xvar = self.xvar
		new.yvar = self.yvar
		new.fixes = self.fixes
		new.name=self.name+'+'+sel2.name
		new.pretty_fish=pandas.DataFrame.from_dict({k:list(self.data[self.params.index(k)][:]) for k in self.params},
													 orient='index',columns=self.params)
		new.nFish=self.nFish+sel2.nFish
		new.fish_list=[new,sel2]
		new.fish_list=np.append(new.fish_list,self.fish_list)
		new.fish_list=[x for x in new.fish_list if x.nFish==1 or x.nFish==new.nFish]
		return new

	def __mul__(self, fac):
		"""Multiply Fisher matrix by some factor"""
		new = Fisher()
		new.data = self.data * fac
		new.params = self.params
		new.xvar = self.xvar
		new.yvar = self.yvar
		new.fixes = self.fixes
		new.pretty_fish=pandas.DataFrame.from_items([(self.params[i],list(self.data[i][:])) for i in range(len(self.params))],
													 orient='index',columns=self.params)
		return new

	def __str__(self):
		print(self.pretty_fish)
		return('')

	def merit(self,param1,param2):

		# cov_matrix=self.cov()
		# ind1=self.pindex(param1)
		# ind2=self.pindex(param2)
		# sigma_x=np.sqrt(cov_matrix[ind1][ind1])
		# sigma_y=np.sqrt(cov_matrix[ind2][ind2])
		# rho=cov_matrix[ind1][ind2]
		# return(1./(sigma_x*sigma_y*np.sqrt(1-rho**2)))
		dx,dy,p=self.dxdyp(param1,param2)
		a,b,_=setell(dx,dy,p)
		return(1./(6.17*a*b))


	def plot(self,param1,param2,x_limits,y_limits,bestfit1=None,bestfit2=None,alpha = 0.9,color_list=None):
		xo=bestfit1 if bestfit1 is not None else self.cosmo_truths[param1]
		yo=bestfit2 if bestfit2 is not None else self.cosmo_truths[param2]
		if color_list is not None and not isinstance(color_list,(tuple,list)):
			color_list=[color_list]
		if color_list is None or len(color_list)!=len(self.fish_list):
			print('Either your color_list is the wrong size or you did not define it, taking defaults...')
			color_list=[blues,reds,purples,yellows,oranges,darkoranges,greens,greys,lightblues,
						lightreds,lightgreens,lightyellows,lightgreys]

		i=0
		patches=[]
		merits=[x.merit(param1,param2) for x in self.fish_list]

		for fish in np.array(self.fish_list)[np.argsort(merits)]:
			dx,dy,p=fish.dxdyp(param1,param2)
			print(dx,dy,p)
			plotellsp(xo, yo, dx, dy, p, colors=color_list[i], alpha=alpha)
			patches.append(plt.plot([],[],'s',ms=10,label=fish.name+': FOM=%.1f'%np.sort(merits)[i],color=color_list[i][0])[0])
			i+=1
		if x_limits is not None:
			plt.xlim(x_limits)
		if y_limits is not None:
			plt.ylim(y_limits)
		plt.xlabel(param1,fontsize=14)
		plt.ylabel(param2,fontsize=14)

		plt.figlegend(handles=patches,fontsize=14,bbox_to_anchor=(.75,.85))

def main():
	# xvar, yvar = strspl('w_0 w_a')
	# xo, yo = -1, 0  # Best fit w0, wa
	# alpha = 0.9
	# fish=Fisher(inroot='SN-IVS-o')
	# fish.load('/Users/jpierel/rodney/Fisher/data/')
	# print('xvar',fish.xvar)
	# print('yvar',fish.yvar)
	# print(fish.params)
	# dx, dy, p = fish.dxdyp(xvar, yvar)
	# print(dx,dy,p)
	# print(fish.data)
	# plotellsp(xo, yo, dx, dy, p,  alpha=alpha)
	# plt.xlim((-2,0))
	# plt.ylim((-5,5))
	# plt.show()
	# sys.exit()

	zl=np.random.normal(.5,.15,size=5000)
	zs=np.random.normal(2,.75,size=5000)
	goods=np.where(zs>zl)[0]
	wfirst_dict = {
	 			'N':10,   # number of Lensed SNe Ia with good time delays
	 			'dTL':5,  # % lens modeling uncertainty for each
	 			'dTT':1,  # % time delay measurement uncertainty for each
	 			'zl':zl[goods][:4000],'zs':zs[goods][:4000]
	 		}
	wfirst = Survey(**wfirst_dict)
	# #wfirst.survey_grid(['Om0','Ode0'],{'w':[-1.5,-.8],'w0':[-1.5,-.5],'wa':[-3,3],
	# 	#'h':[.62,.74],'Om0':[0,1],'Ode0':[0,1]},constants={'w':-1},dTc=.64,npoints=100)
	# #wfirst.plot_survey_contour()
	wfirst.survey_fisher(['h','Ode0','Ok','w0','wa'])
	wfirst.fisher_matrix.prior('Ode0',.00001)
	wfirst.fisher_matrix.prior('Ok',.00001)
	wfirst.fisher_matrix.prior('h',.00001)
	wfirst.fisher_matrix.plot('w0','wa',x_limits=(-1.6,-.4),y_limits=(-4,4))
	plt.show()

if __name__=='__main__':
	main()
