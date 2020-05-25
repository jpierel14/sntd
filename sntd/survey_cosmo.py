from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import quad
from matplotlib import ticker, rcParams

from copy import deepcopy
import sncosmo
import nestle
import corner

from .util import weighted_quantile

def EA(z1, z2, Om0, w0, wa,iszero=False):
	"""The integral of the inverse of the normalized 
	Hubble parameter, from z1 to z2:  int_z1^z2{1/E(z')dz}
	Eq 8 of Coe & Moustakas 2009.
	Assuming a flat universe."""

	Ode0 = 1-Om0 # assuming flat universe
	#if w0>-0.3 or w0<-2:
	#    result = -np.inf
	#elif np.abs(wa)>5:
	#    result = -np.inf
	#else:
	Ez2 = lambda z: ( Om0*(1+z)**3 + 
					 Ode0*(1+z)**(3*(1+w0+wa)) * 
					 np.exp(-3*wa*z/(1+z)) )
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


def Eratio(zl, zs, Om0, w0, wa):
	"""The time delay expansion function ratio 
	(script E, eq 12 in Coe & Moustakas 2009).
	Assuming a flat universe."""

	EL,ES = EA(zl,zs, Om0, w0, wa,iszero=True)
	#ES = EA(0,zs, Om0, w0, wa,iszero=True)
	ELS = ES-EL#EA(zl, zs, Om0, w0, wa)
	
	#if EL<=0 or ES<=0 or ELS<=0:
	#    Erat = -np.inf
	#else:
	Erat = EL * ES / ELS
	return(Erat)

def Evariance(Erat, dTc):
	"""The variance on the E ratio due to the 
	time delay and lens potential uncertainty. 
	dTc : percent uncertainty on time delay distance ratio 
		(combining time delay measurement + lens modeling error)  
		(This is a percentage:  enter "5" for 5% precision)
	"""
	return( Erat**2 * (dTc/100.)**2 )


def loglikelihoodE(w0wa, zl, zs, dTc=0.1, 
				   Om0=0.3, w0true=-1., watrue=0.): 
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
	Eratio_true = Eratio(zl, zs, Om0, w0true, watrue)
	Eratio_w0wa = Eratio(zl, zs, Om0, w0wa[0], w0wa[1])

	#if np.isfinite(Eratio_w0wa):
	loglike = -0.5 * ( (Eratio_true-Eratio_w0wa)**2 / 
					  Evariance(Eratio_w0wa, dTc))
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
	def __init__(self, dTL=5, dTT=5, zl=0.3, zs=0.8, name='mySurvey', **kwargs):
		if not isinstance(zl,(list,tuple,np.ndarray)):
			zl=[zl]
		if not isinstance(zs,(list,tuple,np.ndarray)):
			zs=[zs]
		assert(len(zl)==len(zs))

		self.N = len(zl)
		self.dTL = dTL
		self.dTT = dTT
		self.zl = zl
		self.zs = zs
		
		self.likelihood = None
		self.nestle_result = None
		self.w0true = -1
		self.watrue = 0
		self.Om0true = 0.3
		self.name=name

	@property
	def dTc(self):
		return np.sqrt(self.dTL**2 + self.dTT**2) / np.sqrt(self.N)   
	
	def compute_likelihood(self, w0list, walist):
		loglE = np.array([[loglikelihoodE([w0, wa], self.zl, self.zs, self.dTc)
							for wo in w0list] 
							  for wa in walist])
		
		likelihood = 10**loglE
		like_rescaled = rescale_likelihood(likelihood)
		self.likelihood = like_rescaled

	def survey_nestle(self,vparam_names,bounds,constants={},
			npoints=100, method='single',
			maxiter=None, maxcall=None, modelcov=False, rstate=None,
			verbose=False, warn=True, **kwargs):
		

		ppfs={}
		if bounds is not None:
			for key, val in bounds.items():
				a, b = val
				f = sncosmo.utils.Interp1D(0., 1., np.array([a, b]))
				ppfs[key] = f

		required_params=['w0','wa']
		for param in required_params:
			if param not in ppfs.keys() and param not in constants.keys():
				raise RuntimeError('%s not in bounds or constants'%param)
		vparam_names=np.append(['w0','wa'],[x for x in vparam_names if x not in ['w0','wa']])
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
			for param in required_params:
				if param in vparam_names:
					newParams[param]=parameters[iparam_names.index(param)]
				else:
					newParams[param]=constants[param]

			return(np.sum(loglikelihoodE([newParams['w0'],newParams['wa']], 
					self.zl, self.zs, dTc=self.dTc, 
					   Om0=self.Om0true, w0true=self.w0true, watrue=self.watrue)))

		res = nestle.sample(likelihood, prior_transform, ndim, npdim=npdim,
							npoints=npoints, method=method, maxiter=maxiter,
							maxcall=maxcall, rstate=rstate,#decline_factor=.5,
							callback=(nestle.print_progress if verbose else None))

		self.nestle_result=res
		self.nestle_vparam_names=vparam_names
		self.cosmology_fit={p:weighted_quantile(res.samples[:,iparam_names.index(p)],[.16,.5,.84],res.weights) for p in iparam_names}


	def plot_survey_contour(self):
		if self.nestle_result is None:
			print('Cannot plot nestle result without running nestle first.')
			return
		
		fig=corner.corner(self.nestle_result.samples,weights=self.nestle_result.weights,
			labels=[r'$w_0$',r'$w_a$'],
			truths=[self.w0true,self.watrue],
			no_fill_contours=True,fill_contours=False,levels=[.5])

	




