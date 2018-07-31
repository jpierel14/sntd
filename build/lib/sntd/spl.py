# This file is somehow the configuration of PyCS for the particular lens to be analysed.
# All the fine-tuning of methods happens here.
# This is done by defining the curve shifting functions that will be used in your scripts.

import pycs
import pycs.regdiff # (needs pymc)

# A simple attempt to get a multi-purpose free-knot spline method :
def spl(lcs):
        '''
        #best fit
	spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=150)
	for l in lcs:
		l.resetml()
        
	spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=100)
	
        spline = pycs.spl.topopt.opt_fine(lcs, nit=5, knotstep=50)
        '''
	spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
	for l in lcs:
		l.resetml()
		spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=60)
		spline = pycs.spl.topopt.opt_fine(lcs, nit=5, knotstep=100)
        
	return spline

# The dispersion method :
rawdispersionmethod = lambda lc1, lc2 : pycs.disp.disps.linintnp(lc1, lc2, interpdist = 30.0)
dispersionmethod = lambda lc1, lc2 : pycs.disp.disps.symmetrize(lc1, lc2, rawdispersionmethod)
def disp(lcs):
	return pycs.disp.topopt.opt_full(lcs, rawdispersionmethod, nit=10, verbose=False)

# The regression difference method : (needs pymc) : uncomment import on line 6 !
def regdiff(lcs):
	return pycs.regdiff.multiopt.opt_ts(lcs, pd=5, scale=200.0, verbose=False)


# The small scale extrinsic variability, used to generated the synthetic curves:
def Atweakml(lcs):
	return pycs.sim.twk.tweakml(lcs, beta=-1.5, sigma=0.25, fmin=1/500.0, fmax=None, psplot=False)

def Btweakml(lcs):
	return pycs.sim.twk.tweakml(lcs, beta=-1.0, sigma=0.9, fmin=1/500.0, fmax=None, psplot=False)

def Ctweakml(lcs):
	return pycs.sim.twk.tweakml(lcs, beta=-1.0, sigma=1.5, fmin=1/500.0, fmax=None, psplot=False)

def Dtweakml(lcs):
	return pycs.sim.twk.tweakml(lcs, beta=-0.0, sigma=4.5, fmin=1/500.0, fmax=None, psplot=False)
