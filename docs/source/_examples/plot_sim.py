"""
===================
Simulate Supernovae
===================

Simulating a multiply-imaged supernova.
"""
		
###############################################################
# Create a simulated multiply-imaged supernova that we can then fit,
# with no microlensing included in the simulation. Note that your final
# printed information will be different, as this is a randomly generated
# supernova. The function being used in these examples is 
# :py:func:`~sntd.simulation.createMultiplyImagedSN` . 
#
# **No Microlensing**
   
import sntd
import numpy as np

myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.4,z_lens=.53, bands=['F110W','F160W'],
             zp=[26.8,26.2], cadence=5., epochs=35.,time_delays=[20., 70.], magnifications=[10,5],
 objectName='My Type Ia SN',telescopename='HST',av_host=False)
print(myMISN)
myMISN.plot_object()


###############################################################
# Specify the distributions you want to use for any model
# parameter by providing a function that returns the parameter
# in any way you want. 

def x1_func():
    return(np.random.normal(1,.5))
def c_func():
    return(np.random.normal(-.05,.02))
param_funcs={'x1':x1_func,'c':c_func}
myMISN2 = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.33,z_lens=.53, bands=['F110W','F125W'],
              zp=[26.8,26.2], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[7,3.5],
              objectName='My Type Ia SN',telescopename='HST',sn_params=param_funcs)
print(myMISN2)
myMISN2.plot_object()


###############################################################
# Specify the distributions you want to use for dust
# parameters by providing a function that returns the parameter
# in any way you want. 

def hostav_func():
    return(np.random.normal(.5,.1))
def lensav_func():
    return(np.random.normal(.7,.2))
param_funcs={'host':hostav_func,'lens':lensav_func}
myMISN3 = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.33,z_lens=.53, bands=['F110W','F125W'],
              zp=[26.8,26.2], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[7,3.5],
              objectName='My Type Ia SN',telescopename='HST',av_dists=param_funcs)
print(myMISN3)
myMISN3.plot_object()
