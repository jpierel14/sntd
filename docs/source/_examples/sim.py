"""
===========
Simulations
===========

Create a simulated multiply-imaged supernova that we can then fit,
with no microlensing included in the simulation. Note that your final
printed information will be different, as this is a randomly generated
supernova.
"""
		
###############################################################
# **No Microlensing**
   
import sntd

myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.33,z_lens=.53, bands=['F110W','F125W'],
             zp=[26.8,26.2], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[7,3.5],
 objectName='My Type Ia SN',telescopename='HST')
print(myMISN)
myMISN.plot_object()



