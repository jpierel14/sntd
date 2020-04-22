"""
=======================
Simulating Microlensing
=======================

Simulate a microlensing microcaustic.
"""

import sntd
import numpy as np
       
myML=sntd.realizeMicro(nray=50,kappas=1,kappac=.3,gamma=.4)
time,dmag=sntd.microcaustic_field_to_curve(field=myML,time=np.arange(0,100,1),zl=.5,zs=1.33,plot=True)

###############################################################
# **Including Microlensing in Simulations**
# Now we can take the simulated microcaustic 
# and use it to include microlensing in a 
# multiply-imaged supernova simulation.



myMISN2 = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.2,z_lens=.5, bands=['F110W','F160W'],
                   zp=[26.8,26.2], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[7,3.5],
       objectName='My Type Ia SN',telescopename='HST', microlensing_type='AchromaticMicrolensing',microlensing_params=myML)
myMISN2.plot_object(showMicro=True)