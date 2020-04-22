"""
=====================
Measuring Time Delays
=====================

There are 3 methods built into SNTD to measure time delays 
(parallel, series, color). They are accessed by the same 
function:`~sntd.fit_data`. Here ``myMISN`` could be generated
in the simulation example of the documentation.
"""

###############################################################
# **Parallel:**

import sntd

ex_1,ex_2=sntd.load_example_data()
myMISN=sntd.table_factory([ex_1,ex_2],telescopename='HST',object_name='example_SN')

fitCurves=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
                params=['x0','t0','x1','c'],constants={'z':1.2},refImage='image_1',
                bounds={'t0':(-20,20),'x1':(-3,3),'c':(-1,1),'mu':(.5,2)},fitOrder=['image_2','image_1'],
                method='parallel',microlensing=None,modelcov=False,npoints=500,maxiter=None)
print(fitCurves.parallel.time_delays)
print(fitCurves.parallel.time_delay_errors)
print(fitCurves.parallel.magnifications)
print(fitCurves.parallel.magnification_errors)
fitCurves.plot_object(showFit=True,method='parallel')
fitCurves.plot_fit(method='parallel')
fitCurves.plot_fit(method='parallel')


#####################################################################################################################################
# Note that the bounds for the 't0' parameter are not absolute, the actual peak time will be estimated (unless t0_guess is defined)
# and the defined bounds will be added to this value. Similarly for amplitude, where bounds are multiplicative
#
# Other methods are called in a similar fashion, with a couple of extra arguments:
# 
# **Series:**

    
fitCurves=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
        params=['x1','c'],constants={'z':zs},refImage='image_1',
        bounds={'td':(-20,20),'mu':(.5,2),'x1':(-3,3),'c':(-1,1)},
        method='series',npoints=500)


print(fitCurves.series.time_delays)
print(fitCurves.series.time_delay_errors)
print(fitCurves.series.magnifications)
print(fitCurves.series.magnification_errors)
fitCurves.plot_object(showFit=True,method='series')
fitCurves.plot_fit(method='series')

############################################################### 
# **Color:**


    
fitCurves=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
                    params=['c'],constants={'z':zs,'x1':fitCurves.images['image_1'].fits.model.get('x1')},refImage='image_1',
                    bounds={'td':(-20,20),'mu':(.5,2),'x1':(-3,3),'c':(-1,1)},
                    method='color',microlensing=None,modelcov=False,npoints=500,maxiter=None)

print(fitCurves.color.time_delays)
print(fitCurves.color.time_delay_errors)
fitCurves.plot_object(showFit=True,method='color')
fitCurves.plot_fit(method='color')

#################################################################################################################################################################
# You can include your fit from the parallel method as a prior on light curve and time delay parameters in the series/color methods with the "fit_prior" command:



fitCurves_parallel=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
                	params=['x0','t0','x1','c'],constants={'z':1.2},refImage='image_1',
                	bounds={'t0':(-20,20),'x1':(-3,3),'c':(-1,1),'mu':(.5,2)},fitOrder=['image_2','image_1'],
               	    method='parallel',microlensing=None,modelcov=False,npoints=500,maxiter=None)
fitCurves_color=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
                	params=['c'],constants={'z':zs,'x1':fitCurves.images['image_1'].fits.model.get('x1')},refImage='image_1',
                	bounds={'td':(-20,20),'mu':(.5,2),'x1':(-3,3),'c':(-1,1)},fit_prior=fitCurves_parallel,
                	method='color',microlensing=None,modelcov=False,npoints=500,maxiter=None)

###################################################################################################################################
# **Fitting Using Extra Propagation Effects**
# 
# You might also want to include other propagation effects in your fitting model, and fit relevant parameters. This can be done by
# simply adding effects to an SNCosmo model, in the same way as if you were fitting a single SN with SNCosmo. First we can add some
# extreme dust in the source and lens frames (your final simulations may look slightly different as **c** is chosen randomly):



myMISN2 = sntd.createMultiplyImagedSN(sourcename='salt2', snType='Ia', redshift=1.45,z_lens=.53, bands=['F110W','F160W'],
              zp=[26.9,26.2], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[10,5],
              objectName='My Type Ia SN',telescopename='HST',av_lens=1.5,
              av_host=1)
print(myMISN2.images['image_1'].simMeta['lensebv'],
     myMISN2.images['image_1'].simMeta['hostebv'], 
     myMISN2.images['image_1'].simMeta['c'])

################################################################################
# Okay, now we can fit the MISN first without taking these effects into account:



fitCurves=sntd.fit_data(myMISN2,snType='Ia', models='salt2',bands=['F110W','F160W'],
                                                     params=['x0','x1','t0','c'],
                                                     constants={'z':1.45},
                                                     bounds={'t0':(-15,15),'x1':(-2,2),'c':(-1,1)},
                                                     showPlots=True)

########################################################################################################################
# We can see that the fitter has done reasonably well, and the time delay is still accurate (True delay is 60 days). 
# However, one issue is that the measured value for **c** (0.805) is vastly different than the actual value (0.098) 
# as it attempts to compensate for extinction without a propagation effect. Now let's add in the propagation effects:


dust = sncosmo.CCM89Dust()
salt2_model=sncosmo.Model('salt2',effects=[dust,dust],effect_names=['lens','host'],effect_frames=['free','rest'])
fitCurves=sntd.fit_data(myMISN2,snType='Ia', models=salt2_model,bands=['F110W','F160W'],
                    params=['x0','x1','t0','c','lensebv','hostebv'],
                    constants={'z':1.45,'lensr_v':3.1,'lensz':0.53,'hostr_v':3.1},
                    bounds={'t0':(-15,15),'x1':(-2,2),'c':(-1,1),'lensebv':(0,1.),'hostebv':(0,1.)},
                    showPlots=True)

########################################################################################################################
# Now the measured value for **c** (0.057) is much closer to reality, and the measured times of peak are somewhat
# more accurate. 
