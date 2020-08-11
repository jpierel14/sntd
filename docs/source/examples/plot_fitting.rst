.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_examples_plot_fitting.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_examples_plot_fitting.py:


===================
Measure Time Delays
===================

A series of examples demonstrating various fitting options 
with SNTD.

There are 3 methods built into SNTD to measure time delays
(parallel, series, color). They are accessed by the same 
function: :py:func:`~sntd.fitting.fit_data` . 
Here ``myMISN`` was generated in the :ref:`examples/plot_sim:Simulating Supernovae` part 
of the documentation, using the :py:func:`~sntd.simulation.createMultiplyImagedSN` 
function. The true delay for all of these fits is 50 days.
You can batch process (with sbatch or multiprocessing) using any or all of these methods as well 
(see :ref:`examples:Batch Processing Time Delay Measurements`)

 -----------------------------------
`Run this notebook with Google Colab <https://colab.research.google.com/github/jpierel14/sntd/blob/master/notebooks/docs_fitting.ipynb>`_.
 -----------------------------------

**Parallel:**


.. code-block:: default

    import sntd

    myMISN=sntd.load_example_misn()

    fitCurves=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
                    params=['x0','t0','x1','c'],constants={'z':1.4},refImage='image_1',cut_time=[-30,40],
                    bounds={'t0':(-20,20),'x1':(-3,3),'c':(-.5,.5),'mu':(.5,2)},fitOrder=['image_2','image_1'],
                    method='parallel',microlensing=None,modelcov=False,npoints=100)
    print(fitCurves.parallel.time_delays)
    print(fitCurves.parallel.time_delay_errors)
    print(fitCurves.parallel.magnifications)
    print(fitCurves.parallel.magnification_errors)
    fitCurves.plot_object(showFit=True,method='parallel')
    fitCurves.plot_fit(method='parallel',par_image='image_1')
    fitCurves.plot_fit(method='parallel',par_image='image_2')




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /examples/images/sphx_glr_plot_fitting_001.png
            :class: sphx-glr-multi-img

    *

      .. image:: /examples/images/sphx_glr_plot_fitting_002.png
            :class: sphx-glr-multi-img

    *

      .. image:: /examples/images/sphx_glr_plot_fitting_003.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    {'image_1': 0, 'image_2': 49.840019971680576}
    {'image_1': array([0, 0]), 'image_2': array([-0.15335578,  0.15516413])}
    {'image_1': 1, 'image_2': 0.5026350011268993}
    {'image_1': array([0, 0]), 'image_2': array([-0.00463769,  0.00505748])}

    <Figure size 970x970 with 16 Axes>



Note that the bounds for the 't0' parameter are not absolute, the actual peak time will be estimated (unless t0_guess is defined)
and the defined bounds will be added to this value. Similarly for amplitude, where bounds are multiplicative

Other methods are called in a similar fashion, with a couple of extra arguments:

**Series:**


.. code-block:: default



    fitCurves=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
            params=['x0','t0','x1','c'],constants={'z':1.4},refImage='image_1',cut_time=[-30,40],
            bounds={'t0':(-20,20),'td':(-20,20),'mu':(.5,2),'x1':(-3,3),'c':(-.5,.5)},
            method='series',npoints=100)
        

    print(fitCurves.series.time_delays)
    print(fitCurves.series.time_delay_errors)
    print(fitCurves.series.magnifications)
    print(fitCurves.series.magnification_errors)
    fitCurves.plot_object(showFit=True,method='series')
    fitCurves.plot_fit(method='series')




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /examples/images/sphx_glr_plot_fitting_004.png
            :class: sphx-glr-multi-img

    *

      .. image:: /examples/images/sphx_glr_plot_fitting_005.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [ 2.13811422e-05  1.99464548e+01 -1.65168083e-01  6.53175153e-02
      4.98454606e+01  5.04624203e-01]
    {'image_1': 0, 'image_2': 49.841793423312915}
    {'image_1': array([0, 0]), 'image_2': array([-0.06908658,  0.07011961])}
    {'image_1': 1, 'image_2': 0.5046481823741065}
    {'image_1': array([0, 0]), 'image_2': array([-0.0016663 ,  0.00161565])}

    <Figure size 1390x1390 with 36 Axes>



**Color:**
By default, this will attempt to fit every combination of colors possible from
the bands present in the data. You can define specific colors using the "fit_colors"
argument.


.. code-block:: default



    
    fitCurves=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
                        params=['t0','c'],constants={'z':1.4,'x1':fitCurves.images['image_1'].fits.model.get('x1')},refImage='image_1',
                        color_param_ignore=['x1'],bounds={'t0':(-20,20),'td':(-20,20),'mu':(.5,2),'c':(-.5,.5)},cut_time=[-30,40],
                        method='color',microlensing=None,modelcov=False,npoints=200,maxiter=None,minsnr=3)

    print(fitCurves.color.time_delays)
    print(fitCurves.color.time_delay_errors)
    fitCurves.plot_object(showFit=True,method='color')
    fitCurves.plot_fit(method='color')




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /examples/images/sphx_glr_plot_fitting_006.png
            :class: sphx-glr-multi-img

    *

      .. image:: /examples/images/sphx_glr_plot_fitting_007.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    {'image_1': 0, 'image_2': 51.52545264404653}
    {'image_1': array([0, 0]), 'image_2': array([-1.23242995,  1.34490708])}

    <Figure size 760x760 with 9 Axes>



You can include your fit from the parallel method as a prior on light curve and time delay parameters in the series/color methods with the "fit_prior" command:


.. code-block:: default




    fitCurves_parallel=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
                    	params=['x0','t0','x1','c'],constants={'z':1.4},refImage='image_1',
                    	bounds={'t0':(-20,20),'x1':(-3,3),'c':(-.5,.5),'mu':(.5,2)},fitOrder=['image_2','image_1'],cut_time=[-30,40],
                   	    method='parallel',microlensing=None,modelcov=False,npoints=100,maxiter=None)
    fitCurves_color=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],cut_time=[-50,30],
                    	params=['t0','c'],constants={'z':1.4,'x1':fitCurves.images['image_1'].fits.model.get('x1')},refImage='image_1',
                    	bounds={'t0':(-20,20),'td':(-20,20),'mu':(.5,2),'c':(-.5,.5)},fit_prior=fitCurves_parallel,
                    	method='color',microlensing=None,modelcov=False,npoints=200,maxiter=None,minsnr=3)

    print(fitCurves_parallel.parallel.time_delays)
    print(fitCurves_parallel.parallel.time_delay_errors)
    print(fitCurves_color.color.time_delays)
    print(fitCurves_color.color.time_delay_errors)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    {'image_1': 0, 'image_2': 49.83642937969682}
    {'image_1': array([0, 0]), 'image_2': array([-0.16216226,  0.14232851])}
    {'image_1': 0, 'image_2': 49.82679825541822}
    {'image_1': array([0, 0]), 'image_2': array([-0.30500095,  0.30527895])}




**Fitting Using Extra Propagation Effects**

You might also want to include other propagation effects in your fitting model, and fit relevant parameters. This can be done by
simply adding effects to an SNCosmo model, in the same way as if you were fitting a single SN with SNCosmo. First we can add some
extreme dust in the source and lens frames (your final simulations may look slightly different as **c** is chosen randomly):


.. code-block:: default




    myMISN2 = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.4,z_lens=.53, bands=['F110W','F160W'],
                  zp=[26.9,26.2], cadence=8., epochs=30.,time_delays=[20., 70.], magnifications=[20,10],
                  objectName='My Type Ia SN',telescopename='HST',av_lens=1.5,
                  av_host=1)
    print('lensebv:',myMISN2.images['image_1'].simMeta['lensebv'],
         'hostebv:',myMISN2.images['image_1'].simMeta['hostebv'], 
         'c:',myMISN2.images['image_1'].simMeta['c'])





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    lensebv: 0.48387096774193544 hostebv: 0.3225806451612903 c: 0.0036396427131899175




Okay, now we can fit the MISN first without taking these effects into account:


.. code-block:: default




    fitCurves_dust=sntd.fit_data(myMISN2,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
                                                         params=['x0','x1','t0','c'],npoints=200,
                                                         constants={'z':1.4},minsnr=1,cut_time=[-30,40],
                                                         bounds={'t0':(-15,15),'x1':(-3,3),'c':(-.3,.3)})
    print(fitCurves_dust.parallel.time_delays)
    print(fitCurves_dust.parallel.time_delay_errors)
    print('c:',fitCurves_dust.images['image_1'].fits.model.get('c'))
    fitCurves_dust.plot_object(showFit=True)



.. image:: /examples/images/sphx_glr_plot_fitting_008.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    {'image_1': 0, 'image_2': 50.38279703244322}
    {'image_1': array([0, 0]), 'image_2': array([-0.6135019 ,  0.68283929])}
    c: 0.7402092975258244

    <Figure size 1000x1000 with 2 Axes>



We can see that the fitter has done reasonably well, and the time delay is still accurate (True delay is 50 days). 
However, one issue is that the measured value for **c** is vastly different than the actual value 
as it attempts to compensate for extinction without a propagation effect. Now let's add in the propagation effects:


.. code-block:: default


    import sncosmo
    dust = sncosmo.CCM89Dust()
    salt2_model=sncosmo.Model('salt2-extended',effects=[dust,dust],effect_names=['lens','host'],effect_frames=['free','rest'])
    fitCurves_dust=sntd.fit_data(myMISN2,snType='Ia', models=salt2_model,bands=['F110W','F160W'],npoints=200,
                        params=['x0','x1','t0','c','lensebv','hostebv'],minsnr=1,cut_time=[-30,40],
                        constants={'z':1.4,'lensr_v':3.1,'lensz':0.53,'hostr_v':3.1},
                        bounds={'t0':(-15,15),'x1':(-3,3),'c':(-.1,.1),'lensebv':(.2,1.),'hostebv':(.2,1.)})

    print(fitCurves_dust.parallel.time_delays)
    print(fitCurves_dust.parallel.time_delay_errors)
    print('c:',fitCurves_dust.images['image_1'].fits.model.get('c'),
          'lensebv:',fitCurves_dust.images['image_1'].fits.model.get('lensebv'),
          'hostebv:',fitCurves_dust.images['image_1'].fits.model.get('hostebv'))
    fitCurves_dust.plot_object(showFit=True)



.. image:: /examples/images/sphx_glr_plot_fitting_009.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    {'image_1': 0, 'image_2': 50.40527972431711}
    {'image_1': array([0, 0]), 'image_2': array([-0.66916388,  0.71217856])}
    c: 0.28127107977898574 lensebv: 0.3963084663702129 hostebv: 0.1491811267077542

    <Figure size 1000x1000 with 2 Axes>



Now the measured value for **c** is much closer to reality, and the measured times of peak are somewhat
more accurate. 


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 4 minutes  38.032 seconds)


.. _sphx_glr_download_examples_plot_fitting.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_fitting.py <plot_fitting.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_fitting.ipynb <plot_fitting.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
