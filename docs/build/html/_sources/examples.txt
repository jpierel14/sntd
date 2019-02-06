**********************
Simulating with SNTD
**********************

No Microlensing
================

Create a simulated multiply-imaged supernova that we can then fit,
with no microlensing included in the simulation. Note that your final
printed light curve will be different, as this is a randomly generated
supernova.

		
.. code-block:: python     
	
     import sntd
     import matplotlib.pyplot as plt
     
     myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.33,z_lens=.53, bands=['F110W','F125W'],
                                     zp=[26.8,26.2], cadence=5., epochs=35.,
                                     time_delays=[10., 70.], magnifications=[7,3.5], objectName='My Type Ia SN',telescopename='HST')
     print(myMISN.images['image1'].table)
     myMISN.plot_object()
     plt.show()

Out::
  
  time               band        flux             fluxerr        zp  zpsys  image 
  ------------------ ----- ------------------ ------------------ ---- ----- -------
  0.0                F110W  80.83704095305713 1.2196369693362414 26.8    AB image_1
  5.147058823529412  F110W  89.91280021212975  1.284158612608907 26.8    AB image_1
  10.294117647058824 F110W  96.27127396739887 1.2466015593879611 26.8    AB image_1
  15.441176470588236 F110W  90.70821457675255   1.22366765470951 26.8    AB image_1
  20.58823529411765  F110W  86.53240018017523 1.2073336440515443 26.8    AB image_1
  25.73529411764706  F110W  76.07501959744152 1.2094036308608453 26.8    AB image_1
  30.88235294117647  F110W  63.17842848471685 1.1769008597903798 26.8    AB image_1
  36.029411764705884 F110W  54.15614809883671 1.1548804427431165 26.8    AB image_1
  41.1764705882353   F110W  44.63934115659113 1.1351349145266738 26.8    AB image_1
  46.32352941176471  F110W 37.517811149554156 1.1509488564354409 26.8    AB image_1
  51.47058823529412  F110W 31.368919986394147 1.1348214534140784 26.8    AB image_1
  56.617647058823536 F110W 26.365621628734704 1.1473493554620837 26.8    AB image_1
  0.0                F125W 44.652886569466006 1.1844546865198073 26.2    AB image_1
  5.147058823529412  F125W  52.12506595181393 1.1882051108553613 26.2    AB image_1
  10.294117647058824 F125W    56.383959775956 1.2055312873637205 26.2    AB image_1
  15.441176470588236 F125W 52.282478317158144 1.2019364467744396 26.2    AB image_1
  20.58823529411765  F125W  51.44217018263837 1.1737524525502405 26.2    AB image_1
  25.73529411764706  F125W  45.41787029296782 1.1346437770894173 26.2    AB image_1
  30.88235294117647  F125W   41.7399546859917  1.138562946339325 26.2    AB image_1
  36.029411764705884 F125W 36.249141725372375 1.1633746045379156 26.2    AB image_1
  41.1764705882353   F125W 30.134015959569542  1.080202890453433 26.2    AB image_1
  46.32352941176471  F125W 26.732438252265634 1.1574450238913505 26.2    AB image_1
  51.47058823529412  F125W 27.135405865133855  1.071456097317496 26.2    AB image_1
  56.617647058823536 F125W 21.536294532971315  1.081999880063472 26.2    AB image_1

Out:

.. image:: examples/noML.png
    :width: 600px
    :align: center
    :height: 600px
    :alt: alternate text
	  
Simulating Microlensing
=======================
Simulate a microlensing microcaustic, and use it to include a microlensing effect in
the simulated supernova.

.. code-block:: python

		import numpy as np
		
		myML=sntd.realizeMicro(nray=50,kappas=1,kappac=.3,gamma=.4)
		time,dmag=sntd.microcaustic_field_to_curve(field=fitData,time=np.arange(0,100,1),zl=.5,zs=1,plot=True)
		plt.show()

Out:

.. image:: examples/micro.png
    :width: 600px
    :align: center
    :height: 600px
    :alt: alternate text

Including Microlensing in Simulations
=====================================
Now we can take the simulated microcaustic and use it to include microlensing in a multiply-imaged supernova simulation.

.. code-block:: python

		myMISN2 = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.33,z_lens=.53, bands=['F110W','F125W'],
                                     zp=[26.8,26.2], cadence=5., epochs=35.,
                                     time_delays=[10., 70.], magnifications=[7,3.5], objectName='My Type Ia SN',telescopename='HST',
				     microlensing_type='AchromaticMicrolensing',microlensing_params=myML)
		myMISN2.plot_object(showMicro=True,showModel=True)

Out:

.. image:: examples/withML.png
    :width: 600px
    :align: center
    :height: 600px
    :alt: alternate text


*******************************
Measuring Time Delays with SNTD
*******************************

Fitting a Multiply-Imaged Supernova
===================================
There are 3 methods built into SNTD to measure time delays (separate, combined, color). They are accessed by the same function:

.. code-block:: python

		fitCurves=sntd.fit_data(myMISN2,snType='Ia', models='salt2-extended',bands=['F110W','F125W'],
                            params=['x0','x1','t0','c'],constants={'z':1.33},
                            bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},method='separate',microlensing=None)
		fitCurves.plot_object(showFit=True,method='separate')
		plt.show()

Out:

.. image:: examples/withML.png
    :width: 600px
    :align: center
    :height: 600px
    :alt: alternate text

Estimating Uncertainty Due to Microlensing
=========================================
