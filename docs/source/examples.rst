*******************
Using Your Own Data
*******************
In order to fit your own data, you must turn your light curve into an astropy table. There is an example multiply-imaged
SN example provided for reference. In this example, we have a doubly-imaged SN with image files (in the sntd/data/examples folder)
'example_image_1.dat' and 'example_image_2.dat'. The only optional column in these files is "image", which sets the name of the key
used to reference this SN image. If you do not provide flux/fluxerr but instead magnitude/magerr SNTD will attemp to translate to
flux/fluxerr, but it's best to simply provide flux from the beginning to avoid conversion errors. First we can read in these tables:

.. code-block:: python
	
	ex_1,ex_2=sntd.load_example_data()
	print(ex_1)

Out:: 

	       time         band        flux        ...  zp  zpsys  image 
	------------------ ----- ------------------ ... ---- ----- -------
	               0.0 F125W  64.59429430606906 ... 26.8    AB image_1
	2.0224719101123596 F125W    62.408324396966 ... 26.8    AB image_1
	 4.044943820224719 F125W  68.10359798573809 ... 26.8    AB image_1
	 6.067415730337078 F125W  71.76160753594853 ... 26.8    AB image_1
	 8.089887640449438 F125W  73.43467553050705 ... 26.8    AB image_1
	10.112359550561798 F125W  74.34296720689296 ... 26.8    AB image_1
	12.134831460674157 F125W  71.73347707161632 ... 26.8    AB image_1
	14.157303370786517 F125W  72.93187923529568 ... 26.8    AB image_1
	16.179775280898877 F125W  70.64111678688164 ... 26.8    AB image_1
	18.202247191011235 F125W  69.31085357488871 ... 26.8    AB image_1
	               ...   ...                ... ...  ...   ...     ...
	38.426966292134836 F160W 19.950527074094737 ... 26.2    AB image_1
	40.449438202247194 F160W 20.963076283234553 ... 26.2    AB image_1
	 42.47191011235955 F160W 21.402880246191344 ... 26.2    AB image_1
	 44.49438202247191 F160W  18.28098879531828 ... 26.2    AB image_1
	 46.51685393258427 F160W 18.947732390210522 ... 26.2    AB image_1
	 48.53932584269663 F160W 15.987591900959364 ... 26.2    AB image_1
	 50.56179775280899 F160W 20.011941798193966 ... 26.2    AB image_1
	 52.58426966292135 F160W 15.516064719260328 ... 26.2    AB image_1
	 54.60674157303371 F160W   17.1543325162061 ... 26.2    AB image_1
	 56.62921348314607 F160W  18.25136177909449 ... 26.2    AB image_1
	58.651685393258425 F160W 17.198071229182016 ... 26.2    AB image_1
	Length = 60 rows

Now, to turn these two data tables into an SNTD curveDict object that will be fit, we use the table_factory function:

.. code-block:: python

	new_MISN=sntd.table_factory([ex_1,ex_2],telescopename='HST',object_name='example_SN')
	print(new_MISN)

Out::

	Telescope: HST
	Object: example_SN
	Number of bands: 2

	------------------
	Image: image_1:
	Bands: set(['F160W', 'F125W'])
	Date Range: 0.00000->58.65169
	Number of points: 60
	------------------
	Image: image_2:
	Bands: set(['F160W', 'F125W'])
	Date Range: 40.44944->119.32584
	Number of points: 80
	------------------

And finally let's fit this SN, which is a Type Ia, with the SALT2 model (your exact time delay may
be slightly different after fitting the example data). For reference, the true delay here is 60 days.


.. code-block:: python

	fitCurves=sntd.fit_data(new_MISN,snType='Ia', models='salt2',bands=['F125W','F160W'],
                        params=['x0','x1','t0','c'],constants={'z':1.33},
                        bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)})
	print(fitCurves.parallel.time_delays)
	fitCurves.plot_object(showFit=True)
	plt.show()


Out::

	{'image_1': 0, 'image_2': 60.2649320870058}

.. image:: examples/example_fit.png
    :width: 600px
    :align: center
    :height: 600px
    :alt: alternate text


****************************************
Batch Processing Time Delay Measurements
****************************************

Parallel processing and batch processing is built into SNTD in order to fit a large number of (likely simulated) MISN. To access this feature,
simply provide a list of MISN instead of a single sntd curveDict object, specifying whether you want to use multiprocessing (split the list across multiple cores)
or batch processing (splitting the list into multiple jobs with sbatch). If you specify batch mode, you need to provide
the partition and number of jobs you want to implement. 

.. code-block:: python

  myMISN1 = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.33,z_lens=.53, bands=['F110W','F125W'],
                   zp=[26.8,26.2], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[7,3.5],
       objectName='My Type Ia SN',telescopename='HST')
  myMISN2 = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.33,z_lens=.53, bands=['F110W','F125W'],
                   zp=[26.8,26.2], cadence=5., epochs=35.,time_delays=[10., 50.], magnifications=[7,3.5],
       objectName='My Type Ia SN',telescopename='HST')
  curve_list=[myMISN1,myMISN2]
  fitCurves=sntd.fit_data(curve_list,snType='Ia', models='salt2-extended',bands=['F110W','F125W'],
                    params=['x0','t0','x1','c'],constants={'z':1.3},refImage='image_1',
                    bounds={'t0':(-20,20),'x1':(-3,3),'c':(-1,1)},fitOrder=['image_2','image_1'],
                    method='parallel',npoints=1000,par_or_batch='batch', batch_partition='myPartition',nbatch_jobs=2)

  for curve in fitCurves:
    print(curve.parallel.time_delays)
  
  fitCurves=sntd.fit_data(curve_list,snType='Ia', models='salt2-extended',bands=['F110W','F125W'],
                    params=['x0','t0','x1','c'],constants={'z':1.3},refImage='image_1',
                    bounds={'t0':(-20,20),'x1':(-3,3),'c':(-1,1)},fitOrder=['image_2','image_1'],
                    method='parallel',npoints=1000,par_or_batch='parallel')
  for curve in fitCurves:
    print(curve.parallel.time_delays)

Out::

  Submitted batch job 5784720
  {'image_1': 0, 'image_2': 60.3528844834}
  {'image_1': 0, 'image_2': 40.34982372733}
  Fitting MISN number 1...
  Fitting MISN number 2...
  {'image_1': 0, 'image_2': 60.32583528844834}
  {'image_1': 0, 'image_2': 40.22834982372733}


If you would like to run multiple methods in a row in batch mode, the recommended way is by providing a list of the methods to the fit_data function. You 
can have it use the parallel fit as a prior on the subsequent fits by setting fit_prior to True instead of giving it a curveDict object.


.. code-block:: python

  
  fitCurves_batch=sntd.fit_data(curve_list,snType='Ia', models='salt2-extended',bands=['F110W','F125W'],
                    params=['x0','t0','x1','c'],constants={'z':1.3},refImage='image_1',fit_prior=True,
                    bounds={'t0':(-20,20),'x1':(-3,3),'c':(-1,1)},fitOrder=['image_2','image_1'],
                    method=['parallel','series','color'],npoints=1000,par_or_batch='batch', batch_partition='myPartition',nbatch_jobs=2)
