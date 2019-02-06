**********************
Core-Collapse Examples
**********************

Getting Started
===============

Read in a lightcurve and the filename of an SED:

		
.. code-block:: python     
	
     from __future__ import print_function
     import snsedextend
     sedFile=snsedextend.example_sed
     myLC=snsedextend.load_example_lc()
     print(myLC)

Out::
  
    name  band     time     mag   magerr
    ------ ---- ----------- ------ ------
    2006aj    U 53788.16733 18.039  0.078
    2006aj    U 53790.12301 17.995  0.047
    2006aj    U 53790.13652 17.982  0.057
    2006aj    U 53796.11875 18.599  0.112
    2006aj    U 53797.11978 18.679  0.054
    2006aj    U 53798.16206 18.738  0.058
    2006aj    B 53788.11683 18.438   0.05
    ...  ...         ...    ...    ...
    2006aj    H    53793.12 16.752  0.047
    2006aj    H   53798.145 16.635  0.037
    2006aj    K    53788.16 17.183   0.12
    2006aj    K   53790.133 16.714  0.094
    2006aj    K   53792.123 16.406  0.145
    2006aj    K    53793.12 16.393   0.08
    2006aj    K   53798.145 16.769   0.08
    Length = 105 rows

Generate a Color Table
======================
Produce a color table from the example data with some assumptions. You can set any parameter that you would like that is used for SNCosmo fitting.

.. code-block:: python

		colorTable=snsedextend.curveToColor(myLC,colors=['U-B', 'r-J', 'r-H', 'r-K'],
		                 snType='Ic', zpsys='vega', bounds={'hostebv': (-1, 1),
				 't0': (53787.94, 53797.94)},constants={'mwr_v': 3.1,
				 'mwebv': '0.1267', 'z': '0.033529863', 'hostr_v': 3.1},
				 dust='CCM89Dust', effect_frames=['rest', 'obs'], effect_names=['host', 'mw'])

Out:

.. literalinclude:: examples/ex.txt
	:language: html
	:lines: 1-3

Now print the result:

.. code-block:: python
		
	print(colorTable)

Out:

.. literalinclude:: examples/ex.txt
	:language: html
	:lines: 4-21

Color Curve Fitting
===================
Now we can fit this color table and get a best model by minimizing BIC.
This function returns a python dictionary with colors as keys and an astropy Table object
with time and color vectors as values::

  curveDict=snsedextend.fitColorCurve(colorTable)

Out:

.. literalinclude:: examples/ex.txt
	:language: html
	:lines: 23-27
		
SED Extrapolation
=================
Now you can provide an SED to be extrapolated, and let it do the work (This is a type Ic). This will return an
sncosmo.Source object and simultanously save the new SED to a file defined by newFileLoc (default current directory).::

  newSED=snsedextend.extendCC(colorTable,curveDict,sedlist=[sedFile],zpsys='vega',showplots=False,verbose=True)

Plotting from Timeseries
========================
You can directly plot an SED file.::
  
  snsedextend.plotSED('SEDS/typeIc/SDSS-004012.SED',day=1,showPlot=True,MINWAVE=4000,MAXWAVE=20000,saveFig=False)

Out:

.. image:: examples/example_plot.png
    :width: 600px
    :align: center
    :height: 400px
    :alt: alternate text
