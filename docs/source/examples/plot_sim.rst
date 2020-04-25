.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_examples_plot_sim.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_examples_plot_sim.py:


=====================
Simulating Supernovae
=====================

Simulating a multiply-imaged supernova.

Create a simulated multiply-imaged supernova that we can then fit,
with no microlensing included in the simulation. Note that your final
printed information will be different, as this is a randomly generated
supernova. The function being used in these examples is 
`~sntd.createMultiplyImagedSN` .  

**No Microlensing**


.. code-block:: default

   
    import sntd
    import numpy as np

    myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.4,z_lens=.53, bands=['F110W','F160W'],
                 zp=[26.8,26.2], cadence=5., epochs=35.,time_delays=[20., 70.], magnifications=[10,5],
     objectName='My Type Ia SN',telescopename='HST',av_host=False)
    print(myMISN)
    myMISN.plot_object()





.. image:: /examples/images/sphx_glr_plot_sim_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Telescope: HST
    Object: My Type Ia SN
    Number of bands: 2

    ------------------
    Image: image_1:
    Bands: ['F110W', 'F160W']
    Date Range: 0.00000->138.97059
    Number of points: 56

    Metadata:
       z:1.4
       t0:20.0
       x0:4.717933761242799e-06
       x1:2.0391629732663894
       c:-0.05282284476175894
       sourcez:1.4
       hostebv:0
       lensebv:0
       lensz:0.53
       mu:10
       td:20.0
    ------------------
    Image: image_2:
    Bands: ['F110W', 'F160W']
    Date Range: 25.73529->175.00000
    Number of points: 57

    Metadata:
       z:1.4
       t0:70.0
       x0:2.3589668806213994e-06
       x1:2.0391629732663894
       c:-0.05282284476175894
       sourcez:1.4
       hostebv:0
       lensebv:0
       lensz:0.53
       mu:5
       td:70.0
    ------------------

    <Figure size 1000x1000 with 2 Axes>



Specify the distributions you want to use for any model
parameter by providing a function that returns the parameter
in any way you want. 


.. code-block:: default


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





.. image:: /examples/images/sphx_glr_plot_sim_002.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Telescope: HST
    Object: My Type Ia SN
    Number of bands: 2

    ------------------
    Image: image_1:
    Bands: ['F125W', 'F110W']
    Date Range: 0.00000->123.52941
    Number of points: 50

    Metadata:
       z:1.33
       t0:10.0
       x0:6.494176446119458e-06
       x1:1.131118738018449
       c:-0.05397004195498919
       sourcez:1.33
       hostebv:0.0967741935483871
       lensebv:0
       lensz:0.53
       mu:7
       td:10.0
    ------------------
    Image: image_2:
    Bands: ['F125W', 'F110W']
    Date Range: 25.73529->175.00000
    Number of points: 56

    Metadata:
       z:1.33
       t0:70.0
       x0:3.247088223059729e-06
       x1:1.131118738018449
       c:-0.05397004195498919
       sourcez:1.33
       hostebv:0.0967741935483871
       lensebv:0
       lensz:0.53
       mu:3.5
       td:70.0
    ------------------

    <Figure size 1000x1000 with 2 Axes>



Specify the distributions you want to use for dust
parameters by providing a function that returns the parameter
in any way you want. 


.. code-block:: default


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



.. image:: /examples/images/sphx_glr_plot_sim_003.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Telescope: HST
    Object: My Type Ia SN
    Number of bands: 2

    ------------------
    Image: image_1:
    Bands: ['F125W', 'F110W']
    Date Range: 0.00000->123.52941
    Number of points: 50

    Metadata:
       z:1.33
       t0:10.0
       x0:7.620855568458378e-06
       x1:1.2103272697854583
       c:-0.00832174205257879
       sourcez:1.33
       hostebv:0.18363487218657767
       lensebv:0.20937423654696036
       lensz:0.53
       mu:7
       td:10.0
    ------------------
    Image: image_2:
    Bands: ['F125W', 'F110W']
    Date Range: 25.73529->175.00000
    Number of points: 58

    Metadata:
       z:1.33
       t0:70.0
       x0:3.810427784229189e-06
       x1:1.2103272697854583
       c:-0.00832174205257879
       sourcez:1.33
       hostebv:0.18363487218657767
       lensebv:0.20937423654696036
       lensz:0.53
       mu:3.5
       td:70.0
    ------------------

    <Figure size 1000x1000 with 2 Axes>




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  8.761 seconds)


.. _sphx_glr_download_examples_plot_sim.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_sim.py <plot_sim.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_sim.ipynb <plot_sim.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
