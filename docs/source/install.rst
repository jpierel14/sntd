************
Installation
************

SNTD works on Python 3.5+ and requires the
following Python packages:

- `numpy <http://www.numpy.org/>`_
- `scipy <http://www.scipy.org/>`_
- `astropy <http://www.astropy.org>`_
- `SNCosmo <http://sncosmo.readthedocs.io>`_
- `cython <https://cython.org/>`_
- `sklearn <https://scikit-learn.org/stable/tutorial/basic/tutorial.html>`_
- `nestle <https://github.com/kbarbary/nestle>`_


Install using pip
=================

Using pip::

    pip install sntd

If you plan to use SNTD to simulate the effects of microlensing,
which uses the `[Wambsganss 1999] <https://www.sciencedirect.com/science/article/pii/S0377042799001648>`_
microlens code, then you will need
a fortran compiler. You can install gfortran `here <https://gcc.gnu.org/wiki/GFortranBinaries>`_.
   

Common Installation Issues
==========================

1. You will need a C compiler (e.g. ``gcc`` or ``clang``) to be
   installed for the installation to succeed due to SNCosmo.

2. Sometimes there is an issue with the install related to numpy,
   particularly if you are installing in a fresh build of python
   (e.g. a new virtual environment). If this happens, try
   installing numpy using pip, and then installing SNTD using pip.

3. On MacOS and a fresh python build, you may have an issue with
   importing sntd because of matplotlib and an error saying python
   was not installed as a framework. A fix for that issue is on
   `stack overflow <https://stackoverflow.com/questions/21784641/installation-issue-with-matplotlib-python>`_.

4. Sometimes you can get an error when using the microlensing component of SNTD involving the 
   `microlens` executable. This involves differences between the machine
   the microlens executable was created on and yours. Just go to the 
   `microlens github repo <https://github.com/psaha/microlens>`_ and follow the instructions for downloading
   and installing the wambsganss `microlens code here <https://github.com/psaha/microlens/blob/master/wambsganss/README>`_,
   for the `with FITS-file output` option. When you have the `microlens` executable, just replace the SNTD executable in
   `sntd/sntd/microlens` with your new executable, and reinstall. 


Install latest development version
==================================

SNTD is being developed `on github
<https://github.com/sntd>`_. To get the latest development
version using ``git``::

    git clone git://github.com/sntd.git
    cd sntd

then::

    tox
    python setup.py install

`tox` can be install with `pip`, and runs the unit test suite (optional).

Optional dependencies
=====================

Several additional packages are recommended for enabling optional
functionality in SNCosmo. These packages are all pip installable.


- `matplotlib <http://www.matplotlib.org/>`_ for plotting functions.

- `pyParz <https://pypi.org/project/pyParz/>`_ for microlensing uncertainty estimation and fitting large numbers of lensed SN.

- `corner <https://github.com/dfm/corner.py>`_ used for plotting joint and marginalized fitting posteriors

- `iminuit <https://iminuit.readthedocs.io/>`_ used for finding the best model quickly when fitting a list of models


