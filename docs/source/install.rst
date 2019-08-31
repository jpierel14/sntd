************
Installation
************

SNTD works on Python 3.5+ and requires the
following Python packages:

- `numpy <http://www.numpy.org/>`_
- `scipy <http://www.scipy.org/>`_
- `astropy <http://www.astropy.org>`_
- `SNCosmo <http://sncosmo.readthedocs.io>`_
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

Install latest development version
==================================

SNTD is being developed `on github
<https://github.com/sntd>`_. To get the latest development
version using ``git``::

    git clone git://github.com/sntd.git
    cd sntd

then::

    python setup.py test
    python setup.py install


Optional dependencies
=====================

Several additional packages are recommended for enabling optional
functionality in SNCosmo.

- `matplotlib <http://www.matplotlib.org/>`_ for plotting
  functions.
- `iminuit <http://iminuit.github.io/iminuit/>`_ for light curve
  fitting using the Minuit minimizer in `sncosmo.fit_lc`.
- `emcee <http://dan.iel.fm/emcee/>`_ for MCMC light curve parameter
  estimation in `sncosmo.mcmc_lc`.
- `pyParz <https://pypi.org/project/pyParz/>`_ for microlensing uncertainty estimation.

These packages are all pip installable.

The `corner <https://github.com/dfm/corner.py>`_ package is also
recommended for plotting results from the samplers `sncosmo.mcmc_lc`
and `sncosmo.nest_lc`, but is not used by any part of sncosmo.
