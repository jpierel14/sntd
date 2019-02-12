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

.. note::

    You will need a C compiler (e.g. ``gcc`` or ``clang``) to be
    installed for the installation to succeed due to SNCosmo.

.. note::

   Sometimes there is an issue with the install related to numpy. If this happens, try installing numpy using pip, and then installing SNTD using pip.
   


Install latest development version
==================================

SNTD is being developed `on github
<https://github.com/sntd>`_. To get the latest development
version using ``git``::

    git clone git://github.com/sntd.git
    cd sntd

then::

    ./setup.py install


Optional dependencies
=====================

.. note::
   If you plan to use SNTD to simulate the effects of microlensing, which uses the Wambsganss 1990 microlens code, then you will need a fortran compiler. You can install gfortran `here
   <https://gcc.gnu.org/wiki/GFortranBinaries>`_.

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
