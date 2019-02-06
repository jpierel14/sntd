************
Installation
************

SNSEDextend works on Python 2.7 and Python 3.4+ and requires the
following Python packages:

- `numpy <http://www.numpy.org/>`_
- `scipy <http://www.scipy.org/>`_
- `astropy <http://www.astropy.org>`_
- `SNCosmo <http://sncosmo.readthedocs.io>`_

Install using pip
=================

Using pip::

    pip install snsedextend

.. note::

    You will need a C compiler (e.g. ``gcc`` or ``clang``) to be
    installed for the installation to succeed due to SNCosmo.

Setting Environment Variable
============================
If you have not already, you will need to define the path to your SNDATA_ROOT folder by setting it as an environment variable. If you do not, the package will set it to your current directory. On Mac or Linux, set it by adding the following line to your .bashrc file::
  
  export SNDATA_ROOT='full/path/to/folder'

On Windows, follow the directions here_.

.. _here: http://www.dowdandassociates.com/blog/content/howto-set-an-environment-variable-in-windows-command-line-and-registry/

Install latest development version
==================================

SNSEDextend is being developed `on github
<https://github.com/snsed>`_. To get the latest development
version using ``git``::

    git clone git://github.com/snsed.git
    cd snsed

then::

    ./setup.py install


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
- `nestle <http://kbarbary.github.io/nestle/>`_ for nested sampling
  light curve parameter estimation in `sncosmo.nest_lc`.

iminuit, emcee and nestle can be installed using pip.

The `corner <https://github.com/dfm/corner.py>`_ package is also
recommended for plotting results from the samplers `sncosmo.mcmc_lc`
and `sncosmo.nest_lc`, but is not used by any part of sncosmo.
