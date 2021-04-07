from setuptools import setup
import os
import glob
import warnings
import sys
import fnmatch
import subprocess
#from setuptools.command.test import test as TestCommand
from distutils.core import setup
import numpy.distutils.misc_util


if sys.version_info < (3, 0):
    sys.exit('Sorry, Python 2 is not supported')


AUTHOR = 'Justin Pierel'
AUTHOR_EMAIL = 'jr23@email.sc.edu'
VERSION = '2.5.5'
LICENSE = 'BSD'
URL = 'sntd.readthedocs.org'


def recursive_glob(basedir, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(basedir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches


PACKAGENAME = 'sntd'
# Add the project-global data
pkgdatadir = os.path.join(PACKAGENAME, 'data')
microdatadir = os.path.join(PACKAGENAME, 'microlens')
simdatadir = os.path.join(PACKAGENAME, 'sim')
batchdatadir = os.path.join(PACKAGENAME, 'batch')
data_files = []
data_files.extend(recursive_glob(pkgdatadir, '*'))
data_files.extend(recursive_glob(microdatadir, '*'))
data_files.extend(recursive_glob(simdatadir, '*'))
data_files.extend(recursive_glob(batchdatadir, '*'))

data_files = [f[len(PACKAGENAME)+1:] for f in data_files]


setup(
    name='sntd',
    setup_requires=['numpy', 'cython'],
    install_requires=['numpy', 'scipy', 'cython', 'sncosmo',
                      'astropy', 'matplotlib', 'nestle', 'pyParz', 'sklearn',
                      'iminuit==1.4.9', 'corner', 'pandas'],
    packages=['sntd'],
    version=VERSION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    package_data={'sntd': data_files},
    include_package_data=True
)
