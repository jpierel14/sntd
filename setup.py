from setuptools import setup
import os,glob,warnings,sys,fnmatch
from setuptools.command.test import test as TestCommand
from distutils.core import setup
import numpy.distutils.misc_util


if sys.version_info < (3,0):
    sys.exit('Sorry, Python 2 is not supported')

class SNTDTest(TestCommand):

    def run_tests(self):
        import sntd
        errno = sntd.test()
        sntd.test_sntd()
        sys.exit(errno)

AUTHOR = 'Justin Pierel'
AUTHOR_EMAIL = 'jr23@email.sc.edu'
VERSION = '1.0.9'
LICENSE = 'BSD'
URL = 'sntd.readthedocs.org'

def recursive_glob(basedir, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(basedir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches

PACKAGENAME='sntd'
# Add the project-global data
pkgdatadir = os.path.join(PACKAGENAME, 'data')
microdatadir = os.path.join(PACKAGENAME, 'microlens')
simdatadir = os.path.join(PACKAGENAME, 'sim')
data_files = []
data_files.extend(recursive_glob(pkgdatadir, '*'))
data_files.extend(recursive_glob(microdatadir, '*'))
data_files.extend(recursive_glob(simdatadir, '*'))

data_files = [f[len(PACKAGENAME)+1:] for f in data_files]


setup(
    name='sntd',
    cmdclass={'test': SNTDTest},
    setup_requires=['numpy'],
    install_requires=['numpy','scipy','sncosmo',
    	'astropy','matplotlib','nestle','pyParz','sklearn',
        'pytest-astropy'],
    packages=['sntd'],
    version=VERSION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    package_data={'sntd':data_files},
    include_package_data=True
)
