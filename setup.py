from setuptools import setup
import os,glob,warnings,sys,fnmatch

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
    setup_requires['numpy'],
    install_requires=['numpy','scipy','sncosmo',
    	'astropy','matplotlib','nestle','pyParz','sklearn'],
    packages=['sntd'],
    author='Justin Pierel',
    author_email='jr23@email.sc.edu',
    version='1.0.1',
    package_data={'sntd':data_files}
)
