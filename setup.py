from setuptools import setup

setup(
    name='sntd',
    setup_requires=['numpy'],
    install_requires=['numpy','scipy','sncosmo',
    	'astropy','matplotlib','nestle','pyParz','sklearn'],
    packages=['sntd'],
    author='Justin Pierel',
    author_email='jr23@email.sc.edu',
    version='1.0.0' 
)
