from setuptools import setup
setup(
    name='sntd',
    install_requires=['cython','numpy','scipy','sncosmo',
    	'astropy','matplotlib','nestle','pyParz','pycs','datetime','iminuit'],
    packages=['sntd'],
    author='Justin Pierel',
    author_email='jr23@email.sc.edu',
    version='0.0.3' 
)
