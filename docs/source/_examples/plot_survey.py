"""
===================
Constrain Cosmology
===================

Simulate cosmological constraints
from a sample of lensed SN.
"""
		
###################################################################
# All of these cosmology tools are based on `Coe & Moustakas 2009 <https://arxiv.org/pdf/0906.4108.pdf>`_.
# and Dan Coe's `Fisher matrix starter paper <https://arxiv.org/pdf/0906.4123.pdf>`_.
#
# **Creating a Survey**
import sntd
import numpy as np

#####################################################################
# Start by defining your survey parameters. In this case we have a survey called "Test Survey" with
# 10 lenses with normally distributed lens and source redshifts,
# 5% lens model uncertainty and 2% time delay uncertainty.

np.random.seed(3)

my_survey=sntd.Survey(dTl=5,dTT=2,zl=np.random.normal(.5,.1,size=10),zs=np.random.normal(1.6,.2,size=10),name='Test Survey')

#####################################################################
# **Gridded Parameter Search**
#
# This will make a smooth contour plot for 2 parameters.

my_survey.survey_grid(vparam_names=['h','Ode0'],
                      bounds={'h':[.65,.75],'Ode0':[0,1]},npoints=50)

my_survey.plot_survey_contour(['h','Ode0'],math_labels=[r'$h$',r'$\Omega_\lambda$'],confidence=[.68,.95],alphas=[.9,.4],show_legend=True)

#####################################################################
# **MCMC-Like Parameter Search**

my_survey.survey_nestle(vparam_names=['h','Ode0'],
                      bounds={'h':[.65,.75],'Ode0':[0,1]},npoints=200)

my_survey.plot_survey_contour(['h','Ode0'],math_labels=[r'$h$',r'$\Omega_\lambda$'],filled=False)

#####################################################################
# **Fisher Matrix Analysis**
# 
# This will make a 5x5 fisher matrix with the given parameters

my_survey.survey_fisher(['h','Ode0','Om0','w0','wa'])

#####################################################################
# Add a prior that assumes perfect knowledge of all other parameters

my_survey.fisher_matrix.prior('Om0',0.0001)
my_survey.fisher_matrix.prior('Ode0',0.0001)
my_survey.fisher_matrix.prior('h',0.0001)
my_survey.fisher_matrix.plot('w0','wa',x_limits=[-1.7,-.3],y_limits=[-4,4])


