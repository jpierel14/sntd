##Tests for pipeline
import sys,os,traceback,shutil,unittest
from copy import deepcopy
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import sntd
import numpy as np
import warnings

warnings.simplefilter('ignore')

_NOSBATCH_ = True
_GOFAST_ = True
_PARONLY_ = False

np.random.seed(3)
class TestMicrolensing(unittest.TestCase):
	"""
	Test SNTD microlensing simulation tools
	"""
	def setUp(self):
		self.myML=sntd.realizeMicro(nray=10,kappas=1,kappac=.3,gamma=.4)
		self.myMISN_ml = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=.5,z_lens=.2, 
			bands=['bessellb','bessellr'],
			   zp=[25,25], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[20,20],
			   objectName='My Type Ia SN',telescopename='HST', microlensing_type='AchromaticMicrolensing',
			   microlensing_params=self.myML)
		
	@unittest.skipIf(_GOFAST_,"Skipping slow `test_fit_lc_Micro_parallel`")
	def test_fit_lc_Micro_parallel(self):
		fitCurves=sntd.fit_data(self.myMISN_ml,snType='Ia', models='salt2-extended',bands=['bessellb','bessellr'],
			params=['x0','x1','t0','c'],constants={'z':.5},bounds={'t0':(-1,1),'x1':(-2,2),'c':(-1,1)},
			method='parallel',microlensing='achromatic',t0_guess={'image_1':10,'image_2':70},
			nMicroSamples=5,maxcall=None,npoints=50,minsnr=0)

	@unittest.skipIf(_GOFAST_,"Skipping slow `test_fit_lc_Micro_series`")
	def test_fit_lc_Micro_series(self):
		fitCurves=sntd.fit_data(self.myMISN_ml,snType='Ia', models='salt2-extended',bands=['bessellb','bessellr'],
			params=['t0','x0','x1','c'],constants={'z':.5},bounds={'t0':(-15,15),'x1':(-2,2),'c':(-1,1),'td':(-15,15),'mu':(.5,2)},
			method='series',microlensing='achromatic',t0_guess={'image_1':10,'image_2':70},
			nMicroSamples=5,maxcall=100,npoints=25,minsnr=0)


class TestSimulation(unittest.TestCase):
	"""
	Test SNTD light curve simulation tools
	"""
	def setUp(self):
		self.myML=sntd.realizeMicro(nray=10,kappas=1,kappac=.3,gamma=.4)

	def test_sim_lc_noMicro(self):
		myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=.5,z_lens=.2, 
			bands=['bessellb','bessellv','bessellr'],
			zp=[25,25,25], cadence=5., epochs=20.,time_delays=[20., 70.], magnifications=[5,5],
			objectName='My Type Ia SN',telescopename='HST')

	

class TestFitting(unittest.TestCase):
	"""
	Test SNTD light curve fitting tools
	"""
	def setUp(self):
		self.myMISN = sntd.load_example_misn()
		# self.myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=.5,z_lens=.2, 
		# 	bands=['F110W','F160W'],
		# 	zp=[25,25], cadence=5., epochs=13.,time_delays=[20., 70.], magnifications=[5,5],
		# 	objectName='My Type Ia SN',telescopename='HST')

	@unittest.skipIf(_PARONLY_,"Skipping non-parallel fit.")
	def test_quality_check(self):
		for method in ['parallel','series','color']:
			fitCurves=sntd.fit_data(self.myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
				params=['x0','x1','t0','c'],bounds={'t0':(-15,15),'x1':(-2,2),'c':(-1,1),'td':(-15,15),'mu':(.5,2)},
				color_param_ignore=['x1'],min_n_bands=1000,min_points_per_band=10000,
				method=method,microlensing=None,maxcall=50,minsnr=0,set_from_simMeta={'z':'z'},
				t0_guess={'image_1':20,'image_2':70})

	def test_parallel_fit(self):
		fitCurves=sntd.fit_data(self.myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
				params=['x0','x1','t0','c'],bounds={'t0':(-15,15),'x1':(-2,2),'c':(-1,1),'td':(-30,30),'mu':(.5,2)},
				color_param_ignore=['x1'],use_MLE=False,refImage='image_1',cut_time=[-40,20],
				method='parallel',microlensing=None,maxcall=None,npoints=25,minsnr=0,
				set_from_simMeta={'z':'z'},t0_guess={'image_1':20,'image_2':65})

	@unittest.skipIf(_PARONLY_,"Skipping non-parallel fit.")
	def test_series_fit(self):
		fitCurves=sntd.fit_data(self.myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
				params=['x0','x1','t0','c'],bounds={'t0':(-15,15),'x1':(-2,2),'c':(-1,1),'td':(-30,30),'mu':(.5,2)},
				color_param_ignore=['x1'],use_MLE=False,refImage='image_1',
				method='series',microlensing=None,maxcall=100,npoints=10,minsnr=0,
				set_from_simMeta={'z':'z'},t0_guess={'image_1':20,'image_2':70})

	@unittest.skipIf(_PARONLY_,"Skipping non-parallel fit.")
	def test_color_fit(self):
		fitCurves=sntd.fit_data(self.myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
				params=['x0','x1','t0','c'],bounds={'t0':(-15,15),'x1':(-2,2),'c':(-1,1),'td':(-30,30),'mu':(.5,2)},
				color_param_ignore=['x1'],use_MLE=False,refImage='image_1',
				method='color',microlensing=None,maxcall=100,npoints=10,minsnr=0,
				set_from_simMeta={'z':'z'},t0_guess={'image_1':20,'image_2':70})
class TestCosmology(unittest.TestCase):
	"""
	Test SNTD cosmology tools.
	"""
	def setUp(self):
		test_cosmo_dict = {
			'N':1,   # number of Lensed SNe Ia with good time delays
			'dTL':2,  # % lens modeling uncertainty for each
			'dTT':.1,  # % time delay measurement uncertainty for each
			'zl':.5,'zs':2
		}
	
		self.test_cosmo = sntd.Survey(**test_cosmo_dict)

	def test_survey_grid(self): 
		self.test_cosmo.survey_grid(['w','Ode0'],{'w':[-1.5,-.5],'Ode0':[0,1]},npoints=2)

	def test_survey_nestle(self):
		self.test_cosmo.survey_nestle(['w','Ode0'],{'w':[-1.5,-.5],'Ode0':[0,1]},npoints=2)

	def test_survey_fisher(self):
		self.test_cosmo.survey_fisher(['w','Ode0'])

class TestBatch(unittest.TestCase):
	def setUp(self):
		self.myMISN = sntd.load_example_misn()

	def test_multiprocessing(self):
		fitCurves=sntd.fit_data([self.myMISN]*5,snType='Ia', models='salt2-extended',bands=['bessellb','bessellr'],
				params=['x0','x1','t0','c'],constants=[{'z':.5}]*5,
				bounds={'t0':(-15,15),'x1':(-2,2),'c':(-1,1),'td':(-15,15),'mu':(.5,2)},
				method='parallel',microlensing=None,maxcall=50,npoints=10,minsnr=0,t0_guess={'image_1':10,'image_2':70},
				verbose=False)

	@unittest.skipIf(_NOSBATCH_,"Skipping sbatch test")
	def test_sbatch(self):
		fitCurves=sntd.fit_data([self.myMISN]*100,snType='Ia', models='salt2-extended',bands=['bessellb','bessellr'],
				params=['x0','x1','t0','c'],constants={'z':.5},
				bounds={'t0':(-15,15),'x1':(-2,2),'c':(-1,1),'td':(-15,15),'mu':(.5,2)},
				method='parallel',wait_for_batch=False,
				par_or_batch='batch',nbatch_jobs=2,microlensing=None,maxcall=5,minsnr=0,t0_guess={'image_1':10,'image_2':70})

	@unittest.skipIf(_NOSBATCH_,"Skipping sbatch test")
	def tearDown(self):
		try:
			shutil.rmtree('batch_output')
		except RuntimeError:
			pass


def test_loader(loader):
    suite = unittest.TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite


if __name__ == '__main__':
    #TEST LIST
    test_cases = 'ALL'
    #test_cases = [TestFitting]

    if test_cases == 'ALL':
        unittest.main()
    else:
        runner = unittest.TextTestRunner()
        runner.run(test_loader(unittest.TestLoader()))
