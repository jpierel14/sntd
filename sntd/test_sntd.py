##Tests for pipeline
import sys,os,traceback
from copy import deepcopy
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import sntd

def test_sntd():
	failed=0
	total=0
	
	try:   
		total+=1 
		print('Testing simulation without microlensing...',end='')
		myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=.5,z_lens=.2, bands=['bessellb','bessellr'],
			  zp=[25,25], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[7,3.5],
			  objectName='My Type Ia SN',telescopename='HST')
		print('Passed!')
	except Exception as e:
		print('Failed')
		print(traceback.format_exc())
		
		failed+=1

	try:    
		total+=1
		print('Testing microlensing realization...',end='')
		myML=sntd.realizeMicro(nray=1,kappas=1,kappac=.3,gamma=.4)
		print('Passed!')
	except Exception as e:
		print('Failed')
		print(traceback.format_exc())
		failed+=1
	
	
	for method in ['parallel','series','color']:
		try:
			total+=1
			print('Testing fitting MISN without microlensing using %s method...'%method,end='')
			fitCurves=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['bessellb','bessellr'],
				params=['x0','x1','t0','c'],constants={'z':.5},bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},
				seriesGrids={'td':(-5,5),'mu':(.8,1.2)},
				refModel=fitCurves.images['image_1'].fits.model if method!='parallel' else None,
				method=method,microlensing=None,maxiter=10,outer_maxiter=1,inner_maxiter=1,outer_npoints=1,inner_npoints=1)
			print('Passed!')
		except Exception as e:
			if method=='parallel':
				print('Failed (this will ruin the next 2 tests)')
			else:
				print('Failed')
			print(traceback.format_exc())
			failed+=1

	try:
		total+=1
		print('Testing simulating MISN with microlensing...',end='')
		myMISN_ml = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=.5,z_lens=.2, bands=['bessellb','bessellr'],
			   zp=[25,25], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[7,3.5],
			   objectName='My Type Ia SN',telescopename='HST', microlensing_type='AchromaticMicrolensing',microlensing_params=myML)
		print('Passed!')
	except Exception as e:
		print('Failed (this will ruin the next test)')
		print(traceback.format_exc())
		failed+=1

	try:
		total+=1

		print('Testing fitting MISN with microlensing using parallel method...',end='')
		fitCurves=sntd.fit_data(myMISN_ml,snType='Ia', models='salt2-extended',bands=['bessellb','bessellr'],
			params=['x0','x1','t0','c'],constants={'z':.5},bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},
			method='parallel',microlensing='achromatic',
			nMicroSamples=1,npoints=1,maxiter=10)

		print('Passed!')
	except Exception as e:
		print('Failed')
		print(traceback.format_exc())
		failed+=1
	try:
		total+=1
		print('Testing example data loading...',end='')
		ex_1,ex_2=sntd.load_example_data()
		new_MISN=sntd.table_factory([ex_1,ex_2],telescopename='HST',object_name='example_SN')
		print('Passed!')
	except Exception as e:
		print('Failed (this will ruin the next test)')
		print(traceback.format_exc())
		failed+=1
	
	try:
		total+=1
		print('Testing fitting a list of MISN...',end='')
		myMISN1 = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=.5,z_lens=.2, bands=['bessellb','bessellr'],
			  zp=[25,25], cadence=3., epochs=35.,time_delays=[10., 30.], magnifications=[10,5],
			  objectName='My Type Ia SN',telescopename='HST')

		myMISN2 = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1,z_lens=.2, bands=['bessellb','bessellr'],
			  zp=[25,25], cadence=3., epochs=35.,time_delays=[10., 20.], magnifications=[10,5],
			  objectName='My Type Ia SN',telescopename='HST')
		
		
		curve_list=[myMISN1,myMISN2]

		fitCurves=sntd.fit_data(curve_list,snType='Ia', models='salt2-extended',bands=['bessellb','bessellr'],
			params=['x0','x1','t0','c'],constants=[{'z':.5},{'z':1}],bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},
			method='parallel',microlensing=None,maxiter=None,npoints=1,verbose=False)
	

	except:
		print('Failed')
		print(traceback.format_exc())
		failed+=1

	print('Passed %i/%i tests.'%(total-failed,total))

	return

if __name__ == '__main__':
	test_sntd()
