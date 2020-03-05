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
			  zp=[25,25], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[20,20],
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
				params=['x0','x1','t0','c'],constants={'z':.5},bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1),'td':(-15,15),'mu':(.5,2)},
				method=method,microlensing=None,maxcall=5,minsnr=0,t0_guess={'image_1':10,'image_2':70})
			print('Passed!')
		except Exception as e:
			print('Failed')
			print(traceback.format_exc())
			failed+=1

	try:
		total+=1
		print('Testing simulating MISN with microlensing...',end='')
		myMISN_ml = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=.5,z_lens=.2, bands=['bessellb','bessellr'],
			   zp=[25,25], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[20,20],
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
			method='parallel',microlensing='achromatic',t0_guess={'image_1':10,'image_2':70},
			nMicroSamples=1,maxcall=5,minsnr=0)

		print('Passed!')
	except Exception as e:
		print('Failed')
		print(traceback.format_exc())
		failed+=1
	try:
		total+=1
		
		print('Testing fitting MISN with microlensing using series method...',end='')
		
		fitCurves=sntd.fit_data(myMISN_ml,snType='Ia', models='salt2-extended',bands=['bessellb','bessellr'],
			params=['x1','c'],constants={'z':.5},bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1),'td':(-15,15),'mu':(.5,2)},
			method='series',microlensing='achromatic',t0_guess={'image_1':10,'image_2':70},
			nMicroSamples=1,maxcall=5,minsnr=0)

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
			  zp=[25,25], cadence=3., epochs=35.,time_delays=[10., 30.], magnifications=[20,20],
			  objectName='My Type Ia SN',telescopename='HST')

		myMISN2 = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1,z_lens=.2, bands=['bessellb','bessellr'],
			  zp=[25,25], cadence=3., epochs=35.,time_delays=[10., 20.], magnifications=[20,20],
			  objectName='My Type Ia SN',telescopename='HST')
		
		
		curve_list=[myMISN1,myMISN2]

		fitCurves=sntd.fit_data(curve_list,snType='Ia', models='salt2-extended',bands=['bessellb','bessellr'],
			params=['x0','x1','t0','c'],constants=[{'z':.5},{'z':1}],t0_guess={'image_1':10,'image_2':70},
			bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},
			method='parallel',microlensing=None,maxcall=5,verbose=False,minsnr=0)
	
		print('Passed!')
	except:
		print('Failed')
		print(traceback.format_exc())
		failed+=1

	print('Passed %i/%i tests.'%(total-failed,total))

	return

if __name__ == '__main__':
	test_sntd()
