import sntd
import matplotlib.pyplot as plt


myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=1.33,z_lens=.53, bands=['F110W','F125W'],
                                     zp=[26.8,26.2], cadence=5., epochs=35.,time_delays=[10., 70.], magnifications=[7,3.5],
                                     objectName='My Type Ia SN',telescopename='HST')

print(myMISN)

fitCurves=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F125W'],
                        params=['x0','x1','t0','c'],constants={'z':1.33},bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},
                        method='separate',microlensing='Achromatic')


fitCurves.plot_object(showFit=True,method='separate')
plt.show()