import sntd,sncosmo,sys
import matplotlib.pyplot as plt


'''
myObj=sntd.curveDict(telescopename='HST',object='Refsdal')
print(myObj.table)
curve=sntd.read_data('testDat.dat',band='F125W',zp=26.2,zpsys='AB')
myObj.add_curve(curve)
print(myObj)

#fitData=sntd.realizeMicro(nray=10,kappas=1,kappac=.3,gamma=.4)
'''
fitData=sntd.realizeMicro(nray=10,kappas=1,kappac=.3,gamma=.4)
modname = sncosmo.SALT2Source(modeldir='salt2-extended')
myMISN = sntd.createMultiplyImagedSN(modname, 'Ia', 1.33,z_lens=.53, bands=['F110W','F125W'],
                                     zp=[26.8,26.2], cadence=5., epochs=35.,skynoiseRange=(.001,.005),gain=70. ,
                                     time_delays=[10., 78.], magnifications=[7,3.5], objectName='My Type Ia SN',
                                     telescopename='HST',minsnr=5.0,microlensing_type='AchromaticMicrolensing',
                                     microlensing_params=fitData)

fitCurves=sntd.fit_data(myMISN,snType='Ia', models=modname,bands=['F110W','F125W'],
                            params=['x0','x1','t0','c'],constants={'z':1.33},
                            bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},showPlots=False,method='separate')

fitCurves=sntd.fit_data(fitCurves,snType='Ia', models=modname,bands=['F110W','F125W'],
                        params=['x0','x1','t0','c'],constants={'z':1.33},combinedGrids={'mu':(-.2,.2),'td':(-10,10)},
                        bounds={'x1':(-2,2),'c':(0,1)},showPlots=False,method='color')
fig=fitCurves.plot_object(savefig=False,showModel=False,method='color',showFit=True,showMicro=False)
plt.show()
sys.exit()
fitCurves=sntd.fit_data(fitCurves,snType='Ia', models=modname,bands=['F110W'],
                        params=['x0','x1','t0','c'],constants={'z':1.33},combinedGrids={'mu':(-.2,.2),'td':(-10,10)},
                        bounds={'x1':(-2,2),'c':(0,1)},showPlots=False,method='combined')

#myMISN.combine_curves(time_delays={'image_1':10,'image_2':78},magnifications={'image_1':7,'image_2':3.5})


