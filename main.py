import sntd,sncosmo,sys,os
import matplotlib.pyplot as plt
import snsedextend
import numpy as np
from astropy.table import Table
from astropy.io import ascii

def get_curve(temp_table,band1,band2,band1zp,band2zp):
    temp_table['time']=temp_table['time'].astype(int)
    common_times=np.unique(temp_table['time'])
    all_dat=[]
    for t in common_times:
        if len(temp_table[temp_table['time']==t])==2:
            all_dat.append(temp_table[temp_table['time']==t])
    color_table=Table(names=['time','B-V','BV_err'],masked=True)
    for t in all_dat:
        row=[t['time'][0],-2.5*np.log10(t['flux'][0]/t['flux'][1])+band1zp-band2zp,
             np.sqrt((1.0857*t['fluxerr'][0]/t['flux'][0])**2+(1.0857*t['fluxerr'][1]/t['flux'][1])**2)]
        color_table.add_row(row)

    temp_time=np.min(color_table['time'])
    color_table['time']-=temp_time
    curveDict=snsedextend.fitColorCurve(color_table)
    curveDict['B-V']['time']+=temp_time
    curveDict['B-V'].rename_column('B-V',band1+'-'+band2)
    color_table['time']+=temp_time
    return(curveDict['B-V'],color_table)

ex_1,ex_2=sntd.load_example_data()
new_MISN=sntd.table_factory([ex_1,ex_2],telescopename='HST',object_name='example_SN')
fitCurves=sntd.fit_data(new_MISN,snType='Ia', models='salt2-extended',bands=['F125W','F160W'],
                                                     params=['x0','x1','t0','c'],constants={'z':1.33},
                                                     bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)})
print(fitCurves.time_delays)
fitCurves.plot_object(showFit=True)
plt.savefig('example_fit.png',format='png',overwrite=True)
sys.exit()

myMISN = sntd.createMultiplyImagedSN('salt2-extended', 'Ia', 1.33,z_lens=.53, bands=['F125W','F160W'],
                                          zp=[26.8,26.2], cadence=2., epochs=90.,
                                          time_delays=[10., 70.], magnifications=[7,3.5], objectName='My Type Ia SN',
                                          telescopename='HST',minsnr=5.0,microlensing_type=None)

for im in myMISN.images:
    ascii.write(myMISN.images[im].table,os.path.join(sntd.util.__dir__,'data','examples','example_'+im+'.dat'))
sys.exit()
curve,tab=get_curve(myMISN.images['image_2'].table,'F125W','F160W',26.8,26.2)
fig=plt.figure()
ax=fig.gca()
ax.plot(curve['time'],curve['F125W-F160W'])
source=sntd.models.BazinSource(myMISN.images['image_2'].table,
                               colorCurve=curve)
mod1=sncosmo.Model('salt2-extended')
mod1.set(z=1.33)
mod1.set(t0=70)
mod2=sncosmo.Model(source)

fig2=plt.figure()
ax2=fig2.gca()
ax2.scatter(tab['time'],tab['B-V'])
time=np.arange(tab['time'].min(),tab['time'].max())
ax2.plot(time,mod1.color('F125W','F160W','ab',time))
ax2.plot(time,mod2.color('F125W','F160W','ab',time))
plt.show()

'''
myObj=sntd.curveDict(telescopename='HST',object='Refsdal')
print(myObj.table)
curve=sntd.read_data('testDat.dat',band='F125W',zp=26.2,zpsys='AB')
myObj.add_curve(curve)
print(myObj)

#fitData=sntd.realizeMicro(nray=10,kappas=1,kappac=.3,gamma=.4)
'''
# fitData=sntd.realizeMicro(nray=10,kappas=1,kappac=.3,gamma=.4)
# modname = sncosmo.SALT2Source(modeldir='salt2-extended')
# myMISN = sntd.createMultiplyImagedSN('salt2-extended', 'Ia', 1.33,z_lens=.53, bands=['F110W','F125W'],
#                                      zp=[26.8,26.2], cadence=5., epochs=35.,
#                                      time_delays=[10., 70.], magnifications=[7,3.5], objectName='My Type Ia SN',
#                                      telescopename='HST',minsnr=5.0,microlensing_type='AchromaticMicrolensing',
#                                      microlensing_params=fitData)
#
# myMISN.plot_object()
# plt.show()
# fitCurves=sntd.fit_data(myMISN,snType='Ia', models=modname,bands=['F110W','F125W'],
#                             params=['x0','x1','t0','c'],constants={'z':1.33},
#                             bounds={'t0':(-15,15),'x1':(-2,2),'c':(0,1)},showPlots=False,method='separate')
#
# fitCurves.plot_object(method='separate',showFit=True)
# fitCurves2=sntd.fit_data(fitCurves,snType='Ia', models=modname,bands=['F110W','F125W'],refImage='image_2',
#                         params=['x0','x1','t0','c'],constants={'z':1.33},combinedGrids={'mu':(-.2,.2),'td':(-10,10)},
#                         bounds={'x1':(-2,2),'c':(0,1)},showPlots=False,method='color')
#
# fig=fitCurves2.plot_object(savefig=False,showModel=False,method='color',showFit=True,showMicro=False)
# plt.show()
#
# fitCurves3=sntd.fit_data(fitCurves,snType='Ia', models=modname,bands=['F110W'],refImage='image_2',
#                         params=['x0','x1','t0','c'],constants={'z':1.33},combinedGrids={'mu':(-.2,.2),'td':(-10,10)},
#                         bounds={'x1':(-2,2),'c':(0,1)},showPlots=False,method='combined')
# fig=fitCurves3.plot_object(savefig=False,showModel=False,method='combined',showFit=True,showMicro=False)
# plt.show()
#myMISN.combine_curves(time_delays={'image_1':10,'image_2':78},magnifications={'image_1':7,'image_2':3.5})


