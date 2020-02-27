import sntd,sncosmo,sys,pickle
import matplotlib.pyplot as plt
import numpy as np


z=1.4
myMISN = sntd.createMultiplyImagedSN(sourcename='salt2-extended', snType='Ia', redshift=z,z_lens=.53, bands=['F110W','F160W'],
                                     zp=[26.8,26.2], cadence=10., epochs=20.,time_delays=[20., 80.], magnifications=[4,8],
                                     objectName='My Type Ia SN',telescopename='HST',scatter=True,minsnr=5)

#sncosmo.plot_lc(myMISN.images['image_1'].table)
#plt.show()
#sncosmo.plot_lc(myMISN.images['image_2'].table)
#plt.show()

# print(myMISN.images['image_1'].simMeta['model'].parameters)
#
# fitCurves=sntd.fit_data(myMISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
#                         params=['c'],constants={'z':z},
#                         bounds={'td':(-20,20),'mu':(.5,2),'x1':(-3,3),'c':(-.5,.5)},
#                         method='color',microlensing=None,modelcov=False,npoints=10,maxiter=None)
from astropy.table import Table
from timeit import default_timer as timer

inds=[1,3]
all_dat=pickle.load(open('/Users/jpierel/Downloads/all_Ia_AP_sim.dat','rb'),encoding='bytes')
for i in range(inds[0],inds[1]):
    im1=all_dat[i][0]
    im2=all_dat[i][1]
    tab1=Table(im1[:-2],names=('time','band','flux','fluxerr'))
    tab2=Table(im2[:-2],names=('time','band','flux','fluxerr'))
    tab1['filter']='F160W'
    tab2['filter']='F160W'
    for j in range(len(tab1)):

        tab1['filter'][j]='F160W' if '6' in tab1['band'][j] else 'F110W'
    for j in range(len(tab2)):
        tab2['filter'][j]='F160W' if '6' in tab2['band'][j] else 'F110W'
    tab1.remove_column('band')
    tab2.remove_column('band')

    tab1['image']=im1[-2]

    tab2['image']=im2[-2]
    tab1['zp']=0
    tab2['zp']=0
    tab1['zpsys']='ab'
    tab2['zpsys']='ab'
    tab1['zp'][tab1['filter']=='F110W']=26.82
    tab2['zp'][tab2['filter']=='F110W']=26.82
    tab1['zp'][tab1['filter']=='F160W']=25.94
    tab2['zp'][tab2['filter']=='F160W']=25.94

    tab1['fluxerr']*=5
    tab2['fluxerr']*=5
    new_MISN=sntd.table_factory([tab1,tab2],telescopename='HST',object_name='example_SN')

    # print(fitCurves.color.table.colnames)
    # print(fitCurves.color.fits.model.parameters)
    # print(fitCurves.color.time_delays)
    # print(fitCurves.color.time_delay_errors)
    # fitCurves.plot_fit(method='color')
    # plt.show()
    # plt.errorbar(fitCurves.color.table['time'],fitCurves.color.table['F110W-F160W'],yerr=fitCurves.color.table['F110W-F160W_err'],fmt='.')
    # fitCurves.color_table('F110W','F160W',time_delays={'image_1':0,'image_2':60})
    # plt.errorbar(fitCurves.color.table['time'],fitCurves.color.table['F110W-F160W'],yerr=fitCurves.color.table['F110W-F160W_err'],fmt='.')
    # plt.plot(np.arange(-30,70,1),fitCurves.color.fits.model.color('F110W','F160W','ab',np.arange(-30,70,1)))
    # plt.show()

    start=timer()

    fitCurves=sntd.fit_data(new_MISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
                            params=['x0','t0','x1','c'],constants={'z':1.4},refImage='image_2',
                            bounds={'t0':(-40,20),'x1':(-3,3),'c':(-1,1)},
                            method='parallel',microlensing=None,modelcov=False,npoints=100,maxiter=None)

    for k in fitCurves.images.keys():
        print(fitCurves.images[k].param_quantiles)
    print('series:',timer()-start)
    print(fitCurves.time_delays)
    print(fitCurves.time_delay_errors)
    print(fitCurves.magnifications)
    print(fitCurves.magnification_errors)
    fitCurves.plot_fit(method='parallel',par_image='image_1')
    plt.show()
    fitCurves.plot_fit(method='parallel',par_image='image_2')
    plt.show()
    sncosmo.plot_lc(fitCurves.images['image_2'].table,model=fitCurves.images['image_2'].fits.model)
    plt.show()
    sncosmo.plot_lc(fitCurves.images['image_1'].table,model=fitCurves.images['image_1'].fits.model)

    #res=fitCurves.images['image_1'].fits.res)
    plt.show()

    start=timer()
    fitCurves=sntd.fit_data(new_MISN,snType='Ia', models='salt2-extended',bands=['F110W','F160W'],
                            params=['x1','c'],constants={'z':1.4},refImage='image_2',
                            bounds={'td':(-40,20),'mu':(.5,2),'x1':(-3,3),'c':(-1,1)},fit_prior=fitCurves,
                            method='series',microlensing=None,modelcov=False,npoints=200,maxiter=None)

    print('parallel:',timer()-start)
    print(fitCurves.series.time_delays)
    print(fitCurves.series.time_delay_errors)
    print(fitCurves.series.magnifications)
    print(fitCurves.series.magnification_errors)
    fitCurves.plot_fit(method='series')
    plt.show()
    sncosmo.plot_lc(fitCurves.series.table,model=fitCurves.series.fits.model)
    plt.show()
    sys.exit()