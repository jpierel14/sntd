'''

import sncosmo
import pickle
data=sncosmo.load_example_data()
model=sncosmo.Model('salt2')
res, fitted_model = sncosmo.fit_lc(data, model,['z', 't0', 'x0', 'x1', 'c'],bounds={'z':(0.3, 0.7)})
#model._source=None
with open('test.pkl','wb') as handle:
    pickle.dump(sncosmo.Model('salt2'),handle,protocol=2)

'''
import timeit,glob,sys,warnings,sncosmo
from multiprocessing import Pool
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter as smooth
from scipy.interpolate import UnivariateSpline as spl
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import vstack
from sntd import io,fitting,simulation
warnings.simplefilter('ignore')
'''
im,curve=ml.realizeMicro(arand=.243,kappas=.2,kappac=0,gamma=0,eps=.6,nray=100,minmass=1,maxmass=1,pixmax=10,pixminx=-10,pixminy=-10,pixdif=20,fracpixd=.15)
fig=plt.figure()
ax=fig.gca()
ax.plot(curve['t'],curve['maglin'])
ax.set_ylim([0,6])
plt.savefig('curve.pdf',format='pdf',overwrite=True)
plt.close()
fig=plt.figure()
cax=plt.imshow(im)
cbar = fig.colorbar(cax, ticks=[np.linspace(1.1*np.min(im),np.max(im),5)])
cbar.ax.set_yticklabels([str(np.round(x,2)) for x in np.linspace((1.1*np.min(im)-1024)/256.,(np.max(im)-1024)/256.,5)])  # vertically oriented colorbar
plt.savefig('lensplane.pdf',format='pdf',overwrite=True)
sys.exit()
'''
#from sntd import models

#dat=sncosmo.load_example_data()
files=glob.glob('data/ref*.dat')
#files=['data/refsdalS1_psfphot.dat']
curves=io.curveDict(telescopename='Hubble',object='Refsdal')
final=None

minFluxes=[]
maxFluxes=[]
for f in files:
    tab=io.read_data(f)
    tab.table=tab.table[tab.table.mask['band']==False]
    tab.table=tab.table[tab.table['band']=='F105W']
    minFluxes.append(np.min(tab.table['flux']))
    maxFluxes.append(np.max(tab.table['flux']))
    #tab.table=tab.table[tab.table['time']<np.min(tab.table['time'])+300]
    curves.add_curve(tab)
    #tab=curves.images[f[-14:-12]].table[curves.images[f[-14:-12]].table.mask['band']==False]
    #tab=tab[tab['band']=='F160W']
    #if not final:
    #    final=tab
    #else:
    #    final=vstack(final,tab)
#print(tab['band'])

#tab=curves.images['S1'].table[curves.images['S1'].table.mask['band']==False]
#tab=tab[tab['band']!='F606W']
#tab=tab[tab['band']!='F140W']
#tab=tab[tab['band']=='F160W']
#print(tab)
#temp=sntd.models.FlexibleSpline(tab,knots=2)
#mod=sncosmo.Model(temp)
#mod.set(t0=0)

#mod='snana-2004gv'
#snType='Ib'
#curves=simulation.createRandMultiplyImagedSN(mod,snType,.1,bands=['bessellv'],zp=25,cadence=10,epochs=15,numImages=4,timeDelayRange=(0,30),objectName='Test',telescopename='HST',microlensing=False)
#lcs=fitting.fit_data(curves,snType='II',models=['SplineSource'],constants={'t0':0},bounds={'dt0':(-5,5),'amplitude':(.99,1.01)},func='spline',combined_or_separate='separate')
#lcs=fitting.fit_data(curves,snType='II',models=['BazinSource'],params=['t0','A','B','fall','rise'],bounds={'A':(0,1000),'B':(-100,100),'fall':(40,100),'rise':(30,60)},combined_or_separate='separate')
#lcs=fitting.fit_data(curves,snType='II',models=['KarpenkaSource'],params=['A','B','fall','rise','t1','t0'],bounds={'A':(0,1000),'B':(0,100),'fall':(5,100),'rise':(0,100),'t1':(-10,10)},combined_or_separate='separate')
lcs=fitting.fit_data(curves,snType='II',models=['NewlingSource'],params=['A','psi','sigma','k','t0','phi'],bounds={'A':(np.min(minFluxes),np.max(maxFluxes)),'psi':(np.min(minFluxes),np.max(maxFluxes)),'sigma':(1,200),'k':(.01,5)},combined_or_separate='separate')


#lcs.plot_object()
#plt.show()
#res,fit=sncosmo.fit_lc(tab,mod,['dt0','amplitude',],bounds={'dt0':(-5,5)},guess_t0=False,guess_amplitude=False)

#sncosmo.plot_lc(tab,model=fit,errors=res)
#plt.show()
#tab=tab[tab['band']=='F160W']
#print(tab)
#temp=models.FlexibleSpline(tab)

#mod=sncosmo.Model(temp)
#sncosmo.plot_lc(tab,model=mod)


#plt.show()
#temp.ban(np.arange(-10,20,1),sncosmo.get_bandpass('F160W').wave)
sys.exit()
#salt2=sncosmo.SALT2Source(modeldir='/Users/jpierel/rodney/SED_Repository/SEDs.P18-UV2IR/Type_Ia/salt2-extended')
#mod=sncosmo.Model(source=salt2)
mod='snana-2004gv'
#mod=sncosmo.Model(source='snana-2006fo')
#mod=sncosmo.Model(source='snana-2004hx')
#['bessellux','bessellb','bessellv','bessellr','besselli']
snType='Ib'
lcs=sim.createRandMultiplyImagedSN(mod,snType,.1,bands=['bessellb','bessellv','bessellr'],zp=25,cadence=10,epochs=15,numImages=4,timeDelayRange=(0,30),objectName='Test',telescopename='HST',microlensing=False)
for k in lcs.images.keys():
    sncosmo.plot_lc(lcs.images[k].table)
    plt.savefig(k+'.pdf',format='pdf',overwrite=True)
    plt.clf()
sntd.spline_fit(lcs)
lcs.plot_object()
lcs.combine_curves()
#fig=plt.figure()
#ax=fig.gca()
#for img in list(set(lcs.combinedCurve['object'])):
#ax.scatter(lcs.combinedCurve['time'][lcs.combinedCurve['band']==b],lcs.combinedCurve['flux'][lcs.combinedCurve['band']==b])

sncosmo.plot_lc(lcs.combined.table)
plt.savefig('combined.pdf',format='pdf',overwrite=True)
plt.close()


#sys.exit()
#sntd.colorFit(lcs)

lcs=sntd.fit_data(lcs,effect_frames=['rest','free'],effect_names=['host','screen'],snType=snType,constants={'hostr_v':3.1,'screenr_v':3.1},
                  bounds={'screenz':(.4,.6),'z':(.08,.12),'hostebv':(0,1),'screenebv':(-1,1)},dust='CCM89Dust',combined_or_separate='combined',combinedGrids={'td':(-5,5),'mu':(-1,1)})
sncosmo.plot_lc(lcs.combined.table,model=lcs.combined.fit.model,errors=lcs.combined.fit.res)
plt.show()
sys.exit()
lcs.plot_object(filename='type'+snType)

sys.exit()
#lcs.plot_microlensing()

sys.exit()
newTime=np.arange(min(curve['time']),max(curve['time'])+.5,.5)
curve=ml.getDiffCurves(newTime)
fig=plt.figure()
ax=fig.gca()
colors={'u':'b','b':'g','v':'orange','r':'r'}

from scipy.interpolate import splrep,splev
for band in [x for x in curve.colnames if x != 'time']:
    #func=interp1d(curve['time'],curve[band])
    func=splrep(curve['time'],curve[band])
    ax.plot(newTime,splev(newTime,func),color=colors[band])
ax.invert_yaxis()
plt.show()
sys.exit()
'''
files=glob.glob('data/*.dat')

curves=sntd.curveDict(telescopename='Hubble',object='Refsdal')
for f in files:
    curves.add_curve(sntd.read_data(f))

print(curves.images['S3'].table)

#print(timeit.timeit("sntd.fit_data(sntd.read_data('refsdalS2_psfphot.dat'),bounds={'z':(1.2, 1.5)})",setup="import sntd",number=1))

#sntd.write_data(curves,'myData.pkl')
'''
temp=sntd.read_data('myData.pkl')
fig=plt.figure()
ax=fig.gca()
tab=temp.images['S4'].table
#tab.sort('time')
colors=['r','b','g','k','orange','indigo']
i=0
for band in temp.images['S4'].bands:
    temp=tab[tab['band']==band]
    temp.sort('time')
    for j in range(len(temp)-1):
        if temp['time'][j]==temp['time'][j+1]:
            temp['time'][j]-=.001
    print(band,len(temp))
    if len(tab[tab['band']==band])<12:
        continue
    #ax.errorbar(temp['time'],temp['flux'],yerr=temp['fluxerr'],fmt='.',color=colors[i])
    #win=len(tab[tab['band']==band])/3 if len(tab[tab['band']==band])/3%2!=0 else len(tab[tab['band']==band])/3+1
    #ax.plot(tab['time'][tab['band']==band],smooth(tab['flux'][tab['band']==band],win,2),color=colors[i])
    myspline=spl(temp['time'],temp['flux'],w=temp['fluxerr'])
    #ax.plot(temp['time'],myspline(temp['time']),color=colors[i])
    micro=temp['flux']-myspline(temp['time'])
    #ax.errorbar(temp['time'],temp['flux']-micro,yerr=temp['fluxerr'],fmt='.',color=colors[i])
    ax.plot(temp['time'],micro,linestyle='dotted',color=colors[i])
    #maxVal=np.max(myspline(np.arange(min(temp['time']),max(temp['time]']),.5)))
    #ax.plot()
    ax.plot([min(temp['time']),max(temp['time'])],[0,0],linestyle='dashed',color='pink')
    #ax.plot()
    i+=1
#ax.invert_yaxis()
#plt.savefig('data.pdf',format='pdf',overwrite=True)
#plt.savefig('residual.pdf',format='pdf',overwrite=True)
plt.show()

'''
from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np

x=np.linspace(8,12,100)
fig,ax=plt.subplots(1,1)
ax.plot(x,norm.pdf(x,10,.2))
#ax.plot(x,norm.pdf(x,0,.6))
plt.show()
'''