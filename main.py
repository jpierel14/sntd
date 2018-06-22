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
import sntd,timeit,glob,sys,warnings,sncosmo
from multiprocessing import Pool
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter as smooth
from scipy.interpolate import UnivariateSpline as spl
import numpy as np
from scipy.interpolate import interp1d

from sntd import ml,simulation as sim,plotting
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