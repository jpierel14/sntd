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
microlensing
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
salt2=sncosmo.SALT2Source(modeldir='/Users/jpierel/rodney/SED_Repository/SEDs.P18-UV2IR/Type_Ia/salt2-extended')
mod=sncosmo.Model(source=salt2)
#mod=sncosmo.Model(source='snana-2004gv')
#mod=sncosmo.Model(source='snana-2006fo')
#mod=sncosmo.Model(source='snana-2004hx')

snType='Ia'
lcs=sim.createRandMultiplyImagedSN(mod,snType,.1,bands=['bessellb','bessellv','bessellr'],zp=100,cadence=2,epochs=10,numImages=4,objectName='Test',telescopename='HST',microlensing=True)
#print(lcs.table)
#print(lcs.images['S1'].table)
#sntd.colorFit(lcs)
sntd.fit_data(lcs,snType=snType,bounds={'z':(.05,.15)})

lcs.plot_object(filename='type'+snType)

