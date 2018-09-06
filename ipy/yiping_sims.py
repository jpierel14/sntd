
# coding: utf-8

# In[2]:

import glob,os,sys,sncosmo,pyParz
import sncosmo
import sntd
from sntd import simulation, fitting,io,models,util
from sntd.plotting import _COLORLIST5
import numpy as np
from scipy import stats
from copy import copy,deepcopy
from astropy.table import vstack
from astropy.io import ascii
from matplotlib import pyplot as plt
from astropy.table import Table


# In[3]:

times=[-3.932323367086956978e+01,-3.832323367086956978e+01,
       -2.832323367086956978e+01,-1.832323367086956978e+01,
       -8.323233670869569778e+00,1.676766329130430222e+00,
       1.167676632913043022e+01,2.167676632913043022e+01,
       3.167676632913043022e+01,4.167676632913043022e+01,
       5.167676632913043022e+01,6.167676632913043022e+01,
       7.167676632913043022e+01,8.167676632913043022e+01,
       9.167676632913043022e+01,1.016767663291304302e+02,
       1.116767663291304302e+02,1.216767663291304302e+02,
       1.316767663291304302e+02,1.416767663291304302e+02]


# In[6]:

#TODO: Make this parallelization a function inside SNTD
#reload(simulation)
nsim = 1000
zp=27
zpsys='AB'
dt_fit_list = [] 
murel_fit_list = [] 
inds={'time':0,'flux':1,'fluxerr':2,'object':3}
default={'zp':zp,'zpsys':'AB','band':'F125W'}

def _td_par(args):  
    lcs=io.curveDict()
    for t in args:
        temp=Table()
        for name in ['time','band','flux','fluxerr','zp','zpsys']:
            if name in ['time','flux','fluxerr']:
                temp[name]=t[inds[name]]
                
            else:
                
                temp[name]=default[name]
            
        curve_i=io.curve(zp=zp,zpsys=zpsys)
        curve_i.table=temp
        
        curve_i.object=t[inds['object']]
        
        curve_i.bands=['F125W']
        lcs.add_curve(curve_i)

    # Part 2: fit each light curve separately to determine lensing parameters
    print("Fitting strongly lensed SN %i"%args[0][-1])
    lcs_tdfit=fitting.fit_data(lcs, snType='IIP', models=['snana-2006iw'],
                               dust='CCM89Dust',effect_frames=['rest'],
                                effect_names=['host'],
                               params=['z','amplitude','t0','hostebv'],
                               bounds={'hostebv':(0,.5),'z':(1.26,1.3),'t0':(-15,15)},
                               combined_or_separate='separate',showPlots=False)
    


    return(np.abs(lcs_tdfit.time_delays['S2']))

simmed=[]
#Currently in order to make unique SNe they need to be created in serial 
#then fit in parallel, I'll look into this
for isim in range(nsim):
    print("Simulating...%i"%isim)
    modname = 'snana-2006iw'
    snType = 'IIP'
    bandlist = ['F110W']
    lcs = simulation.createMultiplyImagedSN(
        modname, snType, 1.28, bands=bandlist,
        zp=zp, skynoiseRange=[0.001,0.005],gain=50.,cadence=1., epochs=10.,
        timeArr=times, time_delays=[0., 56.],
        magnifications=[2.6,1.3], objectName='Test'+snType, 
        telescopename='HST',z_lens=0.51, 
        microlensing_type=None,minsnr=5.0)
    #print(np.max(np.array(lcs.images['S1'].table['flux'])))
    simmed.append([[np.array(lcs.images[k].table['time']),
            np.array(lcs.images[k].table['flux']),
            np.array(lcs.images[k].table['fluxerr']),k,isim] for k in lcs.images.keys()])

print("")
dt_fit_list=pyParz.foreach(simmed,_td_par,None)

np.savetxt('tds_'+snType+'.dat',dt_fit_list) #uncomment to output the histogram data


#plot the resulting histogram
fig=plt.figure()
ax=fig.gca() 
ax.hist(dt_fit_list,rwidth=1,normed=True,label='Measured Time Delay')

ax.plot([56,56],ax.get_ylim(),'r--',label='Actual Time Delay') 
ax.legend(loc='upper left',fontsize=10) 
ax.set_title('Type Ib Doubly-Imaged Time Delay Measurement--N='+str(nsim),size=12) 
ax.set_ylabel('Probability Density',size=12) 
ax.set_xlabel('Relative Time Delay (Days)',size=12) 
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)

plt.savefig('tdhist_'+snType+'.pdf',format='pdf',overwrite=True)
plt.close()


# In[ ]:



