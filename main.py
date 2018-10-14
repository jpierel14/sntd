import glob,os,sys,sncosmo
from sntd import models,fitting
import sntd
from copy import copy
from astropy.table import vstack
from astropy.io import ascii
import numpy as np
#define which photometry we want to analyze (dolphot or hstphot)


photometry='dolphot'

if photometry=='dolphot':
    files=glob.glob(os.path.join('data','pat_refsdal',photometry,'stringent','*.dat'))
elif photometry=='hstphot':
    files=glob.glob(os.path.join('data','pat_refsdal',photometry,'*.dat'))
else:
    raise RuntimeError("Photometry not recognized.")

#Read in the data from pat
bands=[]
lcs=sntd.curveDict(telescopename='Hubble',object='Refsdal_'+photometry)
images=['s1','s2','s3','s4']
for im in images:
    imTab=None
    for f in [x for x in files if im in x]:
        tab=ascii.read(f,comment='=')
        tab.remove_column('flag')
        tab.rename_column('mhjd','time')
        tab=tab[tab['time']>tab['time'][tab['flux']==np.max(tab['flux'])]-200]
        band=f[f.find('_'+photometry)-5:f.find('_'+photometry)]
        bands.append(band)
        tab['band']=band
        tab['zpsys']='ab'
        tab['zp']=25.0
        if imTab is None:
            imTab=copy(tab)
        else:
            imTab=vstack([imTab,tab])

    curve=sntd.table_factory(imTab,telescopename='Hubble',object=os.path.basename(f)[5:7].upper())
    lcs.add_curve(curve)

#Fit with splines first
guess={'F160W':{'S1':57145,'S2':57165,'S3':57114,'S4':57165},'F125W':{'S1':57142,'S2':57160,'S3':57098,'S4':57160}}
for b in lcs.bands:
    print(b)
    if b=='F160W':
        bds={'t0':(-5,5),'amplitude':(.75,1.25),'fall':(50,130),'rise':(20,55),'B':(0,.2)}
    else:
        bds={'t0':(-5,5),'amplitude':(.75,1.25),'fall':(60,180),'rise':(20,55),'B':(0,.2)}
    curves=fitting.fit_data(lcs,bands=[b],snType='II',models=['BazinSource'],params=['amplitude','B','fall','rise','t0'],flip=False,t0_guess=guess[b],bounds=bds,combined_or_separate='separate')
    #delays=fitting.timeDelaysAndMagnifications(curves)
    print(curves.time_delays)
sys.exit()

guess={'F160W':{'S1':57164,'S2':57163,'S3':57115,'S4':57153},'F125W':{'S1':57135,'S2':57157,'S3':57088,'S4':57150}}
for b in lcs.bands:
    if b =='F160W':
        bds={'t0':(-10,10),'amplitude':(.5,1.15),'k':(.5,1.1),'sigma':(80,100),'s':(-.05,.2)}
    else:
        bds={'t0':(-10,10),'amplitude':(.85,1.15),'k':(1,2),'sigma':(70,90),'s':(-.05,.2)}
    curves=fitting.fit_data(lcs,bands=[b],snType='II',models=['PierelSource'],params=['k','sigma','amplitude','t0','s'],flip=False,t0_guess=guess[b],bounds=bds,combined_or_separate='separate')
    delays=fitting.timeDelays(curves)