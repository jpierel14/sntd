import sncosmo,datetime,sys,os,scipy
import numpy as np
from astropy.time import Time
from astropy.table import Table
from astropy.io import ascii
from copy import deepcopy

from .util import __dir__
from .io import curve,curveDict
from .ml import getDiffCurve

__all__=['createRandMultiplyImagedSN']
def _getAbsoluteDist():
    absolutes=ascii.read(os.path.join(__dir__,'sim','data','absolutes.ref'))
    total=float(np.sum(absolutes['N'][absolutes['type']!='Ia']))
    absDict=dict([])
    for row in absolutes:
        if row['type']=='Ia':
            frac=1
        else:
            frac=float(row['N'])/total
        absDict[row['type']]={'dist':(row['mean'],row['sigma']),'frac':frac}
    return(absDict)

def _getAbsFromDist(dist):
    mu,sigma=dist
    return(np.random.normal(mu,sigma))

def createRandMultiplyImagedSN(source,snType,redshift,telescopename='telescope',objectName='object',timeDelayRange=(0,30),muRange=(1,20),
                               numImages=2,cadence=5,epochs=50,bands=['F140W','F105W','F160W'],gain=1000.,skynoiseRange=(1,5),mjdRange=None,zpsys='ab',zp=None,microlensing=False,
                               cc_av=.9,ia_av=.3,dust='CCM89Dust',screenz=None,screenebv=None,effect_frames=['rest','free'],effect_names=['host','screen'],minsnr=0.0):
    #TODO sample from dust prior and add intermediate redshift dust https://github.com/sncosmo/sncosmo/pull/189
    obj=curveDict(telescopename=telescopename,object=objectName)
    if not mjdRange:
        now=np.round(Time(datetime.datetime.now()).mjd,3)
        times=np.linspace(now,now+cadence*epochs,epochs)
    else:
        times=np.linspace(mjdRange[0],mjdRange[-1],epochs)
    bandList=np.array([np.tile(b,len(times)) for b in bands]).flatten()
    ms=sncosmo.get_magsystem(zpsys)

    zpList=[ms.band_flux_to_mag(1,b) for b in bandList] if not zp else [zp for i in range(len(bandList))]
    obs=Table({'time':np.tile(times,len(bands)),'band':bandList,'zpsys':[zpsys for i in range(len(bandList))],'zp':zpList,'skynoise':np.random.uniform(skynoiseRange[0],skynoiseRange[1],len(bandList)),'gain':[gain for i in range(len(bandList))]})

    absolutes=_getAbsoluteDist()
    if dust:
        dust={'SFD98Map':sncosmo.SFD98Map,'CCM89Dust':sncosmo.CCM89Dust,'OD94Dust':sncosmo.OD94Dust,'F99Dust':sncosmo.F99Dust}[dust]()

    else:
        dust=[]
    effects=[dust for i in range(len(effect_names))] if effect_names else []
    effect_names=effect_names if effect_names else []
    effect_frames=effect_frames if effect_frames else []
    if not isinstance(effect_names,(list,tuple)):
        effects=[effect_names]
    if not isinstance(effect_frames,(list,tuple)):
        effects=[effect_frames]
    model=sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)

    if not screenz:
        screenz=redshift/2.
    if not screenebv:
        screenebv=np.random.uniform(-1,1,1)[0]
    model.set(z=redshift,screenz=screenz,screenebv=screenebv)

    if snType in ['IIP','IIL','IIn']:
        absBand='bessellb'
    else:
        absBand='bessellr'
    model.set_source_peakabsmag(_getAbsFromDist(absolutes[snType]['dist']),absBand,zpsys)
    #print(model.source_peakabsmag('bessellr',zpsys))

    #model.set(t0=times[0]+.25*(times[-1]-times[0]))

    if snType=='Ia':
        x0=model.get('x0')

        params={'z':redshift,'t0':times[0]+.3*(times[-1]-times[0]),'x0':x0,'x1':np.random.normal(0.,1.),'c':np.random.normal(0.,.1),'hostebv':ia_av}
    else:
        amp=model.get('amplitude')
        params={'z':redshift,'t0':times[0]+.3*(times[-1]-times[0]),'amplitude':amp,'hostebv':cc_av/3.1}
    params['hostr_v']=3.1
    params['screenr_v']=3.1


    #model.set(**params)
    #print(model.bandflux('bessellb',params['t0'],zp=25,zpsys='ab'))
    #sys.exit()

    #lc=sncosmo.photdata.photometric_data(sncosmo.realize_lcs(obs,model,[params],trim_observations=True)[0])
    #lc=lc.normalized()
    #lc=Table([lc.time,[x.name for x in lc.band],lc.flux,lc.fluxerr,lc.zp,lc.zpsys],names=['time','band','flux','fluxerr','zp','zpsys'])


    for i in range(numImages):
        model=sncosmo.Model(source=source,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
        tempParams=deepcopy(params)
        delay=np.random.uniform(timeDelayRange[0],timeDelayRange[-1])
        mu=np.random.uniform(muRange[0],muRange[-1])
        tempParams['t0']+=delay
        if snType=='Ia':
            tempParams['x0']*=mu
        else:
            tempParams['amplitude']*=mu

        lc=sncosmo.realize_lcs(obs,model,[tempParams],trim_observations=True,scatter=False)[0]
        #do my own scatter, snnosmo's is weird
        for i in range(len(lc)):
            temp=np.random.normal(lc['flux'][i],lc['fluxerr'][i])
            while np.abs(temp)>(np.abs(lc['flux'][i])+np.abs(lc['fluxerr'][i])):
                temp=np.random.normal(lc['flux'][i],lc['fluxerr'][i])
            lc['flux'][i]=temp
        lc=lc[(np.abs(lc['flux']/lc['fluxerr']))>=minsnr]


        #print(lc)
        #temp=deepcopy(lc)

        #temp['time']+=delay
        if microlensing:

            lc,mlCurves=_addML(lc)

        #temp['flux']*=mu
        #temp['fluxerr']*=mu


        tempCurve=curve(zp=zp,zpsys=zpsys)

        tempCurve.table=lc
        tempCurve.bands=list(set(lc['band']))
        tempCurve.simMeta=lc.meta
        tempCurve.simMeta['screenebv']=screenebv
        tempCurve.simMeta['screenz']=screenz
        tempCurve.simMeta['mu']=mu
        tempCurve.simMeta['td']=delay
        if microlensing:
            tempCurve.ml=mlCurves
        else:
            tempCurve.ml=None
        obj.add_curve(tempCurve)




    return(obj)

def _addML(myCurve):
    mlCurves=dict([])
    num=np.random.randint(1,5)
    print(num)
    for b in set(myCurve['band']):
        mlCurve=getDiffCurve(myCurve['time'][myCurve['band']==b],num)
        myCurve['flux'][myCurve['band']==b]/=mlCurve[b]
        myCurve['fluxerr'][myCurve['band']==b]/=mlCurve[b]
        mlCurves[b]=mlCurve[b]

    return(myCurve,mlCurves)