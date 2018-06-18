import sncosmo,datetime,sys,os
import numpy as np
from astropy.time import Time
from astropy.table import Table
from astropy.io import ascii
from copy import deepcopy
from .util import __dir__
from .io import curve,curveDict
from .ml import getDiffCurve
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

def createRandMultiplyImagedSN(model,snType,redshift,telescopename='telescope',objectName='object',timeDelayRange=(0,100),muRange=(1,20),numImages=2,cadence=5,epochs=50,bands=['F140W','F105W','F160W'],gain=1.,skynoise=0.,mjdRange=None,zpsys='ab',zp=None,microlensing=False):
    obj=curveDict(telescopename=telescopename,object=objectName)
    if not mjdRange:
        now=np.round(Time(datetime.datetime.now()).mjd,3)
        times=np.linspace(now,now+cadence*epochs,epochs)
    else:
        times=np.linspace(mjdRange[0],mjdRange[-1],epochs)
    bandList=np.array([np.tile(b,len(times)) for b in bands]).flatten()
    ms=sncosmo.get_magsystem(zpsys)

    zpList=[ms.band_flux_to_mag(1,b) for b in bandList] if not zp else [zp for i in range(len(bandList))]
    obs=Table({'time':np.tile(times,len(bands)),'band':bandList,'zpsys':[zpsys for i in range(len(bandList))],'zp':zpList,'skynoise':np.random.uniform(80,150,len(bandList)),'gain':[gain for i in range(len(bandList))]})


    absolutes=_getAbsoluteDist()

    model.set(z=redshift)
    if snType in ['IIP','IIL','IIn']:
        absBand='bessellb'
    else:
        absBand='bessellr'
    model.set_source_peakabsmag(_getAbsFromDist(absolutes[snType]['dist']),absBand,zpsys)
    #print(model.source_peakabsmag('bessellr',zpsys))
    model.set(t0=times[0]+.25*(times[-1]-times[0]))
    if snType=='Ia':
        x0=model.get('x0')

        params={'z':redshift,'t0':times[0]+.5*(times[-1]-times[0]),'x0':x0,'x1':np.random.normal(0.,1.),'c':np.random.normal(0.,.1)}
    else:
        amp=model.get('amplitude')
        params={'z':redshift,'t0':times[0]+.5*(times[-1]-times[0]),'amplitude':amp}
    lc=sncosmo.realize_lcs(obs,model,[params],trim_observations=True)[0]
    for i in range(numImages):
        temp=deepcopy(lc)
        delay=np.random.uniform(timeDelayRange[0],timeDelayRange[-1])
        temp['time']+=delay
        if microlensing:
            mlCurve=getDiffCurve(temp['time'][temp['band']==temp['band'][0]])
            temp=_addML(mlCurve,temp)

        mu=np.random.uniform(muRange[0],muRange[-1])
        temp['flux']*=mu
        temp['fluxerr']*=mu


        tempCurve=curve(zp=zp,zpsys=zpsys)

        tempCurve.table=temp
        tempCurve.bands=list(set(temp['band']))
        tempCurve.simMeta=deepcopy(lc.meta)

        tempCurve.simMeta['mu']=mu
        tempCurve.simMeta['td']=delay
        if microlensing:
            tempCurve.ml=mlCurve
        else:
            tempCurve.ml=None
        obj.add_curve(tempCurve)




    return(obj)

def _addML(mlCurve,myCurve):
    for b in set(myCurve['band']):
        myCurve['flux'][myCurve['band']==b]/=mlCurve[b]
        myCurve['fluxerr'][myCurve['band']==b]/=mlCurve[b]
    return(myCurve)