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

def createRandMultiplyImagedSN(sourcename,snType,redshift,telescopename='telescope',objectName='object',timeDelayRange=(0,30),muRange=(1,20),
                               numImages=2,cadence=5,epochs=50,bands=['F140W','F105W','F160W'],gain=1000.,skynoiseRange=(1,5),mjdRange=None,zpsys='ab',zp=None,microlensing=False,
                               cc_av=.9,ia_av=.3,dust='CCM89Dust',screenz=None,screenebv=None,effect_frames=['rest','free'],effect_names=['host','screen'],minsnr=0.0,av_host=.3,z_lens=None, av_lens=None,scatter=True):

    #TODO sample from dust prior for host and lens-plane dust.
    obj=curveDict(telescopename=telescopename,object=objectName)
    if not mjdRange:
        now=np.round(Time(datetime.datetime.now()).mjd,3)
        times=np.linspace(now,now+cadence*epochs,epochs)
    else:
        times=np.linspace(mjdRange[0],mjdRange[-1],epochs)
    bandList=np.array([np.tile(b,len(times)) for b in bands]).flatten()
    ms=sncosmo.get_magsystem(zpsys)

    zpList=[ms.band_flux_to_mag(1,b) for b in bandList] if not zp else [zp for i in range(len(bandList))]
    obj.zpDict = dict([(bandList[i], zpList[i]) for i in range(len(bandList))])
    obj.zpsys = zpsys
    obs=Table({'time':np.tile(times,len(bands)),'band':bandList,'zpsys':[zpsys for i in range(len(bandList))],'zp':zpList,'skynoise':np.random.uniform(skynoiseRange[0],skynoiseRange[1],len(bandList)),'gain':[gain for i in range(len(bandList))]})

    absolutes=_getAbsoluteDist()

    # Set up the dust extinction effects in the host galaxy and lens plane
    # TODO : allow different lens-plane dust for each image?
    R_V = 3.1  # TODO: allow alternate dust law
    if dust:
        dust_effect = {'SFD98Map': sncosmo.SFD98Map,
                       'CCM89Dust': sncosmo.CCM89Dust,
                       'OD94Dust': sncosmo.OD94Dust,
                       'F99Dust': sncosmo.F99Dust}[dust]()
        if av_lens:
            effect_frames = ['rest', 'free']
            effect_names = ['host', 'lens']
        else:
            effect_frames = ['rest']
            effect_names = ['host']
    else:
        dust_effect=[]
        effect_names = []
        effect_frames = []

    effects=[dust_effect for i in range(len(effect_names))]
    effect_names=effect_names if effect_names else []
    effect_frames=effect_frames if effect_frames else []
    if not isinstance(effect_names,(list,tuple)):
        effects=[effect_names]
    if not isinstance(effect_frames,(list,tuple)):
        effects=[effect_frames]

    model=sncosmo.Model(source=sourcename, effects=effects,
                        effect_names=effect_names, effect_frames=effect_frames)
    model.set(z=redshift)

    if z_lens is None:
        z_lens=redshift/2.
    if av_lens:# is None:
        # av_lens=np.abs(np.random.normal(0,1))
        ebv_lens = av_lens/R_V
        model.set(lensz=z_lens, lensebv=ebv_lens)

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
    t0 = times[0] + .25 * (times[-1] - times[0])

    if snType=='Ia':
        x0=model.get('x0')
        params={'z':redshift,'t0':t0,'x0':x0,'x1':np.random.normal(0.,1.),
                'c':np.random.normal(0.,.1),'hostebv':av_host/R_V}
    else:
        amp=model.get('amplitude')
        params={'z':redshift,'t0':t0,'amplitude':amp,'hostebv':av_host/R_V}
    params['hostr_v']=R_V
    if av_lens:
        params['lensr_v']=R_V


    for i in range(numImages):
        model=sncosmo.Model(source=sourcename,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
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
        # Draw a random magnification and time delay separately for each image
        mu = np.random.uniform(muRange[0], muRange[-1])
        delay = np.random.uniform(timeDelayRange[0], timeDelayRange[-1])

        # Make a separate model_i for each SN image, so that the magnification
        # can be reflected in the model_i parameters and propagate correctly into
        # flux uncertainties, and the observation epochs will be the same for
        # all images.
        model_i = deepcopy(model)
        params_i = deepcopy(params)
        if snType=='Ia':
            params_i['x0'] *= mu
        else:
            params_i['amplitude'] *= mu
        params_i['t0'] += delay
        lc_i = sncosmo.realize_lcs(obs, model_i, [params_i],
                                   trim_observations=True, scatter=scatter)[0]
        tempCurve=curve(zp=zp, zpsys=zpsys)

        # TODO : implement a more general microlensing apparatus
        if microlensing:
            if lc_i['band'][0].lower().startswith('bessell'):
                bandset = 'bessell'
            elif lc_i['band'][0].lower().startswith('f1'):
                bandset = 'hst'
            mlCurve=getDiffCurve(
                lc_i['time'][lc_i['band']==lc_i['band'][0]],
                bandset=bandset)
            lc_i=_addML(mlCurve,lc_i)
            tempCurve.ml=mlCurve
        else:
            tempCurve.ml=None
        tempCurve.table=lc_i
        tempCurve.bands=list(set(lc_i['band']))
        tempCurve.simMeta=deepcopy(lc_i.meta)
        tempCurve.simMeta['av_lens']=av_lens
        tempCurve.simMeta['z_lens']=z_lens
        tempCurve.simMeta['mu']=mu
        tempCurve.simMeta['td']=delay

        model_i.set(**params_i)
        tempCurve.simMeta['model']=model_i

        obj.add_curve(tempCurve)

    # Store the un-lensed model as a component of the lensed SN object.
    model.set(**params)
    obj.model = model
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