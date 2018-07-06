import sncosmo,sys,math
import numpy as np
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d,splrep,splev,CubicSpline
from sncosmo.utils import integration_grid
from sncosmo.constants import HC_ERG_AA, MODEL_BANDFLUX_SPACING
from collections import Counter
import matplotlib.pyplot as plt
from .util import _findMax

__all__=['SplineSource','BazinSource']


def _param_to_source(source,phase,wave):
    band=None
    for b in np.unique(source.lc['band']):
        temp=sncosmo.get_bandpass(b)
        if temp.minwave()<=min(wave) and temp.maxwave()>=max(wave):
            band=sncosmo.get_bandpass(b)
    if band is None:
        raise RuntimeError("Hmm, your data do not contain the band you want to fit.")
    try:
        zpnorm = 10.**(0.4 * source.lc['zp'][source.lc['band']==band.name.lower()][0])
        band.name=band.name.lower()
    except:
        zpnorm = 10.**(0.4 * source.lc['zp'][source.lc['band']==band.name.upper()][0])
        band.name=band.name.upper()

    #t0=np.mean([np.min(source.lc['time'][source.lc['band']==band.name])-phase[0],np.max(source.lc['time'][source.lc['band']==band.name])-phase[-1]])
    timeInterval=[min(phase),max(phase)]

    finalPhase=np.arange(timeInterval[0],timeInterval[-1],source._tstep)


    wave, dwave = integration_grid(band.minwave(), band.maxwave(),
                                   MODEL_BANDFLUX_SPACING)
    ms = sncosmo.get_magsystem(source.lc['zpsys'][0])
    zpnorm = zpnorm / ms.zpbandflux(band)
    flux=source._param_flux(finalPhase)*HC_ERG_AA/(dwave*np.sum(wave*band(wave))*zpnorm)
    finalWave=np.arange(wave[0]-dwave*10,wave[-1]+dwave*10,dwave)
    finalFlux=np.zeros((len(finalPhase),len(finalWave)))
    for i in range(len(finalPhase)):
        #if finalPhase[i]>= np.min(timeArr) and finalPhase[i] <= np.max(timeArr):
        for j in range(len(finalWave)):
            if finalWave[j]>=wave[0] and finalWave[j]<=wave[-1]:
                finalFlux[i][j]=flux[i]

    if len(flux[np.where(flux==np.max(flux))])==1:
        offset=np.array(finalPhase[flux==np.max(flux)])
    else:
        offset=np.zeros(len(finalPhase))
    out=sncosmo.TimeSeriesSource(np.array(finalPhase)-offset,np.array(finalWave),np.array(finalFlux),zero_before=False)
    return (sncosmo.Model(out))

class BazinSource(sncosmo.Source):

    _param_names=['A','B','fall','rise',]
    param_names_latex=['A','B','t_{fall}','t_{rise}']

    def __init__(self,data,name=None, version=None,tstep=1):
        super(sncosmo.Source, self).__init__() #init for the super class
        self.name = name
        self.version = version
        self._model = {}
        self.lc=_removeDupes(data)
        wave=[]
        for b in np.unique(data['band']):
            wave=np.append(wave,sncosmo.get_bandpass(b).wave)
        wave=np.sort(np.unique(wave))
        wave=np.append([.99*wave[0]],wave)
        wave=np.append(wave,[1.01*wave[-1]])
        self._wave=wave
        self._phase=np.arange(-50,150,1)#np.sort(np.unique(self.lc['time']))

        self._parameters=np.array([1.,0.,15.,5.,])
        self._tstep=tstep


    def _constantBazin(self,length,B):
        return(np.ones(length)*B)

    def _param_flux(self,phase):
        if self._parameters[0]==0:
            bazinFlux=self._constantBazin(len(phase),self._parameters[1])
        else:
            bazinFlux=self._parameters[0]*(np.exp(-phase/self._parameters[2])/(1+np.exp(-phase/self._parameters[3]))) + self._parameters[1]
        if np.inf in bazinFlux or np.any(np.isnan(bazinFlux)):
            return(np.ones(len(phase)))
        return(bazinFlux)



    def _flux(self,phase,wave):
        #if self._parameters[2]<=self._parameters[3]:
        #    return np.ones((len(phase),len(wave)))*(-9999)
        mod=_param_to_source(self,phase,wave)
        return(mod._flux(phase,wave))



class NewlingSource(sncosmo.Source):
    _param_names=['A','psi','sigma','k','phi']
    param_names_latex=['A','\psi','\sigma','k','phi']

    def __init__(self,data,name=None, version=None,tstep=1):
        super(sncosmo.Source, self).__init__() #init for the super class
        self.name = name
        self.version = version
        self._model = {}
        self.lc=_removeDupes(data)
        wave=[]
        for b in np.unique(data['band']):
            wave=np.append(wave,sncosmo.get_bandpass(b).wave)
        wave=np.sort(np.unique(wave))
        wave=np.append([.99*wave[0]],wave)
        wave=np.append(wave,[1.01*wave[-1]])
        self._wave=wave
        self._phase=np.arange(-200,500,5)

        self._parameters=np.array([1.,0.,1.,1.,-1.])
        self._tstep=tstep



    def _param_flux(self,phase):
        #self._parameters[4]=np.min([np.min(phase),self._parameters[4]])
        #self._parameters[4]=-self._parameters[3]*self._parameters[2]

        splPhase=phase[phase>=self._parameters[4]]
        splPhase=splPhase[splPhase<=(self._parameters[3]*self._parameters[2]+self._parameters[4])]

        spl=CubicSpline([self._parameters[4],0.],[0,self._parameters[1]])

        Psi=np.zeros(len(phase))
        for i in range(len(Psi)):
            if phase[i] in splPhase:
                Psi[i]=spl(phase[i])
            elif phase[i]>=(self._parameters[3]*self._parameters[2]+self._parameters[4]):
                Psi[i]=self._parameters[1]
        #print(Psi[0])
        #Psi[phase==splPhase]=spl(splPhase)
        #print(((phase+self._parameters[4])/self._parameters[2]),((phase+self._parameters[4])/self._parameters[2])**self._parameters[3],self._parameters,np.any(np.isnan(phase)))

        newlingFlux=(self._parameters[0]*((phase-self._parameters[4])/self._parameters[2])**self._parameters[3])*np.exp(-(phase-self._parameters[4])/self._parameters[2])*(self._parameters[3]**(-self._parameters[3]))*(np.exp(self._parameters[3]))+Psi
        #print(phase[0],self._parameters[4])
        if np.inf in newlingFlux or np.any(np.isnan(newlingFlux)):
                #print(self._parameters)
                return(np.zeros(len(phase)))
        return(newlingFlux)

    def _flux(self,phase,wave):
        #if self._parameters[2]<=self._parameters[3]:
        #    return np.ones((len(phase),len(wave)))*(-9999)
        mod=_param_to_source(self,phase,wave)
        return(mod._flux(phase,wave))

class KarpenkaSource(sncosmo.Source):
    _param_names=['A','B','t1','rise','fall']
    param_names_latex=['A','B','t_1','t_{rise}','t_{fall}']

    def __init__(self,data,name=None, version=None,tstep=1):
        super(sncosmo.Source, self).__init__() #init for the super class
        self.name = name
        self.version = version
        self._model = {}
        self.lc=_removeDupes(data)
        wave=[]
        for b in np.unique(data['band']):
            wave=np.append(wave,sncosmo.get_bandpass(b).wave)
        wave=np.sort(np.unique(wave))
        wave=np.append([.99*wave[0]],wave)
        wave=np.append(wave,[1.01*wave[-1]])
        self._wave=wave
        self._phase=np.arange(-50,150,1)

        self._parameters=np.array([1.,0.,0.,1.,1.,])
        self._tstep=tstep



    def _param_flux(self,phase):
        karpenkaFlux=(self._parameters[0]*(1+self._parameters[1]*((phase+self._parameters[2])**2)))*(np.exp(-phase/self._parameters[4])/(1+np.exp(-phase/self._parameters[3])))
        if np.inf in karpenkaFlux or np.any(np.isnan(karpenkaFlux)):
            return(np.zeros(len(phase)))
        return(karpenkaFlux)

    def _flux(self,phase,wave):
        #if self._parameters[2]<=self._parameters[3]:
        #    return np.ones((len(phase),len(wave)))*(-9999)
        mod=_param_to_source(self,phase,wave)
        return(mod._flux(phase,wave))

class SplineSource(sncosmo.Source):


    #_param_names = ['dt0','amplitude']

    #param_names_latex=['\Delta \ t_0','A']
    def __init__(self, data,weights=None,name=None, version=None,tstep=1,wstep=10,knots=3,degree=3,func='spline'):
        super(sncosmo.Source, self).__init__() #init for the super class
        self.name = name
        self.version = version
        self._model = {}
        data=_removeDupes(data)
        self.lc=data


        self._phase=np.sort(np.unique(self.lc['time']))
        self._bands=np.unique(self.lc['band'])
        #self._param_names=np.append(['dt0_'+str(i) for i in range(len(np.unique(self.lc['band'])))],['amplitude_'+str(i) for i in range(len(np.unique(self.lc['band'])))])
        #self._parameters = np.append([0. for i in range(len(np.unique(self.lc['band'])))],[1. for i in range(len(np.unique(self.lc['band'])))])#,_findMax(time,flux)[0]])
        #self.param_names_latex=np.append(['\Delta \ t_0 \ ('+b+')' for b in np.unique(self.lc['band'])],['A \ ('+b+')' for b in np.unique(self.lc['band'])])
        self._param_names=['dt0','amplitude']
        self._parameters=np.array([0.,1.])
        self.t0=dict([])
        for b in self._bands:
            t0Est,t0Flux=_findMax(self.lc['time'][self.lc['band']==b],self.lc['flux'][self.lc['band']==b])
            t0Est=t0Est[0] if isinstance(t0Est,np.ndarray) else t0Est
            t0Flux=t0Flux[0] if isinstance(t0Flux,np.ndarray) else t0Flux
            self.t0[str(b)]=[t0Est,t0Flux]

        self._func=func

        self.param_names_latex=['\Delta \ t_0','A']
        self._knots=knots
        self._deg=degree
        wave=[]
        for b in np.unique(data['band']):
            if len(data[data['band']==b])<knots:
                data=data[data['band']!=b]
                continue
            wave=np.append(wave,sncosmo.get_bandpass(b).wave)
        wave=np.sort(np.unique(wave))
        wave=np.append([.99*wave[0]],wave)
        wave=np.append(wave,[1.01*wave[-1]])



        self._wave=wave

        self._steps={'tstep':tstep,'wstep':wstep}
        #self._time = np.linspace(np.min)
        #self._wave = wave
        if not weights:
            self._weights=None
        #    self._weights=np.ones(len(data))


        #if self._parameters[1] not in self._time:
        #    self._t0Set()


    def _dataToSource(self,phase,waveArr):
        finalWave,finalPhase,finalFlux=None,None,[]
        bands=[b for b in np.unique(self.lc['band']) if waveArr[0] >= sncosmo.get_bandpass(b).minwave() and waveArr[-1] <= sncosmo.get_bandpass(b).maxwave()]
        #for b in np.unique(self.lc['band']):
        #    print(sncosmo.get_bandpass(b).minwave(),waveArr[0],sncosmo.get_bandpass(b).maxwave(),waveArr[-1])
        for b in bands:
            #print(type(self._t0),type(self._t0[b]))
            #print(self._t0[b],'test')
            #tup=self._t0[b]#_findMax(self.lc['time'][self.lc['band']==b],self.lc['flux'][self.lc['band']==b])

            t0Est=self.t0[b][0]
            t0Flux=self.t0[b][1]
            timeInterval=[float(np.min(self.lc['time'][self.lc['band']==b])),float(np.max(self.lc['time'][self.lc['band']==b]))]
            phase=phase[phase>=timeInterval[0]]
            phase=phase[phase<=timeInterval[-1]]
            tempTime=np.array(self.lc['time'][self.lc['band']==b])
            tempFlux=np.array(self.lc['flux'][self.lc['band']==b])
            #if b=='F160W':
            #    print(tempTime,tempFlux)
            if self._weights is not None:
                weights=np.ones(len(tempFlux))
            else:
                weights=np.array(self.lc['fluxerr'][self.lc['band']==b])

            #ind=np.where(self._bands==b)[0]
            if t0Est is not None and t0Est not in self.lc['time'][self.lc['band']==b]:
                self._t0=t0Est
                tempTime=np.append(tempTime,self._parameters[0]+np.array([t0Est]))
                tempFlux=np.append(tempFlux,self._parameters[1]*np.array([t0Flux]))
                #tempTime=np.append(tempTime,self._parameters[ind]+np.array([t0Est]))
                #tempFlux=np.append(tempFlux,self._parameters[-1]*np.array([t0Flux]))
                weights=np.append(weights,np.min(weights))
            else:
                tempTime=np.append(tempTime,self._parameters[0]+tempTime[np.where(tempFlux==np.max(tempFlux))])
                tempFlux=np.append(tempFlux,self._parameters[1]*np.max(tempFlux))
                #tempTime=np.append(tempTime,self._parameters[ind]+tempTime[np.where(tempFlux==np.max(tempFlux))])
                #tempFlux=np.append(tempFlux,self._parameters[-1]*np.max(tempFlux))
                weights=np.append(weights,weights[np.where(tempFlux==np.max(tempFlux))])

            srt=np.argsort(tempTime)
            tempTime=tempTime[srt]
            tempFlux=tempFlux[srt]
            weights=weights[srt]
            #if b=='F105W':
            #    print(tempTime)


            band=sncosmo.get_bandpass(b)
            timeArr=np.arange(timeInterval[0],timeInterval[1]+self._steps['tstep'],self._steps['tstep'])
            wave, dwave = integration_grid(band.minwave(), band.maxwave(),
                                           MODEL_BANDFLUX_SPACING)
            #timeArr-=t0Est
            #print(timeArr)
            #print(phase[0],timeInterval[0],phase[-1],timeInterval[-1])
            if phase[0] <timeInterval[0] or phase[-1] > timeInterval[-1]:
                raise RuntimeError('The phase you requested is outside your data bounds.')
            zpnorm = 10.**(0.4 * self.lc['zp'][self.lc['band']==b][0])


            ms = sncosmo.get_magsystem(self.lc['zpsys'][0])
            zpnorm = zpnorm / ms.zpbandflux(b)

            if self._func=='spline':
                spl=splrep(tempTime,tempFlux,k=int(self._knots),w=1./weights,s=len(tempFlux)/2)
                flux=splev(timeArr,spl)*HC_ERG_AA/(dwave*np.sum(wave*band(wave))*zpnorm)
            elif self._func=='chebyshev':
                cheb=np.polynomial.chebyshev.chebfit(tempTime,tempFlux,deg=int(self._deg),w=1./weights)
                flux=np.polynomial.chebyshev.chebval(timeArr,cheb)*HC_ERG_AA/(dwave*np.sum(wave*band(wave))*zpnorm)
            else:
                raise RuntimeError("Don't know that type of flexible function.")


            if finalWave is None:
                finalWave=np.arange(wave[0]-dwave*10,wave[-1]+dwave*10,dwave)
                finalPhase=timeArr
                finalFlux=np.zeros((len(finalPhase),len(finalWave)))
                for i in range(len(finalPhase)):
                    #if finalPhase[i]>= np.min(timeArr) and finalPhase[i] <= np.max(timeArr):
                    for j in range(len(finalWave)):
                        if finalWave[i]>=wave[0] and finalWave[i]<=wave[-1]:
                            finalFlux[i][j]=flux[i]



            source=sncosmo.TimeSeriesSource(np.array(finalPhase),np.array(finalWave),np.array(finalFlux),zero_before=False)


            mod=sncosmo.Model(source)

            #mod.set(t0=t0Est)
            #print(mod.bandflux(b,timeArr))
            #print(self.lc[self.lc['band']==b])
            #mod.set(t0=t0Est)




            #trans = band(wave)
            #print(np.sum(wave*trans),HC_ERG_AA)


            #sncosmo.plot_lc(data=self.lc[self.lc['band']==b],model=mod)
            #plt.show()
            #sys.exit()
        return (mod)



    def _flux(self, phase, wave):
        model=self._dataToSource(phase,wave)

        return(model._flux(phase,wave))




def _removeDupes(data):
    tempTable=Table(names=data.colnames,dtype=[data.dtype[x] for x in data.colnames])

    for b in np.unique(data['band']):
        dupes=[item for item, count in Counter(data['time'][data['band']==b]).items() if count > 1]
        duped=[]
        temp=data[data['band']==b]

        for row in temp:
            if row['time'] not in dupes:
                tempTable.add_row(row)
            else:
                if row['time'] in duped:
                    continue

                row['flux']=np.average(temp['flux'][temp['time']==row['time']],
                                       weights=1./(temp['fluxerr'][temp['time']==row['time']])**2)
                row['fluxerr']=np.sqrt(1./np.sum(1./(temp['fluxerr'][temp['time']==row['time']])**2))
                duped.append(row['time'])
                tempTable.add_row(row)


    return (tempTable)