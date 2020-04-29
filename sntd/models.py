import sncosmo
import numpy as np
from astropy.table import Table
from scipy.interpolate import CubicSpline,interp1d
from sncosmo.utils import integration_grid
from sncosmo.constants import HC_ERG_AA, MODEL_BANDFLUX_SPACING
from scipy.stats import exponnorm
from collections import Counter


__all__=['BazinSource']


def _param_to_source(source,wave,color_curve=None,band1=None,band2=None,ref_color=False):
    band=None
    for b in np.unique(source.lc['band']):
        temp=sncosmo.get_bandpass(b)
        if temp.minwave()<=min(wave) and temp.maxwave()>=max(wave):
            band=sncosmo.get_bandpass(b)
    if band is None:
        raise RuntimeError("Hmm, your data do not contain the band you want to fit.")
    finalPhase=source._phase
    if color_curve is not None and not ref_color:
        zp1=source.lc['zp'][source.lc['band']==band1][0]
        zp2=source.lc['zp'][source.lc['band']==band2][0]
        flux1=sncosmo.Model(source._ts_sources[band1]).bandflux(band1,source._phase,zp1,
                                                                source.lc['zpsys'][source.lc['band']==band1][0])

        temp_flux=flux1/10**(-.4*(color_curve(source._phase)-(zp1-zp2)))

    else:
        temp_flux=np.ones(len(finalPhase))

    
    zpnorm = 10.**(0.4 * source.lc['zp'][source.lc['band']==band.name][0])
    

    wave, dwave = integration_grid(band.minwave(), band.maxwave(),
                                   MODEL_BANDFLUX_SPACING)

    ms = sncosmo.get_magsystem(source.lc['zpsys'][0])
    zpnorm = zpnorm / ms.zpbandflux(band)
    flux=temp_flux*HC_ERG_AA/(dwave*np.sum(wave*band(wave))*zpnorm)
    finalWave=np.arange(wave[0]-dwave*10,wave[-1]+dwave*10,dwave)
    finalFlux=np.zeros((len(finalPhase),len(finalWave)))

    for i in range(len(finalPhase)):
        #if finalPhase[i]>= np.min(timeArr) and finalPhase[i] <= np.max(timeArr):
        for j in range(len(finalWave)):
            if finalWave[j]>=wave[0] and finalWave[j]<=wave[-1]:
                finalFlux[i][j]=flux[i]


    #offset=np.zeros(len(finalPhase))

    out=sncosmo.TimeSeriesSource(np.array(finalPhase),np.array(finalWave),np.array(finalFlux),zero_before=False)
    return (out)






class PierelSource(sncosmo.Source):
    _param_names=['amplitude','k','sigma','s']
    param_names_latex=['A','K','\sigma','Shift']

    def __init__(self,data,name='PierelSource', version=None,tstep=1):
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
        #self._phase=np.arange(0,np.max(data['time'])-np.min(data['time']),1)
        #self._phase=np.arange(-(np.max(data['time'])-np.min(data['time'])),np.max(data['time'])-np.min(data['time']),tstep)
        self._phase=np.arange(-800,800,1)
        self._parameters=np.array([1.,1.,1.,0.])
        self._tstep=tstep
        self._ts_sources={b:_param_to_source(self,self._phase,sncosmo.get_bandpass(b).wave) for b in np.unique(self.lc['band'])}


    def _param_flux(self,phase):
        temp=exponnorm.pdf(phase,self._parameters[1],scale=self._parameters[2])
        if np.max(temp)==0:
            return(np.zeros(len(phase)))
        pierelFlux=self._parameters[0]*exponnorm.pdf(phase,self._parameters[1],loc=-phase[temp==np.max(temp)],scale=self._parameters[2])/np.max(temp)+self._parameters[3]

        #print(phase[0],self._parameters[4])
        if np.inf in pierelFlux or np.any(np.isnan(pierelFlux)):
            #print(self._parameters)
            return(np.zeros(len(phase)))
        return(pierelFlux)

    def _flux(self,phase,wave):

        band=[x for x in np.unique(self.lc['band']) if sncosmo.get_bandpass(x).wave[0]<=wave[0] and sncosmo.get_bandpass(x).wave[-1]>=wave[-1]][0]

        src=self._ts_sources[band]

        return(src._flux(phase,wave)*(self._param_flux(phase)[:,None]))

class BazinSource(sncosmo.Source):

    _param_names=['amplitude','B','fall','rise']
    param_names_latex=['A','B','t_{fall}','t_{rise}']

    def __init__(self,data,name='BazinSource', version=None,tstep=1,colorCurve=None):
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
        #self._phase=np.arange(-(np.max(data['time'])-np.min(data['time'])),np.max(data['time'])-np.min(data['time']),tstep)
        self._phase=np.arange(-300,300,1)
        self._parameters=np.array([1.,0.,30.,15.])
        self._tstep=tstep
        if colorCurve is not None:
            color=[x for x in colorCurve.colnames if x!='time'][0]
            curve=interp1d(colorCurve['time'],colorCurve[color],fill_value=0.,bounds_error=False)
            self._ts_sources=dict([])
            i=0
            for b in [color[0:color.find('-')],color[color.find('-')+1:]]:
                self._ts_sources[b]=_param_to_source(self,sncosmo.get_bandpass(b).wave,curve,
                                                     color[0:color.find('-')],color[color.find('-')+1:],ref_color=i==0)
                i+=1
        else:
            self._ts_sources={b:_param_to_source(self,sncosmo.get_bandpass(b).wave) for b in np.unique(self.lc['band'])}


    def _constantBazin(self,length,B):
        return(np.ones(length)*B)

    def _param_flux(self,phase):
        if self._parameters[0]==0:
            return(self._constantBazin(len(phase),self._parameters[1]))
        elif self._parameters[2]<=self._parameters[3]:
            return(np.zeros(len(phase)))
        else:
            temp=(np.exp(-phase/self._parameters[2])/(1+np.exp(-phase/self._parameters[3])))
        if np.max(temp)==0:
            return(np.zeros(len(phase)))

        bazinFlux=self._parameters[0]*(np.exp(-(phase+np.median(phase[np.where(temp==np.max(temp))[0]]))/self._parameters[2])/\
            (1+np.exp(-(phase+np.median(phase[np.where(temp==np.max(temp))[0]]))/self._parameters[3])))/np.max(temp) + self._parameters[1]
        


        if np.inf in bazinFlux or np.any(np.isnan(bazinFlux)):
            return(np.zeros(len(phase)))
        return(bazinFlux)



    def _flux(self,phase,wave):

        band=[x for x in np.unique(self.lc['band']) if sncosmo.get_bandpass(x).wave[0]<=wave[0] and sncosmo.get_bandpass(x).wave[-1]>=wave[-1]][0]

        src=self._ts_sources[band]

        return(src._flux(phase,wave)*(self._param_flux(phase)[:,None]))



class NewlingSource(sncosmo.Source):
    _param_names=['A','psi','sigma','k','phi']
    param_names_latex=['A','\psi','\sigma','k','phi']

    def __init__(self,data,name='NewlingSource', version=None,tstep=1,flip=False):
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

    def __init__(self,data,name='KarpenkaSource', version=None,tstep=1):
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