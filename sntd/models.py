import sncosmo,sys,math
import numpy as np
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d,splrep,splev
from sncosmo.utils import integration_grid
from sncosmo.constants import HC_ERG_AA, MODEL_BANDFLUX_SPACING
import matplotlib.pyplot as plt
from .fitting import _findMax


from copy import copy

class FlexibleSpline(sncosmo.Source):


    _param_names = ['k']

    _SCALE_FACTOR = 1e-12
    param_names_latex=['k']
    def __init__(self, data,weights=None,name=None, version=None,tstep=1,wstep=10):
        self.name = name
        self.version = version
        self._model = {}
        self._phase=data['time']
        wave=[]
        for b in np.unique(data['band']):
            wave=np.append(wave,sncosmo.get_bandpass(b).wave)
        wave=np.sort(np.unique(wave))
        wave=np.append([.99*wave[0]],wave)
        wave=np.append(wave,[1.01*wave[-1]])



        self._wave=wave
        self.lc=data
        self._steps={'tstep':tstep,'wstep':wstep}
        #self._time = np.linspace(np.min)
        #self._wave = wave
        if not weights:
            self._weights=None
        #    self._weights=np.ones(len(data))
        self._parameters = np.array([3])#,_findMax(time,flux)[0]])

        #if self._parameters[1] not in self._time:
        #    self._t0Set()

    def _dataToSource(self,phase,waveArr):
        finalWave,finalPhase,finalFlux=None,None,[]
        bands=[b for b in np.unique(self.lc['band']) if waveArr[0] >= sncosmo.get_bandpass(b).minwave() and waveArr[-1] <= sncosmo.get_bandpass(b).maxwave()]
        #for b in np.unique(self.lc['band']):
        #    print(sncosmo.get_bandpass(b).minwave(),waveArr[0],sncosmo.get_bandpass(b).maxwave(),waveArr[-1])
        for b in bands:

            t0Est,t0Flux=_findMax(self.lc['time'][self.lc['band']==b],self.lc['flux'][self.lc['band']==b])
            timeInterval=[float(np.min(self.lc['time'][self.lc['band']==b])),float(np.max(self.lc['time'][self.lc['band']==b]))]
            phase=phase[phase>=timeInterval[0]]
            phase=phase[phase<=timeInterval[-1]]
            tempTime=np.array(self.lc['time'][self.lc['band']==b])
            tempFlux=np.array(self.lc['flux'][self.lc['band']==b])



            if t0Est is not None and t0Est not in self.lc['time'][self.lc['band']==b]:
                tempTime=np.append(tempTime,[t0Est])
                tempFlux=np.append(tempFlux,[t0Flux])
                srt=np.argsort(tempTime)
                tempTime=tempTime[srt]
                tempFlux=tempFlux[srt]
            if not self._weights:
                weights=np.ones(len(tempFlux))
            else:
                weights=self._weights
            spl=splrep(tempTime,tempFlux,k=int(self._parameters[0]),w=1/weights)
            band=sncosmo.get_bandpass(b)
            timeArr=np.arange(timeInterval[0],timeInterval[1]+self._steps['tstep'],self._steps['tstep'])
            wave, dwave = integration_grid(band.minwave(), band.maxwave(),
                                           MODEL_BANDFLUX_SPACING)
            #timeArr-=t0Est
            #print(timeArr)
            #print(phase[0],timeInterval[0],phase[-1],timeInterval[-1])
            if phase[0] <timeInterval[0] or phase[-1] > timeInterval[-1]:
                raise RuntimeError('The phase you requested is outside your data bounds.')
            zpnorm = 10.**(0.4 * self.lc['zp'][0])


            ms = sncosmo.get_magsystem(self.lc['zpsys'][0])
            zpnorm = zpnorm / ms.zpbandflux(b)

            flux=splev(timeArr,spl)*HC_ERG_AA/(dwave*np.sum(wave*band(wave))*zpnorm)

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
        print('here')
        model=self._dataToSource(phase,wave)
        return(model._flux(phase,wave))


    def _bandflux_rvar_single(self, band, phase):
        """Model relative variance for a single bandpass."""

        # Raise an exception if bandpass is out of model range.
        if (band.minwave() < self._wave[0] or band.maxwave() > self._wave[-1]):
            raise ValueError('bandpass {0!r:s} [{1:.6g}, .., {2:.6g}] '
                             'outside spectral range [{3:.6g}, .., {4:.6g}]'
                             .format(band.name, band.wave[0], band.wave[-1],
                                     self._wave[0], self._wave[-1]))

        x1 = self._parameters[1]

        # integrate m0 and m1 components
        wave, dwave = integration_grid(band.minwave(), band.maxwave(),
                                       MODEL_BANDFLUX_SPACING)
        trans = band(wave)
        m0 = self._model['M0'](phase, wave)
        m1 = self._model['M1'](phase, wave)
        tmp = trans * wave
        f0 = np.sum(m0 * tmp, axis=1) * dwave / HC_ERG_AA
        m1int = np.sum(m1 * tmp, axis=1) * dwave / HC_ERG_AA
        ftot = f0 + x1 * m1int

        # In the following, the "[:,0]" reduces from a 2-d array of shape
        # (nphase, 1) to a 1-d array.
        lcrv00 = self._model['LCRV00'](phase, band.wave_eff)[:, 0]
        lcrv11 = self._model['LCRV11'](phase, band.wave_eff)[:, 0]
        lcrv01 = self._model['LCRV01'](phase, band.wave_eff)[:, 0]
        scale = self._model['errscale'](phase, band.wave_eff)[:, 0]

        v = lcrv00 + 2.0 * x1 * lcrv01 + x1 * x1 * lcrv11

        # v is supposed to be variance but can go negative
        # due to interpolation.  Correct negative values to some small
        # number. (at present, use prescription of snfit : set
        # negatives to 0.0001)
        v[v < 0.0] = 0.0001

        result = v * (f0 / ftot)**2 * scale**2

        # treat cases where ftot is negative the same as snfit
        result[ftot <= 0.0] = 10000.

        return result

    def bandflux_rcov(self, band, phase):
        """Return the *relative* model covariance (or "model error") on
        synthetic photometry generated from the model in the given restframe
        band(s).
        This model covariance represents the scatter of real SNe about
        the model.  The covariance matrix has two components. The
        first component is diagonal (pure variance) and depends on the
        phase :math:`t` and bandpass central wavelength
        :math:`\lambda_c` of each photometry point:
        .. math::
           (F_{0, \mathrm{band}}(t) / F_{1, \mathrm{band}}(t))^2
           S(t, \lambda_c)^2
           (V_{00}(t, \lambda_c) + 2 x_1 V_{01}(t, \lambda_c) +
            x_1^2 V_{11}(t, \lambda_c))
        where the 2-d functions :math:`S`, :math:`V_{00}`, :math:`V_{01}`,
        and :math:`V_{11}` are given by the files ``errscalefile``,
        ``lcrv00file``, ``lcrv01file``, and ``lcrv11file``
        respectively and :math:`F_0` and :math:`F_1` are given by
        .. math::
           F_{0, \mathrm{band}}(t) = \int_\lambda M_0(t, \lambda)
                                     T_\mathrm{band}(\lambda)
                                     \\frac{\lambda}{hc} d\lambda
        .. math::
           F_{1, \mathrm{band}}(t) = \int_\lambda
                                     (M_0(t, \lambda) + x_1 M_1(t, \lambda))
                                     T_\mathrm{band}(\lambda)
                                     \\frac{\lambda}{hc} d\lambda
        As this first component can sometimes be negative due to
        interpolation, there is a floor applied wherein values less than zero
        are set to ``0.01**2``. This is to match the behavior of the
        original SALT2 code, snfit.
        The second component is block diagonal. It has
        constant covariance between all photometry points within a
        bandpass (regardless of phase), and no covariance between
        photometry points in different bandpasses:
        .. math::
           CD(\lambda_c)^2
        where the 1-d function :math:`CD` is given by the file ``cdfile``.
        Adding these two components gives the *relative* covariance on model
        photometry.
        Parameters
        ----------
        band : `~numpy.ndarray` of `~sncosmo.Bandpass`
            Bandpasses of observations.
        phase : `~numpy.ndarray` (float)
            Phases of observations.
        Returns
        -------
        rcov : `~numpy.ndarray`
            Model relative covariance for given bandpasses and phases.
        """

        # construct covariance array with relative variance on diagonal
        diagonal = np.zeros(phase.shape, dtype=np.float64)
        for b in set(band):
            mask = band == b
            diagonal[mask] = self._bandflux_rvar_single(b, phase[mask])
        result = np.diagflat(diagonal)

        # add kcorr errors
        for b in set(band):
            mask1d = band == b
            mask2d = mask1d * mask1d[:, None]  # mask for result array
            kcorrerr = self._colordisp(b.wave_eff)
            result[mask2d] += kcorrerr**2

        return result

