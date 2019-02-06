import os,sys,math,subprocess,sncosmo,abc
from textwrap import dedent

import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table
from scipy.interpolate import splrep,splev
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import WMAP9 as cosmo
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.mlab as mlab
from scipy.interpolate import interp1d,interp2d
from sncosmo.models import _ModelBase
import extinction

from .util import __dir__,__current_dir__
from .mldata import MicrolensingData

__all__=['_mlProp','_mlFlux','realizeMicro','microcaustic_field_to_curve',
         'AchromaticMicrolensing','AchromaticSplineMicrolensing','ChromaticSplineMicrolensing','PeakAchromaticMicrolensing',
         '_CCM89Dust','_OD94Dust','_F99Dust']
#def identifyML(lc):
def realizeMicro(arand=.25,debug=0,kappas=.75,kappac=.15,gamma=.76,eps=.6,nray=300,minmass=10,maxmass=10,power=-2.35,pixmax=5,pixminx=0,pixminy=0,pixdif=10,fracpixd=.3,iwrite=0,verbose=False):
    types=['%.3f','%i','%.2f','%.2f','%.2f','%.3f','%i','%.6f','%.6f','%.3f','%.3f','%.3f','%.3f','%.3f','%.3f','%i']
    inData=[arand,debug,kappas,kappac,gamma,eps,nray,minmass,maxmass,power,pixmax,pixminx,pixminy,pixdif,fracpixd,iwrite]
    inputFile=np.loadtxt(os.path.join(__dir__,'microlens','default_input'),dtype='str',delimiter='tab')
    outFile=[]
    #dt=np.dtype([('a',np.float64),('b',np.unicode_),('c',np.unicode_)])
    for i in range(len(inputFile)-1):
        #dat=inputFile[i].split()

        #if len(dat)<3:
        #    outFile.append(dat[0])
        #    break

        dat=str(inData[i])

        #print(np.array([dat[0],dat[1],' '.join(dat[2:])]))


        #outFile.append([dat[0],dat[1],' '.join(dat[2:])])
        outFile.append(dat)

    #outFile=np.array(outFile)

    #print(outFile)
    #np.savetxt(os.path.join(__dir__,'microlens','input'),outFile,fmt=['%3.3f','%s','%s'],delimiter='tab')
    thefile=open(os.path.join(__dir__,'microlens','input'),'w')

    for i in range(len(outFile)-1):
        #thefile.write((types[i]+'\t\t%s\t\t%s\n')%(float(outFile[i][0]),outFile[i][1],outFile[i][2]))
        thefile.write((types[i]+'\n')%(float(outFile[i])))
    thefile.write(outFile[-1])
    thefile.close()

    num=np.loadtxt(os.path.join(__dir__,'microlens','jobnum'),dtype='str')
    try:
        os.remove(os.path.join(__dir__,'microlens','IRIS'+str(num)))
    except:
        pass
    try:
        os.remove(os.path.join(__dir__,'microlens','IRIS'+str(num)+'.fits'))
    except:
        pass
    os.chdir(os.path.join(__dir__,'microlens'))
    if verbose:
        subprocess.call(r'./microlens')
    else:
        with open(os.devnull,'w') as f:
            subprocess.call(r'./microlens',stdout=f)


    num=np.loadtxt(os.path.join(__dir__,'microlens','jobnum'),dtype='str')
    #lensPlane=np.array(fits.open(os.path.join(__dir__,'microlens','IRIS'+str(num)+'.fits'))[0].data,dtype=np.float64)
    try:
        lensPlane=fits.open(os.path.join(__dir__,'microlens','IRIS'+str(num)+'.fits'))[0].data
    except:
        print('There was an error with the inputs of your microcaustic.')
        sys.exit()
    #curve=ascii.read(os.path.join(__dir__,'microlens','out_line'),names=('t','xval','yval','pixvalue','maglin','xpix','ypix'))
    os.chdir(__current_dir__)
    return(lensPlane)



def microcaustic_field_to_curve(field,time,zl,zs,velocity=(10**4)*(u.kilometer/u.s),M=(1*u.solMass).to(u.kg),loc='Random',plot=False):

    D=cosmo.angular_diameter_distance_z1z2(zl,zs)*cosmo.angular_diameter_distance(zs)/cosmo.angular_diameter_distance(zl)
    D=D.to(u.m)
    einsteinRadius=np.sqrt(4*const.G*M*D/const.c**2)
    einsteinRadius=einsteinRadius.to(u.kilometer)
    try:
        velocity.to(u.kilometer/u.s)
    except:
        print('Assuming velocity is in km/s.')
        velocity*=(u.kilometer/u.s)
    try:
        M.to(u.kg)
    except:
        print('Assuming mass is in kg.')
    #mlimage=fits.getdata(field)
    h,w=field.shape

    height=10*einsteinRadius.value
    width=10*einsteinRadius.value
    #print(10*einsteinRadius)
    #center=(width/2,height/2)
    pixwidth=width/w
    pixheight=height/h
    if pixwidth!=pixheight:
        print('Hmm, you are not using squares...')
        sys.exit()
    maxRadius=((np.max(time)*u.d).to(u.s))*velocity
    maxRadius=maxRadius.value
    maxx=int(math.floor(maxRadius/pixwidth))
    maxy=int(math.floor(maxRadius/pixheight))
    mlimage=field[maxx:-maxx][maxy:-maxy]


    if loc=='Random' or not isinstance(loc,(list,tuple)):
        loc=(int(np.random.uniform(maxx,w-maxx)),int(np.random.uniform(maxy,h-maxy)))


    tempTime=np.array([((x*u.d).to(u.s)).value for x in time])
    snSize=velocity.value*tempTime/pixwidth



    dmag=mu_from_image(field,loc,snSize,'disk',plot,time)

    return(time,dmag)




def createCircularMask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

def createGaussMask(h,w,center=None,radius=None):
    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])
        #Set up the 2D Gaussian:
    delta = 0.025
    x = np.arange(-3.0, 3.0, delta)
    y = np.arange(-3.0, 3.0, delta)
    X, Y = np.meshgrid(x, y)
    sigma = 1.0
    Z = mlab.bivariate_normal(X, Y, sigma, sigma, 0.0, 0.0)
    #Get Z values for contours 1, 2, and 3 sigma away from peak:
    z1 = mlab.bivariate_normal(0, 1 * sigma, sigma, sigma, 0.0, 0.0)
    z2 = mlab.bivariate_normal(0, 2 * sigma, sigma, sigma, 0.0, 0.0)
    z3 = mlab.bivariate_normal(0, 3 * sigma, sigma, sigma, 0.0, 0.0)
    #plt.figure()
    #plot Gaussian:
    #im = plt.imshow(Z, interpolation='bilinear', origin='lower',
                    #extent=(-50,50,-50,50),cmap=cm.gray)
    #Plot contours at whatever z values we want:
    #CS = plt.contour(Z, [z1, z2, z3], origin='lower', extent=(-50,50,-50,50),colors='red')
    #plt.show()

class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def mu_from_image(image, center,sizes,brightness,plot,time):
    h, w = image.shape
    mu = []
    if plot:
        fig=plt.figure(figsize=(10,10))

        ax=fig.gca()
        plt.imshow(-(image-1024)/256., aspect='equal', interpolation='nearest', cmap=cm.bwr,norm=MidpointNormalize(vmin=-2,vmax=2,midpoint=0),
                  vmin=-2, vmax=2, origin='lower')

        ax.set_xticklabels([0,0,2,4,6,8,10],fontsize=14)
        ax.set_yticklabels([0,0,2,4,6,8,10],fontsize=14)
        ax.set_xlabel('$R_E$',fontsize=18,labelpad=0)
        ax.set_ylabel('$R_E$',fontsize=18)
    #for r,a in zip([snSize[l),snSize[150],snSize[-1]],[.4,.5,.7]):
    #
    #print(np.mean(image),np.std(image))
    #print(np.mean((image-1024)/256.))
    image=10**(.4*(image-1024)/256.)
    i=0
    alphas=[1,.5,.7]
    for r in sizes:
        if r in [sizes[int(len(sizes)/5)],sizes[int(len(sizes)/2)],sizes[int(len(sizes)-1)]]:

            circle = Circle(center, r, color='#004949', alpha=alphas[i])
            i+=1
            if plot:
                ax.add_patch(circle)
        if brightness=='disk':
            mask = createCircularMask(h,w,center=center,radius=r)
            try:
                totalMag=float((image[mask]).sum())/float(mask.sum())
            except:
                totalMag=0
            if totalMag==0:
                mu.append(1024)
            else:
                mu.append(totalMag)
        else:
            mask1,mask2,mask3=createGaussMask(h,w,center=center,radius=r/3)
            scale=np.array([.68,.27,.05])
            totalMags=[]
            for mask in [mask1,mask2,mask3]:
                try:
                    if mask.sum()==0:
                        totalMags.append(0)
                        continue
                    tempMag=float(image[mask].sum())/float(mask.sum())
                except RuntimeError:
                    tempMag=0
                totalMags.append(tempMag)
            if np.max(totalMags==0):
                mu.append(1024)
            else:
                mu.append(np.dot(np.array(totalMags),scale))



    mu = np.array(mu)
    mu/=np.mean(mu)
    dmag=-2.5*np.log10(mu)
    if plot:
        cbaxes = fig.add_axes([.82, 0.33, 0.04, 0.55])

        cb = plt.colorbar(cax = cbaxes)
        cb.ax.set_ylabel('Magnification (Magnitudes)',fontsize=18,rotation=270,labelpad=25)
        cb.ax.invert_yaxis()
        cb.ax.tick_params(labelsize=14)


        ax_divider = make_axes_locatable(ax)
        ax_ml = ax_divider.append_axes("bottom", size="25%", pad=.7)
        for tick in ax_ml.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in ax_ml.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        ax_ml.plot(time,dmag,ls='-',marker=' ', color='#004949')
        ax_ml.set_ylabel(r'$\Delta m$ (mag)',fontsize=18)
        ax_ml.set_xlabel('Time from Explosion (days)',fontsize=18)
        ax_ml.invert_yaxis()
        #ax.plot(sizes[10:-10],dmag[10:-10])
        #plt.savefig('sntd_microlensing.pdf',format='pdf',overwrite=True)
        #plt.show()
        #plt.clf()
        #plt.close()

    return(mu)

class AchromaticMicrolensing(sncosmo.PropagationEffect):
    """ Simulated microlensing magnification, read in from an external
    data file.  The input data file must provide a column for SN phase
    and magnification (no wavelength dependence).
    """
    _param_names = []
    param_names_latex = []
    _minwave = 0.
    _maxwave = 10.**6

    #def __init__(self, mlfilename, magformat='multiply', **kwargs):
    def __init__(self, time,dmag, sigma=None,magformat='multiply', **kwargs):
        """Read in the achromatic microlensing data file.

        magformat : str
        Format of the magnification column.  May be ``multiply`` or ``add,``
        where ``multiply`` means the magnification column provides a
        multiplicative magnification factor, mu, so the effect is applied to
        the source as flux * mu, and ``add`` means the magnification column
        provides an additive magnitude, DeltaM=-2.5*log10(mu).

        Keyword arguments are passed on to astropy.table.Table.read().
        """
        self._magformat=magformat
        self._parameters = np.array([])
        #mldata = read_mldatafile(mlfilename, magformat=magformat, **kwargs)
        mldata=MicrolensingData(data={'phase':time,'magnification':dmag},magformat=magformat)
        self.mu = mldata.magnification_interpolator() #Now always a multiplicative mu
        if sigma is not None:
            self.sigma=interp1d(time,sigma,fill_value=0.,kind='cubic',bounds_error=False)
        else:
            self.sigma=interp1d(time,np.zeros(len(time)),fill_value=0.,kind='cubic',bounds_error=False)

    def propagate(self,phase, wave, flux):
        """Propagate the magnification onto the model's flux output."""
        mu = np.expand_dims(self.mu(phase), 1)
        return flux * mu



class PeakAchromaticMicrolensing(sncosmo.PropagationEffect):
    """ Simulated microlensing magnification, read in from an external
    data file.  The input data file must provide a column for SN phase
    and magnification (no wavelength dependence).
    """
    _param_names = ['A','D']
    param_names_latex = ['A','$\Delta$']
    _minwave = 0.
    _maxwave = 10.**6

    #def __init__(self, mlfilename, magformat='multiply', **kwargs):
    def __init__(self, time, magformat='multiply', **kwargs):
        """Read in the achromatic microlensing data file.

        magformat : str
        Format of the magnification column.  May be ``multiply`` or ``add,``
        where ``multiply`` means the magnification column provides a
        multiplicative magnification factor, mu, so the effect is applied to
        the source as flux * mu, and ``add`` means the magnification column
        provides an additive magnitude, DeltaM=-2.5*log10(mu).

        Keyword arguments are passed on to astropy.table.Table.read().
        """
        self._magformat=magformat
        self._parameters = np.array([0.,1.])
        self._time=time
        #mldata = read_mldatafile(mlfilename, magformat=magformat, **kwargs)
        mldata=MicrolensingData(data={'phase':time,'magnification':self._parameters[0]*time+self._parameters[1]},magformat=magformat)
        self.mu = mldata.magnification_interpolator() #Now always a multiplicative mu


    def propagate(self,phase, wave, flux):
        """Propagate the magnification onto the model's flux output."""
        self.update_mu()
        mu = np.expand_dims(self.mu(phase), 1)

        return flux * mu

    def update_mu(self):
        m,delta=self._parameters
        mldata=MicrolensingData(data={'phase':self._time,'magnification':m*self._time+delta},magformat='multiply')
        self.mu = mldata.magnification_interpolator() #Now always a multiplicative mu




class AchromaticSplineMicrolensing(sncosmo.PropagationEffect):
    """Average of randomly anchored splines, to mimic microlensing.
    We create a mock microlensing difference curve, giving the change in
    magnitude as a function of time (no variation with wavelength).
    A set of `nspl` splines are generated, each passing through `nanchor`
    anchor points, evenly spaced in time, and with y values (representing
    Delta magnitude) randomly drawn from a normal distribution with varance
    equal to `sigmadm` squared.  The final delta magnitude curve is the mean
    of the set of these `nspl` random spline curves.
    Caveat emptor: this is just a crude hack. It looks like a reasonable
    approximation for achromatic SN microlensing, but it is not actually
    derived from a real lensing simulation.
    """
    _param_names = []
    param_names_latex = []
    _minwave = 0.
    _maxwave = 10.**6

    def __init__(self, nanchor=10, sigmadm=2.0, nspl=10):
        # self._parameters = np.array([nanchor, sigmadm, nspl])
        self._parameters = np.array([])
        self._nanchor = nanchor
        self._nspl = nspl
        self._sigmadm = sigmadm

        # Define a delta mag curve as an average of random splines
        # The time dimension spans from 0 to 1, but will be rescaled
        # when propagated onto a model, so that it stretches from the model
        # minphase to maxphase.
        splineset = []
        tarray = np.linspace(0.0, 1.0, 100)
        for i in range(nspl):
            time_anchors = np.linspace(0.0, 1.0, nanchor)
            deltam_anchors = np.random.normal(0, sigmadm, len(time_anchors))
            spl1d = Spline1d(time_anchors, deltam_anchors)
            splineset.append(spl1d(tarray))
        splmean = np.mean(np.array(splineset), 0)
        self._deltamag = interp1d(tarray, splmean)


    def propagate(self, phase,wave, flux):
        """Propagate the magnification onto the model's flux output."""
        # magnify the flux
        deltamag = np.expand_dims(self._deltamag(phase), 1)
        return flux * 10**(-0.4 * deltamag)


class ChromaticSplineMicrolensing(sncosmo.PropagationEffect):
    """Average of randomly anchored splines, to mimic microlensing.
    We create a mock microlensing difference curve as is done for
    AchromaticSplineMicrolensing, giving the change in
    magnitude as a function of time. Then we add variation with wavelength
    by adding a two-dimensional polynomial to the delta mag surface.

    Caveat emptor: this is just a crude hack. It looks like a reasonable
    approximation for SN microlensing including wavelength variation, but it
    is not actually derived from a real lensing simulation.

    Note that the "microlensing" this produces does not have any range of
    SN phase in which the microlensing is achromatic.
    """
    _param_names = []
    param_names_latex = []
    _minwave = 200.   # Angstroms
    _maxwave = 25000. # Angstroms

    def __init__(self, nanchor=100, sigmadm=2.0, nspl=100):
        # self._parameters = np.array([nanchor, sigmadm, nspl])
        self._parameters = np.array([])
        self._nanchor = nanchor
        self._nspl = nspl
        self._sigmadm = sigmadm

        nsteps = 100
        tarray = np.linspace(0., 1., nsteps)
        wavearray = np.linspace(0., 1., nsteps)
        time_anchors = np.linspace(0., 1., nanchor)
        wave_anchors = np.linspace(0., 1., nanchor)

        # First surface: make a 1d delta-mag curve that is the mean of a
        # set of random splines in the time dimension, then extend it
        # without variation into the wavelength dimension.
        splfitarray = []
        for i in range(nspl):
            deltam_anchors = np.random.normal(
                0, sigmadm, len(time_anchors))
            spl1d = Spline1d(time_anchors, deltam_anchors)
            splfitarray.append(spl1d(tarray))
        splmean = np.mean(np.array(splfitarray), 0)
        splmean_surface = np.tile(splmean, nsteps).reshape((nsteps, nsteps))

        # Second surface: a 2-D polynomial grid across time and wavelength,
        # defined with a random covariance matrix:
        # each component in the cov matrix is drawn from a normal dist. with
        # sigma = 1/4th of sigmadeltam.
        # WARNING: right now this setup is basically totally unsupported
        # by actual microlensing simulations.
        cov = np.random.normal(0, sigmadm / 4., 4).reshape(2, 2)
        polygrid_surface = np.polynomial.polynomial.polygrid2d(
            tarray, wavearray, cov)

        deltam_surface = splmean_surface + polygrid_surface
        self._deltamag = interp2d(wavearray, tarray, deltam_surface)


    def propagate(self, phase,wave, flux):
        """Propagate the magnification onto the model's flux output."""
        # magnify the flux
        wavefraction = (wave-self._minwave)/(self._maxwave-self._minwave)
        deltamag = self._deltamag(wavefraction, phase)
        return flux * 10**(-0.4 * deltamag)


def _mlFlux(self,time, wave):
    """Replacement for sncosmo Array flux function."""
    a = 1. / (1. + self._parameters[0])
    phase = (time - self._parameters[1]) * a
    minphase = (self.mintime() - self._parameters[1]) * a
    maxphase = (self.maxtime() - self._parameters[1]) * a
    restwave = wave * a

    # Note that below we multiply by the scale factor to conserve
    # bolometric luminosity.
    f = a * self._source._flux(phase, restwave)

    # Pass the flux through the PropagationEffects.
    for effect, frame, zindex in zip(self._effects, self._effect_frames,
                                     self._effect_zindicies):
        if frame == 'obs':
            effect_wave = wave
            effect_phase=phase*(1./a)
        elif frame == 'rest':
            effect_wave = restwave
            effect_phase=phase
        else:  # frame == 'free'
            effect_a = 1. / (1. + self._parameters[zindex])
            effect_wave = wave * effect_a
            effect_phase=phase/a*(1.+self._parameters[zindex])


        f = effect.propagate(effect_phase,effect_wave, f)

    return f

def _mlProp(_ModelBase):
    """Abstract base class for propagation effects.

    Derived classes must define _minwave (float), _maxwave (float).
    """

    __metaclass__ = abc.ABCMeta

    def minwave(self):
        return self._minwave

    def maxwave(self):
        return self._maxwave

    @abc.abstractmethod
    def propagate(self, wave, flux, phase=None):
        pass

    def _headsummary(self):
        summary = """\
        class           : {0}
        wavelength range: [{1:.6g}, {2:.6g}] Angstroms""" \
            .format(self.__class__.__name__, self._minwave, self._maxwave)
        return dedent(summary)

class _CCM89Dust(sncosmo.PropagationEffect):
    """Cardelli, Clayton, Mathis (1989) extinction model dust."""
    _param_names = ['ebv', 'r_v']
    param_names_latex = ['E(B-V)', 'R_V']
    _minwave = 1000.
    _maxwave = 33333.33

    def __init__(self):
        self._parameters = np.array([0., 3.1])

    def propagate(self, phase,wave, flux):
        """Propagate the flux."""
        ebv, r_v = self._parameters
        return extinction.apply(extinction.ccm89(wave, ebv * r_v, r_v), flux)


class _OD94Dust(sncosmo.PropagationEffect):
    """O'Donnell (1994) extinction model dust."""
    _param_names = ['ebv', 'r_v']
    param_names_latex = ['E(B-V)', 'R_V']
    _minwave = 909.09
    _maxwave = 33333.33

    def __init__(self):
        self._parameters = np.array([0., 3.1])

    def propagate(self, phase,wave, flux):
        """Propagate the flux."""
        ebv, r_v = self._parameters
        return extinction.apply(extinction.odonnell94(wave, ebv * r_v, r_v),
                                flux)


class _F99Dust(sncosmo.PropagationEffect):
    """Fitzpatrick (1999) extinction model dust with fixed R_V."""
    _minwave = 909.09
    _maxwave = 60000.

    def __init__(self, r_v=3.1):
        self._param_names = ['ebv']
        self.param_names_latex = ['E(B-V)']
        self._parameters = np.array([0.])
        self._r_v = r_v
        self._f = extinction.Fitzpatrick99(r_v=r_v)

    def propagate(self, phase,wave, flux):
        """Propagate the flux."""
        ebv = self._parameters[0]
        return extinction.apply(self._f(wave, ebv * self._r_v), flux)