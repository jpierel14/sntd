import os
import sys
import math
import subprocess
import sncosmo
import abc
from textwrap import dedent

import numpy as np
from astropy.io import fits

from astropy import units as u
from astropy import constants as const
from astropy.cosmology import WMAP9 as cosmo
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.mlab as mlab
from scipy.interpolate import interp1d, interp2d
from sncosmo.models import _ModelBase
import extinction

from .util import _filedir_, _current_dir_
from .mldata import MicrolensingData

__all__ = ['_mlProp', '_mlFlux', 'realizeMicro', 'microcaustic_field_to_curve',
           'AchromaticMicrolensing',
           'ChromaticFilterMicrolensing',
           '_CCM89Dust', '_OD94Dust', '_F99Dust']
# def identifyML(lc):


def realizeMicro(arand=.25, debug=0, kappas=.75, kappac=.15, gamma=.76, eps=.6, nray=300, minmass=10, maxmass=10, power=-2.35,
                 pixmax=5, pixminx=0, pixminy=0, pixdif=10, fracpixd=.3, iwrite=0, verbose=False):
    """
    Creates a microcaustic realization based on Wambsganss 1990 microlens code. All parameters are optional as they
    have defaults, see Wambsganss documentation for details on parameters.
    """
    types = ['%.3f', '%i', '%.2f', '%.2f', '%.2f', '%.3f', '%i', '%.6f',
             '%.6f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%i']
    inData = [arand, debug, kappas, kappac, gamma, eps, nray, minmass,
              maxmass, power, pixmax, pixminx, pixminy, pixdif, fracpixd, iwrite]
    inputFile = np.loadtxt(os.path.join(
        _filedir_, 'microlens', 'default_input'), dtype='str', delimiter='tab')
    outFile = []
    for i in range(len(inputFile)-1):

        dat = str(inData[i])

        outFile.append(dat)

    thefile = open(os.path.join(_filedir_, 'microlens', 'input'), 'w')

    for i in range(len(outFile)-1):
        thefile.write((types[i]+'\n') % (float(outFile[i])))
    thefile.write(outFile[-1])
    thefile.close()

    num = '001'
    try:
        os.remove(os.path.join(_filedir_, 'microlens', 'IRIS'+str(num)))
    except:
        pass
    try:
        os.remove(os.path.join(_filedir_, 'microlens', 'IRIS'+str(num)+'.fits'))
    except:
        pass
    try:
        os.remove(os.path.join(_filedir_, 'microlens', 'dat.'+str(num)))
    except:
        pass
    os.chdir(os.path.join(_filedir_, 'microlens'))
    if verbose:
        subprocess.call(r'./microlens')
    else:
        with open(os.devnull, 'w') as f:
            subprocess.call(r'./microlens', stdout=f)

    np.savetxt(os.path.join(_filedir_, 'microlens', 'jobnum'), [num], fmt='%s')

    try:
        fitsFile = fits.open(os.path.join(
            _filedir_, 'microlens', 'IRIS'+str(num)+'.fits'))
        lensPlane = fitsFile[0].data
        fitsFile.close()

    except RuntimeError:
        print('There was an error with the inputs of your microcaustic.')
        sys.exit()
    os.chdir(_current_dir_)
    return(lensPlane)


def microcaustic_field_to_curve(field, time, zl, zs, velocity=(10**4)*(u.kilometer/u.s), M=(1*u.solMass).to(u.kg), width_in_einstein_radii=10,
                                loc='Random', plot=False, ax=None, showCurve=True, rescale=True):
    """
    Convolves an expanding photosphere (achromatic disc) with a microcaustic to generate a magnification curve.

    Parameters
    ----------
    field:  :class:`numpy.ndarray`
        An opened fits file of a microcaustic, can be generated by realizeMicro
    time: :class:`numpy.array`
        Time array you want for microlensing magnification curve, explosion time is 0
    zl: float
        redshift of the lens
    zs: float
        redshift of the source
    velocity: float* :class:`astropy.units.Unit`
        The average velocity of the expanding photosphere
    M: float* :class:`~astropy.units.Unit`
        The mass of the deflector
    width_in_einstein_radii: float
        The width of your map in units of Einstein radii
    loc: str or tuple
        Random is defualt for location of the supernova, or pixel (x,y) coordiante can be specified
    plot: bool
        If true, plots the expanding photosphere on the microcaustic
    ax: `~matplotlib.pyplot.axis`
        An optional axis object to plot on. If you want to show the curve, this should be a list
        like this: [main_ax,lower_ax]
    showCurve: bool
        If true, the microlensing curve is plotted below the microcaustic
    rescale: bool
        If true, assumes image needs to be rescaled: (x-1024)/256
    Returns
    -------
    time: :class:`numpy.array`
        The time array for the magnification curve
    dmag: :class:`numpy.array`
        The magnification curve.
    """

    D = cosmo.angular_diameter_distance_z1z2(
        zl, zs)*cosmo.angular_diameter_distance(zs)/cosmo.angular_diameter_distance(zl)
    D = D.to(u.m)
    try:
        M.to(u.kg)
    except:
        print('Assuming mass is in kg.')
    einsteinRadius = np.sqrt(4*const.G*M*D/const.c**2)
    einsteinRadius = einsteinRadius.to(u.kilometer)
    try:
        velocity.to(u.kilometer/u.s)
    except:
        print('Assuming velocity is in km/s.')
        velocity *= (u.kilometer/u.s)

    h, w = field.shape

    height = width_in_einstein_radii*einsteinRadius.value
    width = width_in_einstein_radii*einsteinRadius.value

    pixwidth = width/w
    pixheight = height/h
    if pixwidth != pixheight:
        print('Hmm, you are not using squares...')
        sys.exit()
    maxRadius = ((np.max(time)*u.d).to(u.s))*velocity
    maxRadius = maxRadius.value
    maxx = int(math.floor(maxRadius/pixwidth))
    maxy = int(math.floor(maxRadius/pixheight))
    mlimage = field[maxx:-maxx][maxy:-maxy]

    if loc == 'Random' or not isinstance(loc, (list, tuple)):
        loc = (int(np.random.uniform(maxx, w-maxx)),
               int(np.random.uniform(maxy, h-maxy)))

    tempTime = np.array([((x*u.d).to(u.s)).value for x in time])
    snSize = velocity.value*tempTime/pixwidth

    dmag = mu_from_image(field, loc, snSize, 'disk', plot,
                         time, ax, showCurve, rescale, width_in_einstein_radii)

    return(time, dmag)


def createCircularMask(h, w, center=None, radius=None):

    if center is None:  # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None:  # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask


def createGaussMask(h, w, center=None, radius=None):
    if center is None:  # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None:  # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])
        # Set up the 2D Gaussian:
    delta = 0.025
    x = np.arange(-3.0, 3.0, delta)
    y = np.arange(-3.0, 3.0, delta)
    X, Y = np.meshgrid(x, y)
    sigma = 1.0
    Z = mlab.bivariate_normal(X, Y, sigma, sigma, 0.0, 0.0)
    # Get Z values for contours 1, 2, and 3 sigma away from peak:
    z1 = mlab.bivariate_normal(0, 1 * sigma, sigma, sigma, 0.0, 0.0)
    z2 = mlab.bivariate_normal(0, 2 * sigma, sigma, sigma, 0.0, 0.0)
    z3 = mlab.bivariate_normal(0, 3 * sigma, sigma, sigma, 0.0, 0.0)
    # plt.figure()
    # plot Gaussian:
    # im = plt.imshow(Z, interpolation='bilinear', origin='lower',
    # extent=(-50,50,-50,50),cmap=cm.gray)
    # Plot contours at whatever z values we want:
    #CS = plt.contour(Z, [z1, z2, z3], origin='lower', extent=(-50,50,-50,50),colors='red')
    # plt.show()


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


def mu_from_image(image, center, sizes, brightness, plot, time, ax, showCurve, rescale, width_in_einstein_radii):
    h, w = image.shape
    mu = []
    if rescale:
        image = 10**(.4*(image-1024)/256.)
    if plot:
        if ax is None:
            fig = plt.figure(figsize=(10, 10))

            ax = fig.gca()

        ax.imshow(-2.5*np.log10(image), aspect='equal', interpolation='nearest', cmap=cm.bwr,
                  norm=MidpointNormalize(vmin=-2, vmax=2, midpoint=0),
                  vmin=-2, vmax=2, origin='lower')

        # ,tuple(np.linspace(0,width_in_einstein_radii,5).astype(str)))
        ax.set_xticks(tuple(np.linspace(0, image.shape[0], 5)))
        # ,tuple(np.linspace(0,width_in_einstein_radii,5)))
        ax.set_yticks((np.linspace(0, image.shape[0], 5)))
        ax.set_xticklabels(np.linspace(
            0, width_in_einstein_radii, 5), fontsize=14)
        ax.set_yticklabels(np.linspace(
            0, width_in_einstein_radii, 5), fontsize=14)
        ax.set_xlabel('$R_E$', fontsize=18, labelpad=0)
        ax.set_ylabel('$R_E$', fontsize=18)

    i = 0
    alphas = [1, .5, .7]
    for r in sizes:
        if r in [sizes[int(len(sizes)/5)], sizes[int(len(sizes)/2)], sizes[int(len(sizes)-1)]]:

            circle = Circle(center, r, color='#004949', alpha=alphas[i])
            i += 1
            if plot:
                ax.add_patch(circle)
        if brightness == 'disk':
            mask = createCircularMask(h, w, center=center, radius=r)
            try:
                totalMag = float((image[mask]).sum())/float(mask.sum())
            except:
                totalMag = 0
            if totalMag == 0:
                mu.append(1024)
            else:
                mu.append(totalMag)
        else:
            mask1, mask2, mask3 = createGaussMask(
                h, w, center=center, radius=r/3)
            scale = np.array([.68, .27, .05])
            totalMags = []
            for mask in [mask1, mask2, mask3]:
                try:
                    if mask.sum() == 0:
                        totalMags.append(0)
                        continue
                    tempMag = float(image[mask].sum())/float(mask.sum())
                except RuntimeError:
                    tempMag = 0
                totalMags.append(tempMag)
            if np.max(totalMags == 0):
                mu.append(1024)
            else:
                mu.append(np.dot(np.array(totalMags), scale))

    mu = np.array(mu)
    mu /= np.median(mu)
    dmag = -2.5*np.log10(mu)
    if plot:
        if ax is None:

            cbaxes = fig.add_axes([.82, 0.33, 0.04, 0.55])
            cb = plt.colorbar(cax=cbaxes)
            cb.ax.set_ylabel('Magnification (Magnitudes)',
                             fontsize=18, rotation=270, labelpad=25)
            cb.ax.invert_yaxis()
            cb.ax.tick_params(labelsize=14)

        if showCurve:
            ax_divider = make_axes_locatable(ax)
            ax_ml = ax_divider.append_axes("bottom", size="25%", pad=.7)
            for tick in ax_ml.xaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            for tick in ax_ml.yaxis.get_major_ticks():
                tick.label.set_fontsize(14)
            ax_ml.plot(time, dmag, ls='-', marker=' ', color='#004949')
            ax_ml.set_ylabel(r'$\Delta m$ (mag)', fontsize=18)
            ax_ml.set_xlabel('Time from Explosion (days)', fontsize=18)
            ax_ml.invert_yaxis()
        # ax.plot(sizes[10:-10],dmag[10:-10])
        # plt.savefig('sntd_microlensing.pdf',format='pdf',overwrite=True)
        # plt.show()
        # plt.clf()
        # plt.close()

    return(mu)


class AchromaticMicrolensing(sncosmo.PropagationEffect):
    """
    An achromatic microlensing object, defined filter to filter.
    """
    _param_names = []
    param_names_latex = []
    _minwave = 0.
    _maxwave = 10.**6

    # def __init__(self, mlfilename, magformat='multiply', **kwargs):
    def __init__(self, time, dmag, magformat='multiply', **kwargs):
        """
        Parameters
            ----------
            time: :class:`~list` or :class:`~numpy.array`
                A time array for your microlensing
            dmag: :class:`~list` or :class:`~numpy.array`
                microlensing magnification
        magformat : str
            Format of the magnification column.  May be ``multiply`` or ``add,``
            where ``multiply`` means the magnification column provides a
            multiplicative magnification factor, mu, so the effect is applied to
            the source as flux * mu, and ``add`` means the magnification column
            provides an additive magnitude, DeltaM=-2.5*log10(mu).

        """
        self._magformat = magformat
        self._parameters = np.array([])
        mldata = MicrolensingData(
            data={'phase': time, 'magnification': dmag}, magformat=magformat)
        self.mu = mldata.magnification_interpolator()  # Now always a multiplicative mu

    def propagate(self, phase, wave, flux):
        """
        Propagate the magnification onto the model's flux output.
        """
        mu = np.expand_dims(self.mu(phase), 1)
        return flux * mu


def getBandNorm(band):
    band = sncosmo.get_bandpass(band)
    wave, dwave = sncosmo.utils.integration_grid(band.minwave(), band.maxwave(),
                                                 sncosmo.constants.MODEL_BANDFLUX_SPACING)
    trans = band(wave)
    f = np.ones((1, len(wave)))

    return np.sum(wave * trans * f, axis=1) * dwave / sncosmo.constants.HC_ERG_AA


class ChromaticFilterMicrolensing(sncosmo.PropagationEffect):
    """
    A chromatic microlensing object, defined filter to filter.
    """
    _param_names = []
    param_names_latex = []
    _minwave = 0.
    _maxwave = 10.**6

    def __init__(self, times, dmags, bands, magformat='multiply', **kwargs):
        """
        Parameters
            ----------
            times: 1D or 2D :class:`~list` or :class:`~numpy.array`
                A list of the time arrays for your microlensing, with nrows=len(bands), ncols=len(dmags)
            dmags: 1D or 2D :class:`~list` or :class:`~numpy.array`
                same as times, but for microlensing magnification
            bands: :class:`~list` or :class:`~numpy.array`
                list of bands defining microlensing
        magformat : str
            Format of the magnification column.  May be ``multiply`` or ``add,``
            where ``multiply`` means the magnification column provides a
            multiplicative magnification factor, mu, so the effect is applied to
            the source as flux * mu, and ``add`` means the magnification column
            provides an additive magnitude, DeltaM=-2.5*log10(mu).
        """
        if not isinstance(bands, (list, tuple, np.ndarray)):
            bands = [bands]
        if len(bands) != len(times):
            times = [times]*len(bands)
        if len(bands) != len(dmags):
            dmags = [dmags]*len(bands)
        self.bandNorms = [getBandNorm(b) for b in bands]
        self._magformat = magformat
        self._parameters = np.array([])

        ml_list = [MicrolensingData(data={'phase': times[i], 'magnification':dmags[i]}, magformat=magformat).magnification_interpolator() for i in
                   range(len(bands))]

        self.bandwaves = [[sncosmo.get_bandpass(band).wave[0],
                           sncosmo.get_bandpass(band).wave[-1]] for band in bands]
        self.bandtimes = [[t[0], t[-1]] for t in times]
        self.wave = np.arange(np.max([self._minwave, np.min(self.bandwaves)]), np.min([self._maxwave, np.max(self.bandwaves)])+.01,
                              sncosmo.constants.MODEL_BANDFLUX_SPACING)
        self.phase = np.arange(np.min(times), np.max(times)+.01, .1)

        all_mu = np.ones((len(self.phase), len(self.wave)))
        for i in range(len(bands)):
            indsp = np.where(np.logical_and(
                self.phase >= self.bandtimes[i][0], self.phase <= self.bandtimes[i][1]))[0]
            indsw = np.where(np.logical_and(
                self.wave >= self.bandwaves[i][0], self.wave <= self.bandwaves[i][1]))[0]
            for ind in indsp:
                all_mu[ind, indsw] = ml_list[i](
                    self.phase[ind])*np.ones(len(indsw))  # /self.bandNorms[i]

        self.mu = interp2d(self.phase, self.wave, np.array(
            all_mu).T, fill_value=1, bounds_error=True)

    def propagate(self, phase, wave, flux):
        """
        Propagate the magnification onto the model's flux output.
        """

        return flux * self.mu(phase, wave).T


def _mlFlux(self, time, wave):
    """Replacement for sncosmo Array flux function."""
    a = 1. / (1. + self._parameters[0])
    phase = (time - self._parameters[1]) * a
    minphase = (self.mintime() - self._parameters[1]) * a
    maxphase = (self.maxtime() - self._parameters[1]) * a
    restwave = wave * a
    # Note that below we multiply by the scale factor to conserve
    # bolometric luminosity.

    f = a * self._source._flux(np.round(phase, 1).astype(np.double),
                               np.round(restwave, 1).astype(np.double))

    # Pass the flux through the PropagationEffects.
    for effect, frame, zindex in zip(self._effects, self._effect_frames,
                                     self._effect_zindicies):
        if frame == 'obs':
            effect_wave = wave
            effect_phase = phase*(1./a)
        elif frame == 'rest':
            effect_wave = restwave
            effect_phase = phase
        else:  # frame == 'free'
            effect_a = 1. / (1. + self._parameters[zindex])
            effect_wave = wave * effect_a
            effect_phase = phase/a*(1.+self._parameters[zindex])

        f = effect.propagate(effect_phase, effect_wave, f)

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

    def propagate(self, phase, wave, flux):
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

    def propagate(self, phase, wave, flux):
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

    def propagate(self, phase, wave, flux):
        """Propagate the flux."""
        ebv = self._parameters[0]
        return extinction.apply(self._f(wave, ebv * self._r_v), flux)
