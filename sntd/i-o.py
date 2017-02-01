#!/Users/jpierel/anaconda3/envs/astro2/bin python2
from pycs.gen.lc import lightcurve
from astropy.table import Table,vstack
import numpy as np
from .util import anyOpen
import warnings

__props__={
    ('mjd', 'mjdobs', 'jd', 'time', 'date', 'mjd_obs','mhjd'):'time',
    ('filter', 'band', 'flt', 'bandpass'): 'band',
    ('flux', 'f'):'flux',
    ('flux_error', 'fluxerr', 'fluxerror', 'fe', 'flux_err'):'fluxerr',
    ('zero_point','zp', 'zpt', 'zeropoint'):'zp',
    ('zpsys', 'magsys', 'zpmagsys'):'zpsys',
    ('mag','magnitude'):'mag',
    ('magerr','magerror','magnitudeerror','magnitudeerr'):'magerr',
    ('ra','Right Ascension'):'ra',
    ('dec','declination'):'dec',
    ('magcol','magnitudeCol','magColNumber'):'magcol',
    ('errcol','magerrcol'):'errcol',
    ('fluxcol','fluxColNum'):'fluxcol',
    ('fluxerrcol','fluxerrorcol'):'fluxerrcol',
    ('startline','firstline'):'startline',
    ('jdcol','jdcolnum'):'jdcol'
}

class curve(lightcurve):
    """
    A class, inheriting from pycs lightcurve superclass, that now also has an
    astropy.table.Table version of the data file for SNCosmo commands and flux/fluxerr
    arrays.
    """
    def __init__(self,table=None):
        """
        Constructor for curve class, which inherits from the PyCS lightcurve
        superclass. Following their format, we initialize a small lightcurve with
        only 5 points, etc. See pycs.gen.lc module for more info on the inherited properties.

        @TODO: implement more general read and right functions

        :param table: astropy.table.Table data structure for use with SNCosmo
        """

        self.table=table
        """@type: ~astropy.table.Table
        @ivar: A table containing the data, used for SNCosmo functions, etc.
        """
        #We just start with default test values
        self.fluxes=np.array([0.0,0.25,0.5,0.75,1.0])
        """@type: 1D float array
        @ivar: Fluxes, preferably connected to magnitudes that exist in the superclass
        """
        self.fluxerrs=np.array([0.1,0.1,0.15,0.05,0.2])
        """@type: 1D float array
        @ivar: Errors of fluxes, as floats
        """

    def calcmagshiftfluxes(self,inverse=False):
        """
        Returns an array of fluxes that you have to add to the curve.fluxes if you want to take into account for the magshift.
        Be careful with the order of things ! Always take into account this magshift *before* the magshift / microlensing.

        :param inverse: Returns the fluxes to subtract from magshiftedfluxes to get back the original fluxes.
        :return: A 1-D float array of length len(self), containing the shifts in flux from self.magshift
        """
        if self.magshift != 0.0:
            # We need to include a check :
            if inverse:
                shifts = -((10**(self.magshift*np.ones(len(self))/2.5)-1)*_mag_to_flux(self.mags))
            else:
                shifts = ((10 ** (-self.magshift*np.ones(len(self)) / 2.5) - 1) * _mag_to_flux(self.mags))
            if np.all(np.isnan(shifts) == False) == False:  # If there is a nan in this...
                print("Ouch, negative flux !")
                return np.zeros(len(self))
            else:
                return shifts
        else:
            return np.zeros(len(self))

    def getfluxes(self,noml=False):
        """
        A getter method returning fluxes, but taking magshift, fluxshift, and microlensing into consideration ("copy" of
        lightcurve.getmags()

        @TODO: More advanced mag_to_flux integration or add ml.calcmlfluxes

        :param noml: to ignore microlensing effects while calculating the flux with fluxshift/magshift
        :return: 1-D float array of length len(self), a copy of shifted fluxes (doesn't effect self.fluxes)
        """
        if self.magshift != 0.0:
            if (self.ml != None) and (noml == False):
                return self.fluxes + self.calcmagshiftfluxes() + self.fluxshift + _mag_to_flux(self.ml.calcmlmags(self))
            else:
                return self.fluxes + self.calcmagshiftfluxes() + self.fluxshift
        else:
            if (self.ml != None) and (noml == False):
                return self.fluxes + self.fluxshift + _mag_to_flux(self.ml.calcmlmags(self))
            else:
                return self.fluxes + self.fluxshift

    def getfluxerrs(self):
        """
        A getter method returning fluxerrs

        @TODO: Figure out if magshift/ml effect fluxerrs

        :return: A copy of self.fluxerrs
        """
        return self.fluxerrs.copy()

    def gettable(self):
        """
        @TODO add info here
        :return:
        """
        temp=self.table.copy()
        temp[_get_default_prop_name('jds')]=self.getjds()
        temp[_get_default_prop_name('mags')] = self.getmags()
        temp[_get_default_prop_name('magerrs')] = self.getmagerrs()
        temp[_get_default_prop_name('flux')] = self.getfluxes()
        temp[_get_default_prop_name('fluxerr')] = self.getfluxerrs()
        return temp

    def applyml(self):
        """
        Simply extends the superclass version of applyml to add the ml effect to the fluxes as well as the mags
        :return: None, but changes self.mags,self.fluxes, and removes current microlensing
        """
        if self.ml == None:
            raise(RuntimeError, "Hey, there is no ml associated to this lightcurve !")

        if self.fluxshift != 0.0:
            raise(RuntimeError, "Apply the fluxshift before applying the ML !")
        self.fluxes+=_mag_to_flux(self.ml.calcmlmags(self))
        super(curve,self).applyml()

    def cutmask(self):
        """
        Simply extends the superclass version of cutmask to include flux
        :return: None, but changes self.jds, self.mags, self.magerrs, self.labels, self.properties, self.mask, and
        removes microlensing if it exists
        """
        self.fluxes=self.fluxes[self.mask]
        self.fluxerrs=self.fluxerrs[self.mask]
        super(curve,self).cutmask()

    def validate(self, verbose=False):
        """
        Simply extends the superclass version of validate to include flux
        :param verbose: Default=False -> won't print "Validation done!"
        :return: None
        """
        ndates=len(self)
        if len(self.fluxes) != ndates or len(self.fluxerrs) != ndates:
            raise(RuntimeError, "The length of your flux/fluxerr arrays are inconsistent with rest of curve!")
        super(curve,self).validate(verbose)

    def sort(self):
        """
        Simply extends the superclass version of sort to include flux
        :return: None, but changes self.jds, self.mags, self.magerrs, self.labels, self.properties, self.mask, and
        removes microlensing if it exists
        """
        sortedindices=np.argsort(self.jds)
        self.fluxes=self.fluxes[sortedindices]
        self.fluxerrs=self.fluxerrs[sortedindices]
        super(curve,self).sort()

    def montecarlofluxes(self,f=1.0,seed=None):
        self.commentlist.append("Monte Carlo on fluxes !")  # to avoid confusions.
        rs = np.random.RandomState(seed)  # we create a random state object, to control the seed.
        self.fluxes += rs.standard_normal(self.fluxes.shape) * f * self.fluxerrs  # and here we do the actual bootstrapping !

    def merge(self,otherlc):
        """
        This is just an extension of the existing merge function in the lightcurve
        class, which just merges the new features (i.e. Table and flux data) as well
        :param otherlc: The other curve object you're merging
        :return: A single merged curve object
        """
        concfluxes=np.concatenate([self.getfluxes(),otherlc.getfluxes()])
        concfluxerrs=np.concatenate([self.getfluxerrs(),otherlc.getfluxerrs()])
        concTable=vstack(self.gettable(),otherlc.gettable())
        self.fluxes=concfluxes
        self.fluxerrs=concfluxerrs
        self.table=concTable
        super(curve,self).merge(otherlc)


def _mag_to_flux(mags):
    return 10.0**(mags/-2.5)

def _get_default_prop_name(prop):
    try:
        temp=__props__[[item for item in __props__.keys() if prop in item or prop.lower() in [x.lower() for x in item] or prop.upper() in [y.upper() for y in item]][0]]
    except:
        return prop
    return temp
temp=curve()
temp.table=Table()
print(temp.table)
