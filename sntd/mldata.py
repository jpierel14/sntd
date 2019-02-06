# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Convenience functions for microlensing data."""
from __future__ import division

from collections import OrderedDict
import copy
from scipy.interpolate import interp1d, interp2d

import numpy as np
from astropy.table import Table
from os import path

from sncosmo.utils import alias_map

__all__ = []

MLDATA_ALIASES = OrderedDict([
    ('time', {'phase', 'day', 'time', 't', 'date',
               'jd', 'mjd', 'mjdobs', 'mjd_obs'}),
    ('magnification', {'magnification', 'mu', 'mag', 'deltamag',
                       'dm', 'deltam'}),
    ('wavelength', {'wavelength', 'lambda', 'wave', 'w', 'lam'}),
    ])

ACHROMATIC_MLDATA_REQUIRED_ALIASES = ('phase', 'magnification')
CHROMATIC_MLDATA_REQUIRED_ALIASES = ('phase', 'magnification', 'wavelength')

class MicrolensingData(object):
    """Internal standardized representation of microlensing magnification
    data table. Applicable both for the chromatic microlensing case, where
    magnification depends on wavelength, and achromatic microlensing, where it
    does not.

    For achromatic microlensing, this class has attributes ``phase`` and
    ``magnification``, which are both numpy arrays of the same length sorted
    by ``phase``.

    For chromatic microlensing, a ``wavelength`` attribute is added as another
    numpy array.

    Parameters
    ----------
    data : `~astropy.table.Table`, dict, `~numpy.ndarray`
        Astropy Table, dictionary of arrays or structured numpy array
        containing the "correct" column names.
    chromatic : bool
        True = magnification depends on wavelength and phase.
        False = magnification depends only on phase.
    magformat : str
        Format of the magnification column.  May be ``multiply`` or ``add,``
        where ``multiply`` means the magnification column provides a
        multiplicative magnification factor, mu, so the effect is applied to
        the source as flux * mu, and ``add`` means the magnification column
        provides an additive magnitude, DeltaM=-2.5*log10(mu).
    """

    def __init__(self, data, chromatic=False, magformat='multiply'):
        self.chromatic = chromatic

        # get column names in input data
        if isinstance(data, Table):
            colnames = data.colnames
        elif isinstance(data, np.ndarray):
            colnames = data.dtype.names
        elif isinstance(data, dict):
            colnames = data.keys()
        else:
            raise ValueError('unrecognized data type')

        if self.chromatic:
            mapping = alias_map(colnames, MLDATA_ALIASES,
                                required=CHROMATIC_MLDATA_REQUIRED_ALIASES)
        else:
            mapping = alias_map(colnames, MLDATA_ALIASES,
                                required=ACHROMATIC_MLDATA_REQUIRED_ALIASES)

        self.magnification = np.asarray(data[mapping['magnification']])
        magform = magformat[:2].lower()
        if magform not in ['ad','mu']:
            raise RuntimeError("``magformat`` must be ``multiply`` or ``add``")
        if magform=='ad':
            self.magnification = 10**(-0.4*self.magnification)
        self.time = np.asarray(data[mapping['time']])

        if self.chromatic:
            self.wavelength = np.asarray(data[mapping['wavelength']])

        # ensure columns are equal length
        if isinstance(data, dict):
            if not (len(self.time) == len(self.magnification)):
                raise ValueError("unequal column lengths")
            if self.chromatic:
                if not (len(self.time) == len(self.wavelength)):
                    raise ValueError("unequal column lengths")

    def sort_by_time(self):
        if not np.all(np.ediff1d(self.time) >= 0.0):
            idx = np.argsort(self.time)
            self.time = self.time[idx]
            self.magnification = self.magnification[idx]
            if self.chromatic:
                self.wavelength = self.wavelength[idx]

    def __len__(self):
        return len(self.time)

    def __getitem__(self, key):
        newdata = copy.copy(self)
        newdata.time = self.time[key]
        newdata.magnifciation = self.magnification[key]
        if self.chromatic:
            newdata.wavelength = self.wavelength[key]
        return newdata


    def magnification_interpolator(self):
        """Return an interpolation function that provides the microlensing
        magnification at any phase (and wavelength, if microlensing is
        chromatic).
        """
        
        if self.chromatic:
            return interp2d(self.time, self.wavelength, self.magnification,
                            bounds_error=False, fill_value=1.0, kind='cubic')
        return interp1d(self.time, self.magnification, bounds_error=False,
                        fill_value=1.0, kind='cubic')


def microlensing_data(data):
    if isinstance(data, MicrolensingData):
        return data
    else:
        return MicrolensingData(data)


def read_mldatafile(datafilename, magformat='multiply', **kwargs):
    """Read in microlensing data from a file.
    NAN values in the magnification array are converted to 1 if
    the magnification format is multiplicative, and 0 if additive.

    magformat : str
        Format of the magnification column.  May be ``multiply`` or ``add,``
        where ``multiply`` means the magnification column provides a
        multiplicative magnification factor, mu, so the effect is applied to
        the source as flux * mu, and ``add`` means the magnification column
        provides an additive magnitude, DeltaM=-2.5*log10(mu).

    """
    #TODO: parse header info for ``magformat`` and ``chromatic``

    datafilepath = path.abspath(path.expanduser(datafilename))
    ext = path.splitext(path.basename(datafilepath))[1].lower()
    if 'format' in kwargs:
        datatable = Table.read(datafilename, **kwargs)
    elif ext in ['.txt', '.text', '.dat']:
        datatable = Table.read(datafilename, format='ascii', **kwargs)
    elif ext in ['.fits']:
        datatable = Table.read(datafilename, format='fits', **kwargs)
    else:
        datatable = Table.read(datafilename, **kwargs)

    if 'col1' in datatable.colnames:
        datatable.rename_column('col1', 'time')
        if 'col3' in datatable.colnames:
            datatable.rename_column('col3', 'magnification')
            datatable.rename_column('col2', 'wavelength')
        else:
            datatable.rename_column('col2', 'magnification')

    if magformat.lower().startswith('mu'):
        inan = np.isnan(datatable['magnification'])
        datatable['magnification'][inan] = 1.0
    else:
        datatable['magnification'] = np.nan_to_num(datatable['magnification'])

    return MicrolensingData(datatable)
