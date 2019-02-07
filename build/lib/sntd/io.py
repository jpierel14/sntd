from collections import OrderedDict as odict

import numpy as np
import os,string,sncosmo,sys
from astropy.io import ascii
from astropy.table import Table,vstack,Column
from scipy.stats import mode
from copy import deepcopy,copy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
try:
    import pickle
except:
    import cPickle

from .util import *

__all__=['curve','curveDict','read_data','write_data','table_factory']

_comment_char={'#','='}
_meta__={'@','$','%','!','&'}

def _sntd_deepcopy(obj):
    newCurve=curveDict()
    for k in obj:
        setattr(newCurve,k,obj[k])
    return(newCurve)

class curveDict(dict):
    """The main object for SNTD. This organizes a MISN, containing
        the multiple light curves in self.images, all the fits, etc.
        See documentation for a flowchart of this object.
    """
    def __deepcopy__(self, memo):
        return deepcopy(dict(self))


    def __init__(self,telescopename="Unknown",object="Unknown"):
        """Constructor for curveDict class. Inherits from the dictionary class, and is the main object of organization used by SNTD.
        Parameters
        ----------
        telescopename : str
            Name of the telescope that the data were gathered from
        object : str
            Name of object of interest
        Returns
        -------
        MISN : `sntd.curveDict`
        """
        super(curveDict, self).__init__() #init for the super class
        self.meta = {'info': ''}
        """@type: dict
            @ivar: The metadata for the curveDict object, intialized with an empty "info" key value pair. It's
            populated by added _metachar__ characters into the header of your data file.
        """
        self.bands=set()
        """
        @type: list
        @ivar: The list of bands contained inside this curveDict
        """
        self.table=None
        """
        @type:~astropy.table.Table
        @ivar: The astropy table containing all of the data in your data file
        """
        self.telescopename=telescopename
        """
        @type: str
        @ivar: Name of the telescope that the data were gathered from
        """
        self.object=object
        """
        @type: str
        @ivar: Object of interest
        """
        self.images=dict([])

        self.combined=curve()
        self.color=curve()



    #these three functions allow you to access the curveDict via "dot" notation
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    __getattr__ = dict.__getitem__



    def __getstate__(self):
        """
        A function necessary for pickling
        :return: self
        """
        return self

    def __setstate__(self, d):
        """
        A function necessary for pickling
        :param d: A value
        :return: self.__dict__
        """
        self.__dict__ = d

    def __str__(self):
        """
        A replacement for the print method of the class, so that when you run print(curveDict()), this is how it shows
        up.
        """
        print('Telescope: %s'%self.telescopename)
        print('Object: %s'%self.object)
        print('Number of bands: %d' %len(self.bands))
        print('')
        for c in np.sort([x for x in self.images.keys()]):
            print('------------------')
            print('Image: %s:'%c)
            print('Bands: {}'.format(self.images[c].bands))
            print('Date Range: %.5f->%.5f' % (
            min(self.images[c].table[_get_default_prop_name('time')]),
            max(self.images[c].table[_get_default_prop_name('time')])))
            print('Number of points: %d' %len(self.images[c].table))
            if self.images[c].simMeta.keys():
                print('')
                print('Metadata:')
                print('\n'.join('   {}:{}'.format(*t) for t in zip(self.images[c].simMeta.keys(),self.images[c].simMeta.values()) if isinstance(t[1],(str,float,int))))
        return '------------------'

    def add_curve(self,myCurve):
        """Adds a curve object to the existing curveDict (i.e. adds
        an image to a MISN)"""
        self.bands.update([x for x in myCurve.bands if x not in self.bands])
        myCurve.object='image_'+str(len(self.images)+1)

        myCurve.table.add_column(Column([myCurve.object for i in range(len(myCurve.table))],name='image'))




        tempCurve=_sntd_deepcopy(myCurve)


        self.images[myCurve.object]=tempCurve

        if self.table:
            for row in tempCurve.table:
                self.table.add_row(row)
        else:
            self.table=copy(tempCurve.table)
        return(self)

    def combine_curves(self,time_delays=None,magnifications=None,referenceImage='image_1'):
        """Takes the multiple images in self.images and combines
            the data into a single light curve using defined
            time delays and magnifications or best (quick) guesses.
        Parameters
        ----------
        time_delays : dict
            Dictionary with image names as keys and relative time
            delays as values (e.g. {'image_1':0,'image_2':20})
        magnifications : dict
            Dictionary with image names as keys and relative
            magnifications as values (e.g.
            {'image_1':1,'image_2':1.1})
        """
        if len(self.images) <2:
            print("Not enough curves to combine!")
            return(self)
        if not time_delays:
            time_delays=_guess_time_delays(self,referenceImage) #TODO fix these guessing functions
        if not magnifications:
            magnifications=_guess_magnifications(self,referenceImage)
        self.combined.table=Table(names=self.table.colnames,dtype=[self.table.dtype[x] for x in self.table.colnames])
        for k in np.sort(list(self.images.keys())):
            temp=deepcopy(self.images[k].table)
            temp['time']-=time_delays[k]
            temp['flux']/=magnifications[k]
            temp.meta=dict([])

            self.combined.table=vstack([self.combined.table,temp])

        self.combined.table.sort('time')
        self.combined.bands=self.bands
        self.combined.meta['td']={k:float(time_delays[k]) for k in time_delays.keys()}
        self.combined.meta['mu']={k:float(magnifications[k]) for k in magnifications.keys()}
        return(self)


    def color_table(self,band1,band2,time_delays=None,magnifications=None,image=[]):
        image=list(image) if not isinstance(image,(list,tuple)) else image
        names=['time','image','zpsys']
        dtype=[self.table.dtype[x] for x in names]
        names=np.append(names,[band1+'-'+band2,band1+'-'+band2+'_err'])
        dtype=np.append(dtype,[dtype[0],dtype[0]])
        self.color.table=Table(names=names,dtype=dtype)
        if not time_delays:
            time_delays=_guess_time_delays(self) #TODO fix these guessing functions
        if not magnifications:
            magnifications=_guess_magnifications(self)
        for im in [x for x in self.images.keys() if x not in image]:
            temp2=copy(self.images[im].table[self.images[im].table['band']==band2])
            temp1=copy(self.images[im].table[self.images[im].table['band']==band1])
            temp1=temp1[temp1['flux']>0]
            temp2=temp2[temp2['flux']>0]
            temp1['flux']/=magnifications[im]
            temp2['flux']/=magnifications[im]
            temp1['time']-=time_delays[im]
            temp2['time']-=time_delays[im]


            temp2['mag']=-2.5*np.log10(temp2['flux'])+temp2['zp']
            temp2['magerr']=1.0857*temp2['fluxerr']/temp2['flux']
            temp1['mag']=-2.5*np.log10(temp1['flux'])+temp1['zp']
            temp1['magerr']=1.0857*temp1['fluxerr']/temp1['flux']
            temp1_remove=[i for i in range(len(temp1)) if temp1['time'][i] not in temp2['time']]
            temp1.remove_rows(temp1_remove)
            temp2_remove=[i for i in range(len(temp2)) if temp2['time'][i] not in temp1['time']]
            temp2.remove_rows(temp2_remove)

            temp1['magerr']=np.sqrt(temp2['magerr']**2+temp1['magerr']**2)


            temp1['mag']-=temp2['mag']
            temp1.rename_column('mag',band1+'-'+band2)
            temp1.rename_column('magerr',band1+'-'+band2+'_err')
            to_remove=[x for x in temp1.colnames if x not in names]
            temp1.remove_columns(to_remove)


            self.color.table=vstack([self.color.table,copy(temp1)])
        self.color.table.sort('time')
        return(self)


    def plot_object(self, bands='all', savefig=False,
                    filename='mySN', orientation='horizontal',method='separate',
                    showModel=False,showFit=False,showMicro=False):
        """Plot the multiply-imaged SN light curves and show/save to a file.
            Each subplot shows a single-band light curve, for all images of the SN.
        Parameters
        ----------
        bands : str or list of str
            'all' = plot all bands; or provide a list of bands to plot
        savefig : bool
            boolean to save or not save plot
        filename : str
            if savefig is True, this is the output filename
        orientation : str
            'horizontal' = all subplots are in a single row
            'vertical' = all subplots are in a single column
        method : str
            Plots the result of separate, combined, or color curve method
        showModel : bool
            If true, the underlying model before microlensing is plotted
            as well
        showFit : bool
            If true and it exists, the best fit model from
            self.images['image'].fits.model is plotted
        showMicro : bool
            If true and it exists, the simulated microlensing is plotted
            as well
        Returns
        -------
        figure : `~matplotlib.pyplot.figure`
        """




        colors=['r','g','b','k','m']
        i=0
        leg=[]

        if method=='combined':
            if bands == 'all':
                bands = set(self.combined.table['band'])

            nbands = len(bands)
            if orientation.startswith('v'):
                ncols = 1
                nrows = nbands
            else:
                ncols = nbands
                nrows = 1
            fig,axlist=plt.subplots(nrows=nrows, ncols=ncols,
                                    sharex=True, sharey=False,figsize=(10,10))
            if nbands==1:
                axlist = [axlist]
            for lc in np.sort([x for x in self.images.keys()]):
                temp=self.combined.table[self.combined.table['image']==lc]
                for b, ax in zip(bands, axlist):
                    if b==list(bands)[0]:
                        leg.append(
                            ax.errorbar(temp['time'][
                                            temp['band']==b],
                                        temp['flux'][
                                            temp['band']==b],
                                        yerr=temp['fluxerr'][
                                            temp['band']==b],
                                        markersize=4, fmt=colors[i]+'.'))
                    else:
                        ax.errorbar(temp['time'][
                                        temp['band']==b],
                                    temp['flux'][
                                        temp['band']==b],
                                    yerr=temp['fluxerr'][
                                        temp['band']==b],
                                    markersize=4, fmt=colors[i]+'.')
                    if showFit:
                        ax.plot(np.arange(np.min(temp['time'][temp['band']==b]),np.max(temp['time'][temp['band']==b]),1),
                                self.combined.fits.model.bandflux(b,np.arange(np.min(temp['time'][temp['band']==b]),np.max(temp['time'][temp['band']==b]),1),
                                                                  zp=temp['zp'][temp['band']==b][0],
                                                                  zpsys=temp['zpsys'][temp['band']==b][0]),color='y')

                    ax.text(0.95, 0.95, b.upper(), fontsize='large',
                            transform=ax.transAxes, ha='right', va='top')


                i+=1
        elif method =='color':
            if bands=='all':
                if len([x for x in self.color.table.colnames if '-' in x and '_' not in x])!=1:
                    print("Want to plot color curves but need 2 bands specified.")
                    sys.exit(1)
                else:
                    colname=[x for x in self.color.table.colnames if '-' in x and '_' not in x][0]
                    bands=[colname[:colname.find('-')],colname[colname.find('-')+1:]]
            elif len(bands) !=2:
                print("Want to plot color curves but need 2 bands specified.")
                sys.exit(1)


            fig=plt.figure(figsize=(10,10))
            ax=fig.gca()
            for lc in np.sort([x for x in self.images.keys()]):
                temp=self.color.table[self.color.table['image']==lc]


                ax.errorbar(temp['time'],temp[bands[0]+'-'+bands[1]],yerr=temp[bands[0]+'-'+bands[1]+'_err'],
                            markersize=4, fmt=colors[i]+'.')

                ax.text(0.95, 0.95, bands[0].upper()+'-'+bands[1].upper(), fontsize='large',
                        transform=ax.transAxes, ha='right', va='top')

                i+=1
            if showFit:
                mod_time=np.arange(np.min(self.color.table['time']),np.max(self.color.table['time']),1)
                modCol=self.color.fits.model.color(bands[0],bands[1],self.color.table['zpsys'][0],mod_time)

                ax.plot(mod_time,modCol,color='y')
        else:
            if bands == 'all':
                bands = self.bands

            nbands = len(bands)
            if orientation.startswith('v'):
                ncols = 1
                nrows = nbands
            else:
                ncols = nbands
                nrows = 1
            fig,axlist=plt.subplots(nrows=nrows, ncols=ncols,
                                    sharex=True, sharey=True,figsize=(10,10))
            if nbands==1:
                axlist = [axlist]
            microAx={b:False for b in bands}
            for lc in np.sort([x for x in self.images.keys()]):
                for b, ax in zip(bands, axlist):
                    if b==list(bands)[0]:
                        leg.append(
                            ax.errorbar(self.images[lc].table['time'][
                                            self.images[lc].table['band']==b],
                                        self.images[lc].table['flux'][
                                            self.images[lc].table['band']==b],
                                        yerr=self.images[lc].table['fluxerr'][
                                            self.images[lc].table['band']==b],
                                        markersize=4, fmt=colors[i]+'.'))
                        if showMicro:
                            ax.set_ylabel('Flux',fontsize='large')
                    else:
                        ax.errorbar(self.images[lc].table['time'][
                                        self.images[lc].table['band']==b],
                                    self.images[lc].table['flux'][
                                        self.images[lc].table['band']==b],
                                    yerr=self.images[lc].table['fluxerr'][
                                        self.images[lc].table['band']==b],
                                    markersize=4, fmt=colors[i]+'.')

                    ax.text(0.95, 0.95, b.upper(), fontsize='large',
                            transform=ax.transAxes, ha='right', va='top')

                    if showFit:
                        time_model = np.arange(self.images[lc].table['time'].min(),
                                               self.images[lc].table['time'].max(),
                                               0.1)
                        ax.plot(time_model,self.images[lc].fits.model.bandflux(b,time_model,
                                    np.mean(self.images[lc].table['zp'][self.images[lc].table['band']==b]),
                                    self.images[lc].zpsys),'k-')
                    if showMicro:
                        time_model=np.arange(self.images[lc].table['time'].min(),
                                             self.images[lc].table['time'].max(),
                                             0.1)
                        if microAx[b] is False:
                            ax_divider = make_axes_locatable(ax)
                            ax_ml = ax_divider.append_axes("bottom", size="25%", pad=.4)
                            #ax_ml.set_xlabel('Days (Observer Frame)',fontsize='large')
                            if b==list(bands)[0]:
                                ax_ml.set_ylabel('Microlensing ($\mu$)',fontsize='large')
                            microAx[b]=ax_ml
                        else:
                            ax_ml=microAx[b]
                        ax_ml.plot(time_model,self.images[lc].simMeta['microlensing_params'](time_model/(1+self.images[lc].simMeta['sourcez'])),color=colors[i],linewidth=3)

                    if showModel:
                        # Plot the underlying model, including dust and lensing
                        # effects other than microlesning, as a black curve for
                        # each simulated SN image
                        time_model = np.arange(self.images[lc].table['time'].min(),
                                               self.images[lc].table['time'].max(),
                                               0.1)
                        time_shifted = time_model - self.images[lc].simMeta['td']


                        flux_magnified = self.model.bandflux(
                            b, time_shifted, self.images[lc].table['zp'][self.images[lc].table['band']==b][0],
                            self.images[lc].zpsys) * \
                                         self.images[lc].simMeta['mu']
                        ax.plot(time_model, flux_magnified, 'k-')


                i+=1


        plt.figlegend(leg,np.sort([x for x in self.images.keys()]), frameon=False,
                      loc='center right', fontsize='medium', numpoints=1)


        if not showMicro:
            fig.text(0.02, .5, 'Flux', va='center',
                rotation='vertical',fontsize='large')
        fig.text(0.5, 0.02, r'Observer-frame time (days)', ha='center',
                 fontsize='large')

        plt.suptitle('Multiply-Imaged SN "'+self.object+'"--'+self.telescopename,fontsize=16)
        if savefig:
            plt.savefig(filename+'.pdf',format='pdf',overwrite=True)
        return fig




class curve(dict):
    """SNTD class that describes each image of a MISN
    """
    def __deepcopy__(self, memo):
        return deepcopy(dict(self))

    def __init__(self,zpsys='AB'):
        #todo: implement more general read and write functions
        """
        Constructor for curve class,
        :param band, zp, zpsys
        """
        #todo populate param documentation
        super(curve,self).__init__()
        self.meta = {'info': ''}
        """@type: dict
            @ivar: The metadata for the curveDict object, intialized with an empty "info" key value pair. It's
            populated by added _metachar__ characters into the header of your data file.
        """
        self.table=None
        """@type: ~astropy.table.Table
        @ivar: A table containing the data, used for SNCosmo functions, etc.
        """
        self.bands=[]
        """@type: string
        @ivar: band names, used
        """
        self.zpsys=zpsys

        self.simMeta=dict([])

        self.fits=None


    #these three functions allow you to access the curveDict via "dot" notation
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    __getattr__ = dict.__getitem__



    def __getstate__(self):
        """
        A function necessary for pickling
        :return: self
        """
        return self

    def __setstate__(self, d):
        """
        A function necessary for pickling
        :param d: A value
        :return: self.__dict__
        """
        self.__dict__ = d

    def validate(self, verbose=False):
        """
        Simply extends the superclass version of validate to include flux
        :param verbose: Default=False -> won't print "Validation done!"

        :return: None
        """
        ndates=len(self)
        if len(self.fluxes) != ndates or len(self.fluxerrs) != ndates or len(self.table)!=ndates:
            raise(RuntimeError, "The length of your flux/fluxerr arrays are inconsistent with rest of curve!")

        super(curve,self).validate(verbose)



def table_factory(table,telescopename="Unknown",object=None):
    """This function will create a new curve object using an astropy table.
    Parameters
    ----------
    table : `~astropy.table.Table`
        Astropy table with all of your data from your data file.
    telescopename : str
        Name of telescope for labeling purposes inside curve object
    object : str
        Name of object for labeling purposes inside curve object
        (e.g. SN2006jf, etc.)
    Returns
    -------
    curve : `~sntd.curve`
    """
    newlc=curve()


    table=standardize_table_colnames(table)
    newlc.bands={x for x in table[_get_default_prop_name('band')]}
    zp_dict=dict([])
    zpsys_dict=dict([])
    for band in newlc.bands:
        zp=set(table[_get_default_prop_name('zp')][table[_get_default_prop_name('band')]==band])
        zpsys=set(table[_get_default_prop_name('zpsys')][table[_get_default_prop_name('band')]==band])
        zp_dict[band]=zp.pop() if len(zp)==1 else np.array(zp)
        zpsys_dict[band]=zpsys.pop() if len(zpsys)==1 else np.array(zpsys)


    newlc.table=table


    newlc.telescopename = telescopename
    newlc.object = object


    return newlc


def _switch(ext):
    switcher = {
        '.pkl': _read_pickle
    }
    return switcher.get(ext, _read_data)


def write_data(curves,filename=None,protocol=-1):
    """Used to write a curveDict object to a pickle
        to be read later
    Parameters
    ----------
    curves : `~sntd.curveDict`
    filename : str
        Name of output file
    protocol : int
        Pickling protocol
    """
    if not filename:
        filename=curves.object
    with open(filename,'wb') as handle:
        try:
            pickle.dump(curves,handle,protocol=protocol)
        except:
            cPickle.dump(curves,handle,protocol=protocol)
    return




def read_data(filename,**kwargs):
    """Used to read a light curve or curve object in pickle format.
        Either way, it'll come out as a curve object.
    Parameters
    ----------
    filename : str
        Name of the file to be read (ascii or pickle)
    Returns
    -------
    curve : `~sntd.curve` or `~sntd.curveDict`
    """
    return(_switch(os.path.splitext(filename)[1])(filename,**kwargs))


def _read_pickle(filename,telescopename="Unknown",object="Unknown",**kwargs):
    try:
        return (pickle.load(open(filename,'rb')))
    except:
        return (cPickle.load(open(filename,'rb')))



def standardize_table_colnames(table):
    for col in table.colnames:
        if col != _get_default_prop_name(col.lower()):
            table.rename_column(col, _get_default_prop_name(col.lower()))
    return table

def _read_data(filename,**kwargs):
    table_done=True
    try:
        table = sncosmo.read_lc(filename, masked=True,**kwargs)
        for col in table.colnames:
            if col != _get_default_prop_name(col.lower()):
                table.rename_column(col, _get_default_prop_name(col.lower()))
    except:
        try:
            table = ascii.read(filename, masked=True,**kwargs)
            for col in table.colnames:
                if col != _get_default_prop_name(col.lower()):
                    table.rename_column(col, _get_default_prop_name(col.lower()))
        except:
            table = Table(masked=True)
            table_done=False


    delim = kwargs.get('delim', None)
    myCurve=curve()
    with anyOpen(filename) as f:
        lines=f.readlines()
        length=mode([len(l.split()) for l in lines])[0][0]#uses the most common line length as the correct length
        for i,line in enumerate(lines):
            if np.any([x in line for x in _comment_char]):
                continue
            line = line.strip(delim)

            if len(line)==0:
                continue
            if np.any([x in line for x in _meta__]):
                pos = line.find(' ')
                if (lines[i + 1][0] in _meta__ or lines[i + 1][0] in _comment_char or len(
                        line) != length):  # just making sure we're not about to put the col names into meta
                    if (pos == -1 or not any([_isfloat(x) for x in line.split()])):
                        if line[-1] not in string.punctuation:
                            line=line+'.'
                    else:
                        myCurve.meta[line[1:pos]] = _cast_str(line[pos:])
                continue
            line=line.split()
            if len(line)!= length:
                raise(RuntimeError,"Make sure your data are in a square matrix.")
            if table.colnames:
                colnames=table.colnames
            else:
                colnames = odict.fromkeys([_get_default_prop_name(x.lower()) for x in line])
                if len(colnames) != len(line):
                    raise (RuntimeError, "Do you have duplicate column names?")
                colnames=odict(zip(colnames,range(len(colnames))))
            startLine=i
            break
    f.close()
    if not table_done:
        lines = [x.strip().split() for x in lines[startLine + 1:] if not np.any([y in x for y in _comment_char])]

        for col in colnames:
            table[col]=np.asarray([_cast_str(x[colnames[col]]) for x in lines])
        colnames=colnames.keys()

    for col in [_get_default_prop_name(x) for x in ['band','zp','zpsys']]:
        if col not in colnames:
            temp=kwargs.get(col,None)
            if temp is not None:
                table[col]=temp
            else:
                print('Column "%s" is not in your file and you did not define it in kwargs.'%col)
                sys.exit(1)
    table=standardize_table_colnames(table)
    bnds = {x for x in table[_get_default_prop_name('band')]}
    table=_norm_flux_mag(table)
    for band in bnds:
        if _isfloat(band[0]):
            band='band_'+band
        try:
            if band[0:5]=='band_':
                sncosmo.get_bandpass(band[5:])
            else:
                sncosmo.get_bandpass(band)
        except:
            print('Skipping band %s, not in registry.' %band)
            table.mask[table[_get_default_prop_name('band')]==band]=True
            continue
        myCurve.bands.append(band)

    myCurve.table=table
    myCurve.zpsys=table['zpsys'][0]
    return myCurve

def _norm_flux_mag(table):
    if 'mag' not in table.colnames:
        if 'flux' not in table.colnames:
            print('No Mag or Flux data.')
            sys.exit(1)
        table=_flux_to_mag(table)
    else:
        if 'flux' not in table.colnames:
            table=_mag_to_flux(table)
    return table


def _flux_to_mag(table):
    table[_get_default_prop_name('mag')] = np.asarray(map(lambda x, y: -2.5 * np.log10(x) + y, table[_get_default_prop_name('flux')],
                                  table[_get_default_prop_name('zp')]))
    table[_get_default_prop_name('magerr')] = np.asarray(map(lambda x, y: 2.5 * np.log10(np.e) * y / x, table[_get_default_prop_name('flux')],
                                     table[_get_default_prop_name('fluxerr')]))
    table[_get_default_prop_name('magerr')][np.isnan(table[_get_default_prop_name('mag')])]=np.nan
    return table


def _mag_to_flux(table):
    table[_get_default_prop_name('flux')] = np.asarray(
        map(lambda x, y: 10 ** (-.4 * (x -y)), table[_get_default_prop_name('mag')],
            table[_get_default_prop_name('zp')]))
    table[_get_default_prop_name('fluxerr')] = np.asarray(
        map(lambda x, y: x * y / (2.5 * np.log10(np.e)), table[_get_default_prop_name('magerr')],
            table[_get_default_prop_name('flux')]))
    return table



