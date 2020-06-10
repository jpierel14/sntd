from collections import OrderedDict as odict

import numpy as np
import os,string,sncosmo,sys,corner
from astropy.io import ascii
from astropy.table import Table,vstack,Column
from scipy.stats import mode
from copy import deepcopy,copy
import matplotlib.pyplot as plt
from sncosmo.snanaio import read_snana_fits
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

class curve(dict):
    """
    SNTD class that describes each image of a MISN
    """
    def __deepcopy__(self, memo):
        return deepcopy(dict(self))

    def __init__(self,zpsys='AB'):
        """
        Constructor for curve class


        """
        #todo populate param documentation
        super(curve,self).__init__()
        self.meta = {'info': ''}
        """@type: :class:`dict`
            The metadata for the curveDict object, intialized with an empty "info" key value pair. It's
            populated by added _metachar__ characters into the header of your data file.
        """
        self.table=None
        """@type: :class:`astropy.table.Table`
            A table containing the data, used for SNCosmo functions, etc.
        """
        self.bands=[]
        """@type: str
            band names used
        """
        self.zpsys=zpsys
        """@type: str
            The zero-point system for this curve object
        """

        self.simMeta=dict([])
        """@type: :class:`dict`
            A dictionary containing simulation metadata if this is a simulated curve object"""

        self.fits=None
        """@type: :class:`~sntd.fitting.newDict`
            Contains fit information from fit_data"""


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


class curveDict(dict):
    """
    The main object for SNTD. This organizes a MISN, containing
    the multiple light curves in self.images, all the fits, etc.
    """
    def __deepcopy__(self, memo):
        return deepcopy(dict(self))


    def __init__(self,telescopename="Unknown",object="Unknown"):
        """
        Constructor for curveDict class. Inherits from the dictionary class,
        and is the main object of organization used by SNTD.

        Parameters
        ----------
        telescopename : str
            Name of the telescope that the data were gathered from
        object : str
            Name of object of interest

        Returns
        -------
        MISN : :class:`~sntd.curveDict`
        """
        super(curveDict, self).__init__() #init for the super class
        self.meta = {'info': ''}
        """@type: :class:`dict`
            The metadata for the curveDict object, intialized with an empty "info" key value pair. It's
            populated by added _metachar__ characters into the header of your data file.
        """
        self.bands=set()
        """
        @type: :class:`list`
            The list of bands contained inside this curveDict
        """
        self.table=None
        """
        @type: :class:`~astropy.table.Table`
            The astropy table containing all of the data in your data file
        """
        self.telescopename=telescopename
        """
        @type: str
            Name of the telescope that the data were gathered from
        """
        self.object=object
        """
        @type: str
            Object of interest
        """
        self.images=dict([])

        self.parallel=curve()
        self.series=curve()
        self.color=curve()
        self.constants={}


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
            min(self.images[c].table[get_default_prop_name('time')]),
            max(self.images[c].table[get_default_prop_name('time')])))
            print('Number of points: %d' %len(self.images[c].table))
            if self.images[c].simMeta.keys():
                print('')
                print('Metadata:')
                print('\n'.join('   {}:{}'.format(*t) for t in zip(self.images[c].simMeta.keys(),self.images[c].simMeta.values()) if isinstance(t[1],(str,float,int))))
        return '------------------'


    def add_curve(self,myCurve,key=None):
        """Adds a curve object to the existing curveDict (i.e. adds
        an image to a MISN)

        Parameters
        ----------
        myCurve: :class:`sntd.curve`
            The curve to add to self.
        key: str
            The key you want to save this as, default is 'image_1,image_2,etc.'

        Returns
        -------
        self: :class:`sntd.curve_io.curveDict`
        """

        self.bands.update([x for x in myCurve.bands if x not in self.bands])
        if key is None:
            myCurve.object='image_'+str(len(self.images)+1)
        else:
            myCurve.object=key

        if 'image' not in myCurve.table.colnames:
            myCurve.table.add_column(Column([myCurve.object for i in range(len(myCurve.table))],name='image'))


        if 'zpsys' not in myCurve.table.colnames:
            print('Assuming AB magsystem...')
            myCurve.table['zpsys']='AB'
        else:
            myCurve.zpsys=myCurve.table['zpsys'][0]

        tempCurve=_sntd_deepcopy(myCurve)


        self.images[myCurve.object]=tempCurve

        if self.table:
            self.table=vstack([tempCurve.table])
        else:
            self.table=copy(tempCurve.table)
        return(self)

    def combine_curves(self,time_delays=None,magnifications=None,referenceImage='image_1',static=False,
                       model=None,minsnr=0):
        """
        Takes the multiple images in self.images and combines
        the data into a single light curve using defined
        time delays and magnifications or best (quick) guesses.

        Parameters
        ----------
        time_delays: :class:`dict`
            Dictionary with image names as keys and relative time
            delays as values (e.g. {'image_1':0,'image_2':20}). Guessed if None.
        magnifications: :class:`dict`
            Dictionary with image names as keys and relative magnifications
            as values (e.g. {'image_1':0,'image_2':20}). Guessed if None.
        referenceImage: str
            The image you want to be the reference (e.g. image_1, image_2, etc.)
        ignore_images: :class:`~list`
            List of images you do not want to include in the color curve.
        static: bool
            Make the color curve, don't shift the data
        model: :class:`~sncosmo.Model`
            If you want to use an sncosmo Model (and the guess_t0_amplitude method) to guess time delays
        minsnr: float
            Cut data that don't meet this threshold before making the color curve.

        Returns
        -------
        self: :class:`sntd.curve_io.curveDict`
        """
        if len(self.images) <2:
            print("Not enough curves to combine!")
            return(self)

        if time_delays is None:
            if model is not None:
                time_delays={}
                magnifications={}
                model=sncosmo.Model(model) if isinstance(model,str) else model
                ref_t0,ref_amp=sncosmo.fitting.guess_t0_and_amplitude(sncosmo.photdata.photometric_data( \
                    self.images[referenceImage].table),model,minsnr)
                self.series.meta['reft0']=ref_t0
                self.series.meta['refamp']=ref_amp
                time_delays[referenceImage]=0
                magnifications[referenceImage]=1
                for k in self.images.keys():
                    if k==referenceImage:
                        continue
                    guess_t0,guess_amp=sncosmo.fitting.guess_t0_and_amplitude(sncosmo.photdata.photometric_data(\
                        self.images[k].table),model,minsnr)
                    time_delays[k]=guess_t0-ref_t0
                    magnifications[k]=guess_amp/ref_amp
            else:
                time_delays=guess_time_delays(self,referenceImage) #TODO fix these guessing functions
        if magnifications is None:
            magnifications=guess_magnifications(self,referenceImage)

        
        self.series.table=Table(names=self.table.colnames,dtype=[self.table.dtype[x] for x in self.table.colnames])
        for k in np.sort(list(self.images.keys())):
            temp=deepcopy(self.images[k].table)
            if not static:
                temp['time']-=time_delays[k]
                temp['flux']/=magnifications[k]
            temp.meta=dict([])

            self.series.table=vstack([self.series.table,temp])

        self.series.table.sort('time')
        self.series.bands=self.bands
        self.series.meta['td']={k:float(time_delays[k]) for k in time_delays.keys()}
        self.series.meta['mu']={k:float(magnifications[k]) for k in magnifications.keys()}


        return(self)


    def color_table(self,band1s,band2s,time_delays=None,referenceImage='image_1',ignore_images=[],
                    static=False,model=None,minsnr=0.0):
        """
        Takes the multiple images in self.images and combines
        the data into a single color curve using defined
        time delays and magnifications or best (quick) guesses.

        Parameters
        ----------
        band1s: str or list
            The first band(s) for color curve(s)
        band2s: str or list
            The second band(s) for color curve(s)
        time_delays: :class:`dict`
            Dictionary with image names as keys and relative time
            delays as values (e.g. {'image_1':0,'image_2':20}). Guessed if None.
        referenceImage: str
            The image you want to be the reference (e.g. image_1, image_2, etc.)
        ignore_images: :class:`~list`
            List of images you do not want to include in the color curve.
        static: bool
            Make the color curve, don't shift the data
        model: :class:`~sncosmo.Model`
            If you want to use an sncosmo Model (and the guess_t0_amplitude method) to guess time delays
        minsnr: float
            Cut data that don't meet this threshold before making the color curve.

        Returns
        -------
        self: :class:`sntd.curve_io.curveDict`
        """

        ignore_images=list(ignore_images) if not isinstance(ignore_images,(list,tuple)) else ignore_images
        names=['time','image','zpsys']
        dtype=[self.table.dtype[x] for x in names]
        names=np.append(names,np.append(np.array([[band1+'-'+band2,band1+'-'+band2+'_err'] for band1,band2 in zip(band1s,band2s)]).flatten(),
                                    np.unique([['flux_%s'%band1,'fluxerr_%s'%band1,'flux_%s'%band2,'fluxerr_%s'%band2,'zp_%s'%band1,'zp_%s'%band2]\
                                     for band1,band2 in zip(band1s,band2s)]).flatten()))
        dtype=np.append(dtype,[dtype[0]]*(len(names)-len(dtype)))
        self.color.table=Table(names=names,dtype=dtype)
        

        if time_delays is None:
            if model is not None:
                time_delays={}
                model=sncosmo.Model(model) if isinstance(model,str) else model
                ref_t0,ref_amp=sncosmo.fitting.guess_t0_and_amplitude(sncosmo.photdata.photometric_data( \
                    self.images[referenceImage].table),model,minsnr)
                self.color.meta['reft0']=ref_t0
                time_delays[referenceImage]=0
                for k in self.images.keys():
                    if k==referenceImage:
                        continue
                    guess_t0,guess_amp=sncosmo.fitting.guess_t0_and_amplitude(sncosmo.photdata.photometric_data( \
                        self.images[k].table),model,minsnr)
                    time_delays[k]=guess_t0-ref_t0
            else:
                time_delays=guess_time_delays(self,referenceImage) #TODO fix these guessing functions

        
        self.color.meta['td']=time_delays
        for im in [x for x in self.images.keys() if x not in ignore_images]:
            
            for band1,band2 in zip(band1s,band2s):
                to_add={}
                temp2=deepcopy(self.images[im].table[self.images[im].table['band']==band2])
                temp1=deepcopy(self.images[im].table[self.images[im].table['band']==band1])
              
              
                temp1=temp1[temp1['flux']>0]
                temp2=temp2[temp2['flux']>0]
                temp1=temp1[temp1['flux']/temp1['fluxerr']>minsnr]
                temp2=temp2[temp2['flux']/temp2['fluxerr']>minsnr]
                if not static:
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

                
                to_add['time']=temp1['time']
                to_add['image']=[im]*len(temp1)
                to_add['zpsys']=temp1['zpsys']
                to_add[band1+'-'+band2]=temp1['mag']
                to_add[band1+'-'+band2+'_err']=temp1['magerr']
                to_add['flux_%s'%band1]=temp1['flux']
                to_add['fluxerr_%s'%band1]=temp1['fluxerr']
                to_add['flux_%s'%band2]=temp2['flux']
                to_add['fluxerr_%s'%band2]=temp2['fluxerr']
                to_add['zp_%s'%band1]=temp1['zp']
                to_add['zp_%s'%band2]=temp2['zp']
                for col in [x for x in names if x not in to_add.keys()]:
                    to_add[col]=[np.nan]*len(temp1)
                
                


                for i in range(len(temp1)):
                    self.color.table.add_row({k:to_add[k][i] for k in to_add.keys()})
        self.color.table.meta={}

        self.color.table.sort('time')

        return(self)

    def clip_data(self,im,minsnr=0,mintime=-np.inf,maxtime=np.inf,peak=0,remove_bands=[],max_cadence=None):
        """
        Clips the data of an image based on various properties.

        Parameters
        ----------
        im: str
            The image to clip
        minsnr: float
            Clip based on a minimum SNR
        mintime: float
            Clip based on a minimum time (observer frame relative to peak)
        maxtime: float
            Clip based on a maximum time (observer frame relative to peak)
        peak: float
            Used in conjunction with min/max time
        remove_bands: list
            List of bands to remove from the light curve
        max_cadence: float
            Clips data so that points are spread by at least max_cadence
        """
        
        self.images[im].table=self.images[im].table[self.images[im].table['flux']/\
                                                    self.images[im].table['fluxerr']>minsnr]
        
        self.images[im].table=self.images[im].table[self.images[im].table['time']>mintime+peak]
        self.images[im].table=self.images[im].table[self.images[im].table['time']<maxtime+peak]
        

        for b in remove_bands:
            self.images[im].table=self.images[im].table[self.images[im].table['band']!=b]

        if max_cadence is not None and isinstance(max_cadence,(int,float)):
            to_remove=[]
            for b in np.unique(self.images[im].table['band']):
                binds=np.where(self.images[im].table['band']==b)[0]
                t=self.images[im].table['time'][binds[0]]
                for i in range(1,len(binds)):
                    if self.images[im].table[binds[i]]['time']<t+max_cadence:
                        to_remove.append(binds[i])
                    else:
                        t=self.images[im].table[binds[i]]['time']
            self.images[im].table.remove_rows(to_remove)

    def quality_check(self,min_n_bands=1,min_n_points_per_band=1,clip=False,method='parallel'):
        """
        Checks the images of a SN to make sure they pass minimum thresholds for fitting.

        Parameters
        ----------
        min_n_bands: int
            The minimum number of bands needed to pass
        min_n_points_per_band: int
            The minimum number of bands in a given band to pass
        clip: bool
            If True, "bad" bands are clipped in place
        method: str
            Should be parallel, series, or color. Checks all images (parallel), or the series
            table (series), or the color table (color)
        Returns
        -------
        self: :class:`sntd.curve_io.curveDict`
        """
        if method=='parallel':
            good_bands=[]
            for im in self.images.keys():
                ngood_bands=0
                for b in np.unique(self.images[im].table['band']):
                    temp_n_for_b=len(self.images[im].table[self.images[im].table['band']==b])
                    if temp_n_for_b<min_n_points_per_band:
                        if clip:
                            self.images[im].table=self.images[im].table[self.images[im].table['band']!=b]
                    else:
                        ngood_bands+=1
                        good_bands.append(b)
                if ngood_bands<min_n_bands:
                    return False
            self.bands=np.unique(good_bands)
        elif method=='series':
            ngood_bands=0
            for b in np.unique(self.series.table['band']):
                temp_n_for_b=len(self.series.table[self.series.table['band']==b])
                if temp_n_for_b<min_n_points_per_band:
                    if clip:
                        self.series.table=self.series.table[self.series.table['band']!=b]
                else:
                    ngood_bands+=1
            if ngood_bands<min_n_bands:
                return False
        elif method=='color':
            if len(self.color.table)<min_n_points_per_band:
                return False
        else:
            print('method unknown for quality_check')
            sys.exit(1)

        return True

    def plot_fit(self,method='parallel',par_image=None):
        """
        Makes a corner plot based on one of the fitting methods

        Parameters
        ----------
        method: str
            parallel, series, or color

        Returns
        -------
        figure object: :class:`~matplotlib.pyplot.figure`
        """
        if method=='parallel':
            if par_image is None:
                par_image=self.parallel.fitOrder[0]
            res=self.images[par_image].fits.res
            samples=res.samples
            try:
                truths=[self.images[par_image].simMeta['model'].get(x) for x in res.vparam_names]
            except:
                truths=None
        elif method=='series':
            res=self.series.fits.res
            samples=res.samples
            
            try:
                truths=[]
                for p in res.vparam_names:
                    if p.startswith('mu_'):
                        im=[x for x in self.images.keys() if x[-1]==p[-1]][0]
                        truths.append(self.images[im].simMeta['model'].parameters[2]/ \
                                      self.images[self.series.refImage].simMeta['model'].parameters[2])

                    elif p.startswith('dt_'):
                        im=[x for x in self.images.keys() if x[-1]==p[-1]][0]
                        truths.append(self.images[im].simMeta['model'].get('t0')-\
                            self.images[self.series.refImage].simMeta['model'].get('t0'))
                    else:
                        im=self.series.refImage
                        truths.append(self.images[im].simMeta['model'].get(p))
            except:
                truths=None

        else:
            res=self.color.fits.res
            samples=res.samples

            try:
                truths=[]
                for p in res.vparam_names:
                    if p.startswith('dt_'):
                        im=[x for x in self.images.keys() if x[-1]==p[-1]][0]
                        truths.append(self.images[im].simMeta['model'].get('t0')-\
                            self.images[self.color.refImage].simMeta['model'].get('t0'))
                    else:
                        im=list(self.images.keys())[0]
                        truths.append(self.images[self.color.refImage].simMeta['model'].get(p))
            except:
                truths=None
        fig = corner.corner(
            samples,
            weights=res.weights,
            labels=res.vparam_names,
            truths=truths,
            quantiles=(0.16,.5, 0.84),
            bins = 30, \
            color='k', \
            show_titles=True, \
            title_fmt = '.2f', \
            smooth1d=False, \
            smooth=True, \
            fill_contours=True, \
            plot_contours =True, \
            plot_density=True, \
            use_mathtext=True,
            title_kwargs={"fontsize": 11},
            label_kwargs={'fontsize':16})
        for ax in fig.get_axes():
            ax.tick_params(axis='both',labelsize=14)
        return(fig)


    def plot_object(self, bands='all', savefig=False,plot3D=False,
                    filename='mySN', orientation='horizontal',method='separate',
                    showModel=False,showFit=False,showMicro=False,**kwargs):
        """Plot the multiply-imaged SN light curves and show/save to a file.
            Each subplot shows a single-band light curve, for all images of the SN.

        Parameters
        ----------
        bands : str or :class:`~list` of :class:`~str`
            'all' = plot all bands; or provide a list of bands to plot
        savefig : bool
            boolean to save or not save plot
        plot3D : bool
            boolean to plot in 3D with plotly
        filename : str
            if savefig is True, this is the output filename
        orientation : str
            'horizontal' = all subplots are in a single row
            'vertical' = all subplots are in a single column
        method : str
            Plots the result of separate, series, or color curve method
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
        colors3d=['red','green','blue','black','purple']
        i=0
        leg=[]

        if method=='series':
            if isinstance(bands,str) and bands == 'all':
                bands = set(self.series.table['band'])

            nbands = len(bands)
            if orientation.startswith('v'):
                ncols = 1
                nrows = nbands
            else:
                ncols = nbands
                nrows = 1
            fig=None
            axlist=None
            if plot3D:
                try:
                    import plotly.graph_objects as go
                    from plotly.subplots import make_subplots
                    fig=go.FigureWidget(make_subplots(rows=nrows,cols=ncols,subplot_titles=list(bands),
                                      specs=[[{'type':'scatter3d'}]*ncols]*nrows))

                    n3dPlots=ncols*nrows


                    axlist=[None]*len(bands)
                except RuntimeError:
                    print('Asked for 3D plot but do not have plotly installed, switching to 2D...')
                    fig=None
                    plot3D=False

            if not plot3D or fig is None:
                fig,axlist=plt.subplots(nrows=nrows, ncols=ncols,
                                    sharex=True, sharey=False,figsize=(10,10))
            if nbands==1:
                axlist = [axlist]
            for lc in np.sort([x for x in self.images.keys()]):
                temp=self.series.table[self.series.table['image']==lc]
                if nrows==1:
                    ccol=0
                    crow=1
                else:
                    crow=0
                    ccol=1

                try:
                    delay=self.series.time_delays[lc]
                    delayerr=self.series.time_delay_errors[lc]
                    mag=self.series.magnifications[lc]
                    magerr=self.series.magnification_errors[lc]
                except:
                    delay=0
                    delayerr=[0,0]
                    mag=1
                    magerr=[0,0]
                for b, ax in zip(bands, axlist):
                    if nrows==1:
                        ccol+=1
                    else:
                        crow+=1

                    if b==list(bands)[0]:

                        if plot3D:

                            fig.add_trace(go.Scatter3d(x=temp['time'][temp['band']==b]+delay,y=temp['time'][temp['band']==b],
                                                   z=temp['flux'][temp['band']==b],
                            error_y=dict(symmetric=False,width=4,array=[delayerr[0]]*len(temp['time'][temp['band']==b]),
                                                       arrayminus=[delayerr[1]]*len(temp['time'][temp['band']==b])),
                            error_z=dict(symmetric=False,width=4,array=temp['flux'][temp['band']==b]*\
                                    np.sqrt((temp['fluxerr'][temp['band']==b]/temp['flux'][temp['band']==b])**2+\
                                            (magerr[0]/mag)**2),
                                 arrayminus=temp['flux'][temp['band']==b]* \
                                       np.sqrt((temp['fluxerr'][temp['band']==b]/temp['flux'][temp['band']==b])**2+ \
                                               (magerr[1]/mag)**2)),mode='markers',
                            marker=dict(color=colors3d[i]),name='Image %s'%lc[-1],**kwargs),row=crow,col=ccol)
                        else:
                            leg.append(
                                ax.errorbar(temp['time'][
                                            temp['band']==b],
                                        temp['flux'][
                                            temp['band']==b],
                                        yerr=temp['fluxerr'][
                                            temp['band']==b],
                                        markersize=4, fmt=colors[i]+'.'))
                    else:
                        if plot3D:

                            fig.add_trace(go.Scatter3d(x=temp['time'][temp['band']==b]+delay,y=temp['time'][temp['band']==b],
                                                       z=temp['flux'][temp['band']==b],
                                                       error_y=dict(symmetric=False,width=4,array=[delayerr[0]]*len(temp['time'][temp['band']==b]),
                                                                    arrayminus=[delayerr[1]]*len(temp['time'][temp['band']==b])),
                                                       error_z=dict(symmetric=False,width=4,array=temp['flux'][temp['band']==b]* \
                                                                                                  np.sqrt((temp['fluxerr'][temp['band']==b]/temp['flux'][temp['band']==b])**2+ \
                                                                                                          (magerr[0]/mag)**2),
                                                                    arrayminus=temp['flux'][temp['band']==b]* \
                                                                               np.sqrt((temp['fluxerr'][temp['band']==b]/temp['flux'][temp['band']==b])**2+ \
                                                                                       (magerr[1]/mag)**2)),mode='markers',
                                                       marker=dict(color=colors3d[i]),showlegend=False,**kwargs),row=crow,col=ccol)
                        else:

                            ax.errorbar(temp['time'][
                                        temp['band']==b],
                                    temp['flux'][
                                        temp['band']==b],
                                    yerr=temp['fluxerr'][
                                        temp['band']==b],
                                    markersize=4, fmt=colors[i]+'.')
                    if showFit:
                        if plot3D:
                            fig.add_trace(go.Scatter3d(x=np.arange(np.min(temp['time'][temp['band']==b]),np.max(temp['time'][temp['band']==b]),1)+delay,
                                            y=np.arange(np.min(temp['time'][temp['band']==b]),np.max(temp['time'][temp['band']==b]),1),
                                           z=self.series.fits.model.bandflux(b,np.arange(np.min(temp['time'][temp['band']==b]),np.max(temp['time'][temp['band']==b]),1),
                                                                             zp=temp['zp'][temp['band']==b][0],
                                                                             zpsys=temp['zpsys'][temp['band']==b][0]),
                                           mode='lines',line=dict(color='yellow',width=8),
                                           showlegend=False,**kwargs),row=crow,col=ccol)
                        else:
                            ax.plot(np.arange(np.min(temp['time'][temp['band']==b]),np.max(temp['time'][temp['band']==b]),1),
                                self.series.fits.model.bandflux(b,np.arange(np.min(temp['time'][temp['band']==b]),np.max(temp['time'][temp['band']==b]),1),
                                                                  zp=temp['zp'][temp['band']==b][0],
                                                                  zpsys=temp['zpsys'][temp['band']==b][0]),color='y')
                    if not plot3D:
                        ax.text(0.95, 0.95, b.upper(), fontsize='large',
                            transform=ax.transAxes, ha='right', va='top')


                i+=1
        elif method =='color':
            n3dPlots=1
            if isinstance(bands,str) and bands=='all':
                if len([x for x in self.color.table.colnames if '-' in x and '_' not in x])!=1:
                    print("Want to plot color curves but need 2 bands specified.")
                    sys.exit(1)
                else:
                    colname=[x for x in self.color.table.colnames if '-' in x and '_' not in x][0]
                    bands=[colname[:colname.find('-')],colname[colname.find('-')+1:]]
            elif len(bands) !=2:
                print("Want to plot color curves but need 2 bands specified.")
                sys.exit(1)

            if plot3D:
                try:
                    import plotly.graph_objects as go
                    from plotly.subplots import make_subplots
                    fig=go.FigureWidget()
                    ax=None

                except RuntimeError:
                    print('Asked for 3D plot but do not have plotly installed, switching to 2D...')
                    fig=None
                    plot3D=False

            if not plot3D or fig is None:
                fig=plt.figure(figsize=(10,10))
                ax=fig.gca()
            for lc in np.sort([x for x in self.images.keys()]):
                temp=self.color.table[self.color.table['image']==lc]
                try:
                    delay=self.color.time_delays[lc]
                    delayerr=self.color.time_delay_errors[lc]

                except:
                    delay=0
                    delayerr=[0,0]

                if plot3D:

                    fig.add_trace(go.Scatter3d(x=temp['time']+delay,y=temp['time'],
                                               z=temp[bands[0]+'-'+bands[1]],
                                               error_y=dict(symmetric=False,width=4,array=[delayerr[0]]*len(temp['time']),
                                                            arrayminus=[delayerr[1]]*len(temp['time'])),
                                               error_z=dict(symmetric=True,width=4,
                                                            array=temp[bands[0]+'-'+bands[1]+'_err']),mode='markers',
                                               marker=dict(color=colors3d[i]),name='Image %s'%lc[-1],**kwargs))
                else:
                    ax.errorbar(temp['time'],temp[bands[0]+'-'+bands[1]],yerr=temp[bands[0]+'-'+bands[1]+'_err'],
                            markersize=4, fmt=colors[i]+'.')

                    ax.text(0.95, 0.95, bands[0].upper()+'-'+bands[1].upper(), fontsize='large',
                        transform=ax.transAxes, ha='right', va='top')

                i+=1
                if showFit:
                    mod_time=np.arange(np.min(self.color.table['time']),np.max(self.color.table['time']),1)
                    modCol=self.color.fits.model.color(bands[0],bands[1],self.color.table['zpsys'][0],mod_time)
                    if plot3D:
                        fig.add_trace(go.Scatter3d(x=mod_time+delay,
                                                   y=mod_time,
                                                   z=modCol,
                                                   mode='lines',line=dict(color='yellow',width=8),
                                                   showlegend=False,**kwargs))
                    else:
                        ax.plot(mod_time,modCol,color='y')
        else:
            if isinstance(bands,str):
                if bands == 'all':
                    bands = self.bands
                else:
                    bands=[bands]


            nbands = len(bands)
            if orientation.startswith('v'):
                ncols = 1
                nrows = nbands
            else:
                ncols = nbands
                nrows = 1

            n3dPlots=nrows*ncols

            fig=None
            axlist=None
            if plot3D:
                try:
                    import plotly.graph_objects as go
                    from plotly.subplots import make_subplots
                    fig=go.FigureWidget(make_subplots(rows=nrows,cols=ncols,subplot_titles=list(bands),
                                      specs=[[{'is_3d': True}]*ncols]*nrows))



                    axlist=[None]*len(bands)
                except RuntimeError:
                    print('Asked for 3D plot but do not have plotly installed, switching to 2D...')
                    fig=None
                    plot3D=False

            if not plot3D or fig is None:
                fig,axlist=plt.subplots(nrows=nrows, ncols=ncols,
                                    sharex=True, sharey=True,figsize=(10,10))
            if nbands==1:
                axlist = [axlist]
            microAx={b:False for b in bands}
            for lc in np.sort([x for x in self.images.keys()]):
                if nrows==1:
                    ccol=0
                    crow=1
                else:
                    crow=0
                    ccol=1

                try:
                    delay=self.parallel.time_delays[lc]
                    delayerr=self.parallel.time_delay_errors[lc]
                    if showFit:
                        mag=self.parallel.magnifications[lc]
                        magerr=self.parallel.magnification_errors[lc]
                    else:
                        mag=1
                        magerr=[0,0]
                except:
                    delay=0
                    delayerr=[0,0]
                    mag=1
                    magerr=[0,0]


                for b, ax in zip(bands, axlist):
                    if nrows==1:
                        ccol+=1
                    else:
                        crow+=1
                    if b==list(bands)[0]:
                        if plot3D:
                            temp=self.images[lc].table
                            fig.add_trace(go.Scatter3d(x=temp['time'][temp['band']==b],y=temp['time'][temp['band']==b]-delay,
                                   z=temp['flux'][temp['band']==b]/mag,
                                   error_y=dict(symmetric=False,width=4,array=[delayerr[0]]*len(temp['time'][temp['band']==b]),
                                                arrayminus=[delayerr[1]]*len(temp['time'][temp['band']==b])),
                                   error_z=dict(symmetric=False,width=4,array=temp['flux'][temp['band']==b]* \
                                                                              np.sqrt((temp['fluxerr'][temp['band']==b]/temp['flux'][temp['band']==b])**2+ \
                                                                                      (magerr[0]/mag)**2),
                                                arrayminus=temp['flux'][temp['band']==b]* \
                                                           np.sqrt((temp['fluxerr'][temp['band']==b]/temp['flux'][temp['band']==b])**2+ \
                                                                   (magerr[1]/mag)**2)),mode='markers',
                                   marker=dict(color=colors3d[i]),name='Image %s'%lc[-1],**kwargs),row=crow,col=ccol)
                        else:
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
                        if plot3D:
                            fig.add_trace(go.Scatter3d(x=temp['time'][temp['band']==b],y=temp['time'][temp['band']==b]-delay,
                               z=temp['flux'][temp['band']==b]/mag,
                               error_y=dict(symmetric=False,width=4,array=[delayerr[0]]*len(temp['time'][temp['band']==b]),
                                            arrayminus=[delayerr[1]]*len(temp['time'][temp['band']==b])),
                               error_z=dict(symmetric=False,width=4,array=temp['flux'][temp['band']==b]* \
                                                                          np.sqrt((temp['fluxerr'][temp['band']==b]/temp['flux'][temp['band']==b])**2+ \
                                                                                  (magerr[0]/mag)**2),
                                            arrayminus=temp['flux'][temp['band']==b]* \
                                                       np.sqrt((temp['fluxerr'][temp['band']==b]/temp['flux'][temp['band']==b])**2+ \
                                                               (magerr[1]/mag)**2)),mode='markers',
                               marker=dict(color=colors3d[i]),name='Image %s'%lc[-1],
                                          showlegend=False,**kwargs),row=crow,col=ccol)
                        else:
                            ax.errorbar(self.images[lc].table['time'][
                                        self.images[lc].table['band']==b],
                                    self.images[lc].table['flux'][
                                        self.images[lc].table['band']==b],
                                    yerr=self.images[lc].table['fluxerr'][
                                        self.images[lc].table['band']==b],
                                    markersize=4, fmt=colors[i]+'.')
                    if showFit:
                        time_model = np.arange(self.images[lc].table['time'].min(),
                                               self.images[lc].table['time'].max(),
                                               0.1)
                        flux_model=self.images[lc].fits.model.bandflux(b,time_model,
                                                                       self.images[lc].table['zp'][self.images[lc].table['band']==b][0],
                                                                       self.images[lc].table['zpsys'][self.images[lc].table['band']==b][0])
                        if plot3D:
                            fig.add_trace(go.Scatter3d(x=time_model,
                                                       y=time_model-delay,
                                                       z=flux_model/mag,
                                                       mode='lines',line=dict(color='yellow',width=8),
                                                       showlegend=False,**kwargs),row=crow,col=ccol)
                        else:
                            ax.plot(time_model,flux_model,'k-')
                    if not plot3D:
                        ax.text(0.95, 0.95, b.upper(), fontsize='large',
                            transform=ax.transAxes, ha='right', va='top')


                    if showMicro and not plot3D:
                        time_model=np.arange(self.images[lc].table['time'].min(),
                                             self.images[lc].table['time'].max(),
                                             0.1)
                        if microAx[b] is False:
                            ax_divider = make_axes_locatable(ax)
                            ax_ml = ax_divider.append_axes("bottom", size="25%", pad=.4)
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
                        if plot3D:
                            fig.add_trace(go.Scatter3d(x=time_shifted,
                                                       y=time_model,
                                                       z=flux_magnified,
                                                       mode='lines',line=dict(color='orange',width=8),
                                                       showlegend=False,**kwargs))
                        else:
                            ax.plot(time_model, flux_magnified, 'r-')


                i+=1

        if plot3D:

            zname = 'Flux' if method !='color' else bands[0]+'-'+bands[1]+' Color'
            tempscene=dict(aspectmode='cube',camera=dict(projection=dict(type='orthographic'),eye=dict(
                x=0,
                y=10,
                z=0,
            )),
                           xaxis=dict(
                               gridcolor='rgb(255, 255, 255)',
                               zerolinecolor='rgb(255, 255, 255)',
                               showbackground=True,
                               backgroundcolor='rgb(230, 230, 230)',
                               title='Time (Observer Frame)',
                               mirror=False,
                               autorange='reversed'


                           ),
                           yaxis=dict(
                               gridcolor='rgb(255, 255, 255)',
                               zerolinecolor='rgb(255, 255, 255)',
                               showbackground=True,
                               backgroundcolor='rgb(230, 230, 230)',
                               title='Corrected Time (Observer Frame)',
                               mirror=False


                           ),
                           zaxis=dict(
                               gridcolor='rgb(255, 255, 255)',
                               zerolinecolor='rgb(255, 255, 255)',
                               showbackground=True,
                               backgroundcolor='rgb(230, 230, 230)',
                               title=zname,
                               mirror=False


                           ))
            
            scenes=dict([])
            for i in range(n3dPlots):
                if i>0:
                    key='scene'+str(i+1)
                    scenes[key]=tempscene
                else:
                    key='scene'
                    scenes[key]=tempscene




            fig.update_layout(
                **scenes
                )


        else:
            plt.figlegend(leg,np.sort(['$'+x+'$' for x in self.images.keys()]), frameon=False,
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


def table_factory(tables,telescopename="Unknown",object_name=None):
    """This function will create a new curve object using an astropy table or tables.

    Parameters
    ----------
    tables : `~astropy.table.Table` or :class:`~list`
        Astropy table with all of your data from your data file, or a list of such tables.
    telescopename : str
        Name of telescope for labeling purposes inside curve object
    object_name : str
        Name of object for labeling purposes inside curve object
        (e.g. SN2006jf, etc.)

    Returns
    -------
    curve : :class:`~sntd.curve`

    """

    new_SN=curveDict(telescopename,object_name)
    if not isinstance(tables,(list,tuple)):
        tables=[tables]
    for table in tables:
        newlc=curve()
        table=standardize_table_colnames(table)
        newlc.bands={x for x in table[get_default_prop_name('band')]}
        zp_dict=dict([])
        zpsys_dict=dict([])
        for band in newlc.bands:
            zp=set(table[get_default_prop_name('zp')][table[get_default_prop_name('band')]==band])
            zpsys=set(table[get_default_prop_name('zpsys')][table[get_default_prop_name('band')]==band])
            zp_dict[band]=zp.pop() if len(zp)==1 else np.array(zp)
            zpsys_dict[band]=zpsys.pop() if len(zpsys)==1 else np.array(zpsys)


        newlc.table=table

        newlc.telescopename = telescopename
        newlc.object = object_name
        if 'image' in newlc.table.colnames:
            im_key=newlc.table['image'][0]
            new_SN.add_curve(newlc,key=im_key)
        else:
            new_SN.add_curve(newlc)

    return new_SN


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

    Returns
    -------
    None
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
    curve : :class:`~sntd.curve` or :class:`~sntd.curveDict`
    """
    return(_switch(os.path.splitext(filename)[1])(filename,**kwargs))


def _read_pickle(filename,telescopename="Unknown",object="Unknown",**kwargs):
    try:
        return (pickle.load(open(filename,'rb')))
    except:
        return (cPickle.load(open(filename,'rb')))



def standardize_table_colnames(table):
    for col in table.colnames:
        if col != get_default_prop_name(col.lower()):
            table.rename_column(col, get_default_prop_name(col.lower()))
    return table

def _read_data(filename,**kwargs):
    table_done=True
    try:
        table = sncosmo.read_lc(filename, masked=True,**kwargs)
        for col in table.colnames:
            if col != get_default_prop_name(col.lower()):
                table.rename_column(col, get_default_prop_name(col.lower()))
    except:
        try:
            table = ascii.read(filename, masked=True,**kwargs)
            for col in table.colnames:
                if col != get_default_prop_name(col.lower()):
                    table.rename_column(col, get_default_prop_name(col.lower()))
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
                colnames = odict.fromkeys([get_default_prop_name(x.lower()) for x in line])
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

    for col in [get_default_prop_name(x) for x in ['band','zp','zpsys']]:
        if col not in colnames:
            temp=kwargs.get(col,None)
            if temp is not None:
                table[col]=temp
            else:
                print('Column "%s" is not in your file and you did not define it in kwargs.'%col)
                sys.exit(1)
    table=standardize_table_colnames(table)
    bnds = {x for x in table[get_default_prop_name('band')]}
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
            table.mask[table[get_default_prop_name('band')]==band]=True
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
    table[get_default_prop_name('mag')] = np.asarray(map(lambda x, y: -2.5 * np.log10(x) + y, table[get_default_prop_name('flux')],
                                  table[get_default_prop_name('zp')]))
    table[get_default_prop_name('magerr')] = np.asarray(map(lambda x, y: 2.5 * np.log10(np.e) * y / x, table[get_default_prop_name('flux')],
                                     table[get_default_prop_name('fluxerr')]))
    table[get_default_prop_name('magerr')][np.isnan(table[get_default_prop_name('mag')])]=np.nan
    return table


def _mag_to_flux(table):
    table[get_default_prop_name('flux')] = np.asarray(
        map(lambda x, y: 10 ** (-.4 * (x -y)), table[get_default_prop_name('mag')],
            table[get_default_prop_name('zp')]))
    table[get_default_prop_name('fluxerr')] = np.asarray(
        map(lambda x, y: x * y / (2.5 * np.log10(np.e)), table[get_default_prop_name('magerr')],
            table[get_default_prop_name('flux')]))
    return table



