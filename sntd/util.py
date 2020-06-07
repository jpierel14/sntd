#!/Users/jpierel/anaconda3/envs/astro2/bin python2

import os,sncosmo,glob,sys,subprocess,time
from astropy.io import ascii
import numpy as np
from collections import OrderedDict as odict
import scipy
import matplotlib.pyplot as plt


from scipy.interpolate import splrep,splev
from copy import copy
import tarfile,os,pickle
from scipy.stats import rv_continuous


__current_dir__=os.path.abspath(os.getcwd())
__filedir__=os.path.abspath(os.path.dirname(__file__))

NORMAL = 0    # use python zip libraries
PROCESS = 1   # use (zcat, gzip) or (bzcat, bzip2)
PARALLEL = 2  # (pigz -dc, pigz) or (pbzip2 -dc, pbzip2)


_props=odict([
    ('time',{'mjd', 'mjdobs', 'jd', 'time', 'date', 'mjd_obs','mhjd','jds'}),
    ('band',{'filter', 'band', 'flt', 'bandpass'}),
    ('flux',{'flux', 'f','fluxes'}),
    ('fluxerr',{'flux_error', 'fluxerr', 'fluxerror', 'fe', 'flux_err','fluxerrs'}),
    ('zp',{'zero_point','zp', 'zpt', 'zeropoint'}),
    ('zpsys',{'zpsys', 'magsys', 'zpmagsys'}),
    ('mag',{'mag','magnitude','mags'}),
    ('magerr',{'magerr','magerror','magnitudeerror','magnitudeerr','magerrs'})
])

_sncosmo_snana= [('snana-2004fe', 'SN Ic', 'CSP-2004fe.SED'),
          ('snana-2004gq', 'SN Ic', 'CSP-2004gq.SED'),
          ('snana-sdss004012', 'SN Ic', 'SDSS-004012.SED'),  # no IAU name
          ('snana-2006fo', 'SN Ic', 'SDSS-013195.SED'),  # PSNID
          ('snana-sdss014475', 'SN Ic', 'SDSS-014475.SED'),  # no IAU name
          ('snana-2006lc', 'SN Ic', 'SDSS-015475.SED'),
          ('snana-2007ms', 'SN II-pec', 'SDSS-017548.SED'),  # type Ic in SNANA
          ('snana-04d1la', 'SN Ic', 'SNLS-04D1la.SED'),
          ('snana-04d4jv', 'SN Ic', 'SNLS-04D4jv.SED'),
          ('snana-2004gv', 'SN Ib', 'CSP-2004gv.SED'),
          ('snana-2006ep', 'SN Ib', 'CSP-2006ep.SED'),
          ('snana-2007Y', 'SN Ib', 'CSP-2007Y.SED'),
          ('snana-2004ib', 'SN Ib', 'SDSS-000020.SED'),
          ('snana-2005hm', 'SN Ib', 'SDSS-002744.SED'),  # PSNID
          ('snana-2006jo', 'SN Ib', 'SDSS-014492.SED'),  # PSNID
          ('snana-2007nc', 'SN Ib', 'SDSS-019323.SED'),
          ('snana-2004hx', 'SN IIP', 'SDSS-000018.SED'),  # PSNID
          ('snana-2005gi', 'SN IIP', 'SDSS-003818.SED'),  # PSNID
          ('snana-2006gq', 'SN IIP', 'SDSS-013376.SED'),
          ('snana-2006kn', 'SN IIP', 'SDSS-014450.SED'),
          ('snana-2006jl', 'SN IIP', 'SDSS-014599.SED'),  # PSNID
          ('snana-2006iw', 'SN IIP', 'SDSS-015031.SED'),
          ('snana-2006kv', 'SN IIP', 'SDSS-015320.SED'),
          ('snana-2006ns', 'SN IIP', 'SDSS-015339.SED'),
          ('snana-2007iz', 'SN IIP', 'SDSS-017564.SED'),
          ('snana-2007nr', 'SN IIP', 'SDSS-017862.SED'),
          ('snana-2007kw', 'SN IIP', 'SDSS-018109.SED'),
          ('snana-2007ky', 'SN IIP', 'SDSS-018297.SED'),
          ('snana-2007lj', 'SN IIP', 'SDSS-018408.SED'),
          ('snana-2007lb', 'SN IIP', 'SDSS-018441.SED'),
          ('snana-2007ll', 'SN IIP', 'SDSS-018457.SED'),
          ('snana-2007nw', 'SN IIP', 'SDSS-018590.SED'),
          ('snana-2007ld', 'SN IIP', 'SDSS-018596.SED'),
          ('snana-2007md', 'SN IIP', 'SDSS-018700.SED'),
          ('snana-2007lz', 'SN IIP', 'SDSS-018713.SED'),
          ('snana-2007lx', 'SN IIP', 'SDSS-018734.SED'),
          ('snana-2007og', 'SN IIP', 'SDSS-018793.SED'),
          ('snana-2007ny', 'SN IIP', 'SDSS-018834.SED'),
          ('snana-2007nv', 'SN IIP', 'SDSS-018892.SED'),
          ('snana-2007pg', 'SN IIP', 'SDSS-020038.SED'),
          ('snana-2006ez', 'SN IIn', 'SDSS-012842.SED'),
          ('snana-2006ix', 'SN IIn', 'SDSS-013449.SED'),
          ('s11-2004hx', 'SN IIL/P', 'S11_SDSS-000018.SED'),
          ('s11-2005lc', 'SN IIP', 'S11_SDSS-001472.SED'),
          ('s11-2005hl', 'SN Ib', 'S11_SDSS-002000.SED'),
          ('s11-2005hm', 'SN Ib', 'S11_SDSS-002744.SED'),
          ('s11-2005gi', 'SN IIP', 'S11_SDSS-003818.SED'),
          ('s11-2006fo', 'SN Ic', 'S11_SDSS-013195.SED'),
          ('s11-2006jo', 'SN Ib', 'S11_SDSS-014492.SED'),
          ('s11-2006jl', 'SN IIP', 'S11_SDSS-014599.SED')]

def get_models_by_sntype(snType):
    mod,types=np.loadtxt(os.path.join(__filedir__,'data','sncosmo','models.ref'),dtype='str',unpack=True)
    modDict={mod[i]:types[i] for i in range(len(mod))}
    mods = [x[0] for x in sncosmo.models._SOURCES._loaders.keys() if x[0] in modDict.keys() and modDict[x[0]][:len(snType)]==snType]
    return(mods)
        

def snana_to_sncosmo(snana_mod):
    for sncosmo_name,typ,snana_name in _sncosmo_snana:
        if snana_name==snana_mod:
            return[typ,sncosmo_name]
    return None

def sncosmo_to_snana(sncosmo_mod):
    for sncosmo_name,typ,snana_name in _sncosmo_snana:
        if sncosmo_name==sncosmo_mod:
            return[typ,snana_name]
    return None

def load_example_data():
    example_files=glob.glob(os.path.join(__filedir__,'data','examples','*.dat'))
    return(ascii.read(example_files[0]),ascii.read(example_files[1]))

def load_example_misn():
    example_file=glob.glob(os.path.join(__filedir__,'data','examples','*.pkl'))
    return(pickle.load(open(example_file[0],'rb')))

def load_batch_fit_names(folder_name='.',verbose=True):
    tars=glob.glob(os.path.join(folder_name,'sntd_fits*.tar.gz'))
    all_names={}
    for tar_fname in tars:
        tar=tarfile.open(tar_fname,'r')
        all_fits=tar.getmembers()
        if verbose:
            print('Found %i fits, loading...'%len(all_fits))
        tar.close()
        all_names[tar_fname]=[x.name for x in all_fits]
    return all_names

def load_batch_fit(fit_name,folder=None,tar_dict=None):
    if tar_dict is None:
        folder='.' if folder is None else folder
        tar_dict=load_batch_fit_names(folder)
    to_return=None
    for tar_fname in tar_dict.keys():
        if fit_name not in tar_dict[tar_fname]:
            continue
        
        tar=tarfile.open(tar_fname,'r')
        f=tar.extractfile(fit_name).read()
        dat=pickle.loads(f)
        tar.close()
        return dat
    print('Did not find your file')
    return

def check_table_quality(table,min_n_bands=1,min_n_points_per_band=1,clip=False):
    ngood_bands=0
    for b in np.unique(table['band']):
        temp_n_for_b=len(table[table['band']==b])
        if temp_n_for_b<min_n_points_per_band:
            if clip:
                table=table[table['band']!=b]
        else:
            ngood_bands+=1
    if ngood_bands<min_n_bands:
        return table,False
    return table,True

def run_sbatch(folder_name,script_name_init,script_name,total_jobs,max_batch_jobs,n_per_node,wait_for_batch,parallelize,ncurves):
    fits_output=tarfile.open(os.path.join(os.path.abspath(folder_name),'sntd_fits.tar.gz'),mode='w')
                
    result=subprocess.call(['sbatch',os.path.join(os.path.abspath(folder_name),
                                                           script_name_init)])
    if wait_for_batch:
        printProgressBar(0,total_jobs)
    ndone=0
    nadded=min(total_jobs,max_batch_jobs)
    saved_fits=0
    tarfit_ind=0
    if parallelize is not None:
        n_per_file=1
    else:
        n_per_file=n_per_node
    
    while True:
        time.sleep(10) #update every 10 seconds
        output=glob.glob(os.path.join(os.path.abspath(folder_name),'sntd_fit*.pkl'))
        saved_fits+=len(output)
        if len(output)>0:
            if int(saved_fits*n_per_file)>=50000*(tarfit_ind+1):
                fits_output.close()
                fits_output=tarfile.open(os.path.join(os.path.abspath(folder_name),'sntd_fits_%i.tar.gz'%tarfit_ind),mode='w')
                tarfit_ind+=1
            for filename in output:
                fits_output.add(filename)
                os.remove(filename)
            if nadded<total_jobs:
                for i in range(math.ceil(len(output)/(n_per_node/n_per_file))):
                    if nadded>total_jobs-1:
                        continue
                    result=subprocess.call(['sbatch',os.path.join(os.path.abspath(folder_name),
                                                             script_name),str(nadded)],stdout=subprocess.DEVNULL)
                    nadded+=1

            if wait_for_batch:
                printProgressBar(saved_fits/(n_per_node/n_per_file),total_jobs)
        if saved_fits>=ncurves:
            break
    fits_output.close()
    if verbose:
        print('Done!')
    return 

def make_sbatch(partition=None,njobs=None,njobstotal=None,python_path=None,init=False,folder=None,parallelize=None,microlensing_cores=None):
    if njobs is None:
        print("Batch mode requires a number of jobs!")
        sys.exit(1)
    if init:
        if njobstotal is None:
            print("Batch mode requires a total number of jobs!")
        n=0
        add=''
        done=False
        while not done:
            try:
                folder_name='batch_output%s'%add
                os.mkdir(folder_name)
                done=True
            except:
                add=str(n)
                n+=1
            if n>50:
                print('Having trouble making batch output folder.')
                sys.exit(1)
    else:
        folder_name=folder


    if python_path is None:
        python_path=subprocess.check_output("which python", shell=True).decode('utf-8').strip('\n')
    if not init:
        with open(os.path.join(__filedir__,'batch','sbatch_job.BATCH')) as f:
            sbatch=f.read()
        if parallelize is None:
            pyfile='run_sntd.py'
        else:
            pyfile='run_sntd_par.py'
    else:
        with open(os.path.join(__filedir__,'batch','sbatch_job_init.BATCH')) as f:
            sbatch=f.read()
        if parallelize is None:
            pyfile='run_sntd_init.py'
        else:
            pyfile='run_sntd_init_par.py'


    sbatch=sbatch.replace('pyjob%j.out',os.path.join(folder_name,'pyjob%j.out'))
    if partition is not None:
        sbatch=sbatch.replace('partition','#SBATCH -p %s'%partition)
    else:
        sbatch=sbatch.replace('partition','')
    sbatch=sbatch.replace('myPython',python_path)
    sbatch=sbatch.replace('run_sntd.py',os.path.join(os.path.abspath(folder_name),pyfile))
    if init:
        sbatch=sbatch.replace('njobstotal','0-%i'%(njobstotal-1))
        sbatch=sbatch.replace('njobs','%i'%njobs)
    if parallelize is not None:
        sbatch=sbatch.replace('ncores',str(parallelize))
    elif microlensing_cores is not None:
        sbatch=sbatch.replace('ncores',str(microlensing_cores))
    else:
        sbatch=sbatch.replace('ncores',str(1))

    if not init:
        with open(os.path.join(folder_name,'sbatch_job.BATCH'),'w') as f:
            f.write(sbatch)
        return('sbatch_job.BATCH',folder_name)
    else:
        with open(os.path.join(folder_name,'sbatch_job_init.BATCH'),'w') as f:
            f.write(sbatch)
        return('sbatch_job_init.BATCH',folder_name)

def plot(plot_type,x,y=None,yerr=None,xerr=None,ax=None,x_lab='',y_lab='',fontsize=18,figsize=(12,12),
         x_name=None,y_name=None,label_name=None,**kwargs):
    if ax is None and plot_type != 'joint':
        fig=plt.figure(figsize=figsize)
        ax=fig.gca()

    if plot_type=='scatter':
        ax.scatter(x,y,**kwargs)
    elif plot_type=='plot':
        ax.plot(x,y,**kwargs)
    elif plot_type=='errorbar':
        ax.errorbar(x,y,xerr=xerr,yerr=yerr,**kwargs)
    elif plot_type=='hist':
        ax.hist(x,**kwargs)
    elif plot_type=='joint':
        g=multivariateGrid(x_name, y_name, label_name, df=x)
        fig=g.ax_joint.__dict__['figure']
        ax=fig.gca()
        fig.set_size_inches(figsize[0],figsize[1])
    else:
        raise RuntimeError('What plot are you trying to do.')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    ax.set_xlabel(x_lab,fontsize=fontsize)
    ax.set_ylabel(y_lab,fontsize=fontsize)
    return(ax)

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

class posterior(rv_continuous):
    "Skewed Normal Distribution"
    def _pdf(self,x,samples,weights):
        pdf,edges=np.histogram(samples,weights=weights,
                               bins=30,density=True)
        
        func=scipy.interpolate.interp1d([(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)],
                                        pdf/np.max(pdf),fill_value=0,bounds_error=False)
        return(func(x))

    def _argcheck(self,*args):
        return True

def guess_magnifications(curves,referenceImage):
    """Guess t0 and amplitude of the model based on the data.

    Assumes the data has been standardized."""
    ref=None
    mags=dict([])
    all_ims=[referenceImage]
    for k in [x for x in curves.images.keys() if x!=referenceImage]:
        all_ims.append(k)
    for k in all_ims:
        if not curves.images[k].fits:
            bestRatio=-np.inf
            bestBand=None
            for b in np.unique(curves.images[k].table['band']):
                ratio=np.abs(np.max(curves.images[k].table['flux'][curves.images[k].table['band']==b])/np.min(curves.images[k].table['flux'][curves.images[k].table['band']==b]))
                if ratio>bestRatio:
                    bestRatio=ratio
                    bestBand=b
            maxTime,maxValue=_findMax(curves.images[k].table['time'][curves.images[k].table['band']==b],curves.images[k].table['flux'][curves.images[k].table['band']==b])


            if k==referenceImage:
                ref=copy(maxValue)
            mags[k]=maxValue/ref

        else:
            if 'x0' in curves.images[k].fits.res.vparam_names:
                amplitude='x0'
            else:
                amplitude='amplitude'
            mag=curves.images[k].fits.model.get(amplitude)
            if k==referenceImage:
                ref=copy(mag)
            mags[k]=mag/ref
    return mags


def guess_time_delays(curves,referenceImage):
    #tds=colorFit(curves)
    #if tds:
    #    return tds
    ref=None
    tds=dict([])
    all_ims=[referenceImage]
    for k in [x for x in curves.images.keys() if x!=referenceImage]:
        all_ims.append(k)
    for k in all_ims:
        if not curves.images[k].fits:
            maxValue,maxFlux=_findMax(curves.images[k].table['time'],curves.images[k].table['flux'])
            if k==referenceImage:
                ref=copy(maxValue)
            tds[k]=maxValue-ref
        else:
            t0=curves.images[k].fits.model.get('t0')
            if k==referenceImage:
                ref=copy(t0)
            tds[k]=t0-ref
    return tds

def _findMax(time,curve):
    #TODO check edge cases
    t0=np.where(curve==np.max(curve))[0][0]
    if t0==0:
        #return(time[0])
        return (time[0],curve[0])

    elif t0==len(time)-1:

        return(time[-1],curve[-1])
    #   return (None,None)

    else:
        fit=splrep(time[t0-1:t0+2],curve[t0-1:t0+2],k=2)

    interptime=np.linspace(time[t0-1],time[t0+1],100)
    flux=splev(interptime,fit)
    return(interptime[flux==np.max(flux)],np.max(flux))


def _findMin(time,curve):
    #TODO check edge cases
    t0=np.where(curve==np.min(curve))[0][0]

    if t0==0:
        #return(time[0])
        return (None,None)

    elif t0==len(time)-1:
        #return(time[-1])
        return (None,None)

    else:
        fit=splrep(time[t0-1:t0+2],curve[t0-1:t0+2],k=2)

    interptime=np.linspace(time[t0-1],time[t0+1],100)
    flux=splev(interptime,fit)
    return(interptime[flux==np.min(flux)],np.max(flux))

def flux_to_mag(table,bandDict,zpsys='AB'):
    """
    Accepts an astropy table of flux data and does the conversion to mags (flux in ergs/s/cm^2/AA)

    :param table: Table containing flux, flux error error, and band columns
    :type table: astropy.Table
    :param bandDict: translates band to sncosmo.Bandpass object (i.e. 'U'-->bessellux)
    :type bandDict: dict
    :param zpsys: magnitude system
    :type zpsys: str,optional
    :returns: astropy.Table object with mag and mag error added
    """
    table=table[table['flux']>0]
    ms=sncosmo.get_magsystem(zpsys)
    table[_get_default_prop_name('mag')] = np.asarray(map(lambda x, y: ms.band_flux_to_mag(x,y), table[_get_default_prop_name('flux')],#/sncosmo.constants.HC_ERG_AA
                                                       np.array([bandDict[z] for z in table[_get_default_prop_name('band')]])))
    table[_get_default_prop_name('magerr')] = np.asarray(map(lambda x, y: 2.5 * np.log10(np.e) * y / x, table[_get_default_prop_name('flux')],
                                                          table[_get_default_prop_name('fluxerr')]))
    table=table[np.abs(table['mag']/table['magerr'])>5]
    return(table)

def _cast_str(s):
    try:
        return int(s)
    except:
        try:
            return float(s)
        except:
            return s.strip()


def get_default_prop_name(prop):
    for key,value in _props.items():
        if {prop.lower()} & value:
            return key
    return prop


def _isfloat(value):
    try:
        float(value)
        return True
    except:
        return False

def anyOpen(filename, mode='r', buff=1024*1024, external=PARALLEL):
    if 'r' in mode and 'w' in mode:
        return None
    if filename.startswith('!'):
        import subprocess
        if 'r' in mode:
            return subprocess.Popen(filename[1:], shell=True, bufsize=buff,
                                    stdout=subprocess.PIPE).stdout
        elif 'w' in mode:
            return subprocess.Popen(filename[1:], shell=True, bufsize=buff,
                                    stdin=subprocess.PIPE).stdin
    elif filename.endswith('.bz2'):
        if external == NORMAL:
            import bz2
            return bz2.BZ2File(filename, mode, buff)
        elif external == PROCESS:
            if not which('bzip2'):
                return anyOpen(filename, mode, buff, NORMAL)
            if 'r' in mode:
                return anyOpen('!bzip2 -dc ' + filename, mode, buff)
            elif 'w' in mode:
                return anyOpen('!bzip2 >' + filename, mode, buff)
        elif external == PARALLEL:
            if not which('pbzip2'):
                return anyOpen(filename, mode, buff, PROCESS)
            if 'r' in mode:
                return anyOpen('!pbzip2 -dc ' + filename, mode, buff)
            elif 'w' in mode:
                return anyOpen('!pbzip2 >' + filename, mode, buff)
    elif filename.endswith('.gz'):
        if external == NORMAL:
            import gzip
            return gzip.GzipFile(filename, mode, buff)
        elif external == PROCESS:
            if not which('gzip'):
                return anyOpen(filename, mode, buff, NORMAL)
            if 'r' in mode:
                return anyOpen('!gzip -dc ' + filename, mode, buff)
            elif 'w' in mode:
                return anyOpen('!gzip >' + filename, mode, buff)
        elif external == PARALLEL:
            if not which('pigz'):
                return anyOpen(filename, mode, buff, PROCESS)
            if 'r' in mode:
                return anyOpen('!pigz -dc ' + filename, mode, buff)
            elif 'w' in mode:
                return anyOpen('!pigz >' + filename, mode, buff)
    elif filename.endswith('.xz'):
        if which('xz'):
            if 'r' in mode:
                return anyOpen('!xz -dc ' + filename, mode, buff)
            elif 'w' in mode:
                return anyOpen('!xz >' + filename, mode, buff)
    else:
        return open(filename, mode, buff)
    return None

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
