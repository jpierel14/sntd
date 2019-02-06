import os,sys,math,subprocess

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

from .util import __dir__,__current_dir__

__all__=['realizeMicro']
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
    thefile=open(os.path.join(__dir__,'microlens','input'),'wb')

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
        plt.savefig('sntd_microlensing.pdf',format='pdf',overwrite=True)
        #plt.show()
        plt.clf()
        plt.close()

    return(mu)

def getDiffCurve(time,num,default=True):

    tab=ascii.read(os.path.join(__dir__,'data','diff'+str(num)+'.dat'))
def getDiffCurve(time, bandset='bessell'):
    """Read in a microlensing difference curve from a data file.
    """
    #TODO : interpolate to account for redshift of the lens and source
    num=np.random.randint(1,5)
    if bandset.lower()=='bessell':
        tab=ascii.read(os.path.join(__dir__,'data','diff'+str(num)+'.dat'))
    elif bandset.lower()=='hst':
        tab = ascii.read(
            os.path.join(__dir__, 'data/hstmicrolensing',
                         'diff' + str(num) + '.dat'))
    else:
        raise RuntimeError("bandset must be 'bessell' or 'hst'")
    outTab=Table()
    outTab['time']=time
    for band in [x for x in tab.colnames if x != 'time']:
        spl=splrep(tab['time']+time[0],10**(-.4*tab[band]))
        outTab[band]=splev(time,spl)


    return(outTab)




