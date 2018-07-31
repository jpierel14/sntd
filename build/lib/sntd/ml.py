import os,sys,tempfile

import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table
from scipy.interpolate import splrep,splev
from .util import __dir__,__current_dir__


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
    os.chdir(os.path.join(__dir__,'microlens'))
    if not verbose:
        #sys.stdout = tempfile.TemporaryFile()
        os.system(r'./microlens > /dev/null')
        os.system(r'./lightcurve >/dev/null')
        #sys.stdout.close()
        #sys.stdout = sys.__stdout__


    else:
        os.system(r'./microlens')
        os.system(r'./lightcurve')
    num=np.loadtxt(os.path.join(__dir__,'microlens','jobnum'),dtype='str')
    #lensPlane=np.array(fits.open(os.path.join(__dir__,'microlens','IRIS'+str(num)+'.fits'))[0].data,dtype=np.float64)
    lensPlane=fits.open(os.path.join(__dir__,'microlens','IRIS'+str(num)+'.fits'))[0].data
    curve=ascii.read(os.path.join(__dir__,'microlens','out_line'),names=('t','xval','yval','pixvalue','maglin','xpix','ypix'))
    os.chdir(__current_dir__)
    return(lensPlane,curve)


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




