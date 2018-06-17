import sys,math
import matplotlib.pyplot as plt
import numpy as np

def plotObject(lcs,bands='all',showfig=False,savefig=True,filename='mySN'):

    colors=['r','g','b','k']
    #markers=['.','^','*','8','s','+','D']
    i=0
    nrows=int(math.ceil(len(lcs.bands)/2.))
    fig,ax=plt.subplots(nrows=nrows,ncols=2,sharex=True,sharey=False)
    leg=[]
    for lc in np.sort(lcs.images.keys()):
        #print(lcs.images[lc].simMeta)
        row=0
        col=0
        if bands=='all':
            bands=lcs.images[lc].bands


        for b in bands:
            if b==bands[0]:
                leg.append(ax[row][col].errorbar(lcs.images[lc].table['time'][lcs.images[lc].table['band']==b],
                                                  lcs.images[lc].table['flux'][lcs.images[lc].table['band']==b],
                                                  yerr=lcs.images[lc].table['fluxerr'][lcs.images[lc].table['band']==b],markersize=4,fmt=colors[i]+'.'))
            else:
                ax[row][col].errorbar(lcs.images[lc].table['time'][lcs.images[lc].table['band']==b],
                                      lcs.images[lc].table['flux'][lcs.images[lc].table['band']==b],
                                      yerr=lcs.images[lc].table['fluxerr'][lcs.images[lc].table['band']==b],markersize=4,fmt=colors[i]+'.')
            if lcs.images[lc].ml:
                ax[row][col].plot(lcs.images[lc].table['time'][lcs.images[lc].table['band']==b],lcs.images[lc].table['flux'][lcs.images[lc].table['band']==b]*lcs.images[lc].ml[b],color=colors[i])
            ax[row][col].annotate(b[-1].upper()+' Filter',size=10,xy=(.55,.87), xycoords='axes fraction')

            if row==0:
                if col==0:
                    col=1
                else:
                    row=1
                    col=0
            else:
                col=1
                row+=1

        i+=1

    if not len(lcs.bands)%2==0:
        fig.delaxes(ax[nrows-1][1])
        ax[nrows-2][1].tick_params(axis='x',labelbottom='on',bottom='on')
        plt.figlegend(leg,np.sort(lcs.images.keys()),loc='lower right',fontsize=16)

    fig.text(0.5, 0.02, r'Time (MJD)', ha='center',fontsize=16)
    fig.text(0.04, .5, 'Flux', va='center', rotation='vertical',fontsize=16)
    plt.suptitle('Multiply-Imaged SN "'+lcs.object+'" on the '+lcs.telescopename,fontsize=18)
    if savefig:
        plt.savefig(filename+'.pdf',format='pdf',overwrite=True)
    if showfig:
        plt.show()
    return
