import sys,math
import matplotlib.pyplot as plt
import numpy as np

#Color-blind friendly palette via
# https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
_COLORLIST15 = [
    "#000000",# 0 black
    "#004949",# 1 darkteal
    "#009292",# 2 teal
    "#ff6db6",# 3 darkpink
    "#ffb6db",# 4 lightpink
    "#490092",# 5 purple
    "#006ddb",# 6 royalblue
    "#b66dff",# 7 orchid
    "#6db6ff",# 8 tarheelblue
    "#b6dbff",# 9 skyblue
    "#920000",#10 brickred
    "#924900",#11 brown
    "#db6d00",#12 darkorange
    "#24ff24",#13 green
    "#ffff6d",#14 yellow
]
_COLORLIST5 = np.array(_COLORLIST15).take([5,2,12,10,0]).tolist()

def plotTimeDelays(lcs,fig=None,ax=None,band=None,color='b',filename='time_delays',offset=1,savefig=False,showfig=True):
    nrows=len(lcs.images.keys())-1
    ncols=nrows
    if fig is None:
        fig,ax=plt.subplots(nrows=nrows,ncols=ncols,sharex=False,sharey=True,figsize=(14,8))
    elif ax is None:
        ax=fig.gca()
    else:
        ncols=ax.shape[1]
        nrows=ax.shape[0]

    images=np.sort(lcs.images.keys())
    ref=None
    for im in images:
        if lcs.time_delays[im]==0:
            ref=im
            break
    if ref is None:
        raise RuntimeError("Don't have a reference image?")
    ind=np.where(images==ref)
    for i in range(nrows):
        for j in range(ncols):
            if j>i:
                try:
                    fig.delaxes(ax[i][j])
                except:
                    pass
            else:
                delay='S'+str(j+1)+'S'+str(i+2) if i+2<5 else 'S'+str(j+1)+'SX'
                #ax[i][j].errorbar(lcs.time_delays['S'+str(i+2)]-lcs.time_delays['S'+str(j+1)],offset,xerr=lcs.time_delay_errors['S'+str(i+2)]-lcs.time_delay_errors['S'+str(j+1)],label="%s"%lcs.bands,fmt='.',color=color)
                ax[i][j].errorbar(lcs.measurements['t0']['S'+str(j+1)][band]-lcs.measurements['t0'][delay[-2:]][band],offset,
                                  xerr=math.sqrt(lcs.measurements['t0_err']['S'+str(j+1)][band]**2+lcs.measurements['t0_err'][delay[-2:]][band]**2),
                                  fmt='.',color=color,linewidth=2,markersize=12)
                ax[i][j].annotate(delay,size=10,xy=(0.05,.85),xycoords='axes fraction')
                ax[i][j].annotate(str(np.round(lcs.measurements['t0']['S'+str(j+1)][band]-lcs.measurements['t0'][delay[-2:]][band],2))+'$\pm$'+str(np.round(math.sqrt(lcs.measurements['t0_err']['S'+str(j+1)][band]**2+lcs.measurements['t0_err'][delay[-2:]][band]**2),2)),
                                  size=8,xy=(.9999*(lcs.measurements['t0']['S'+str(j+1)][band]-lcs.measurements['t0'][delay[-2:]][band]),offset+.3),xycoords='data',color=color)
                #ax[i][j].set_xlim((.8*(lcs.time_delays['S'+str(i+2)]-lcs.time_delays['S'+str(j+1)]),1.2*(lcs.time_delays['S'+str(i+2)]-lcs.time_delays['S'+str(j+1)])))
                ax[i][j].set_ylim((0,offset+1))
                ax[i][j].tick_params(axis='y',labelleft='off')

    if savefig:
        plt.savefig(filename+'.pdf',format='pdf',overwrite=True)
    if showfig:
        plt.show()
    return(ax,fig)





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

            if lcs.images[lc].ml is not None:
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

def display(lclist=[], splist=[],
            title=None, style=None, showlegend=True, showlogo=False, logopos="left", showdates=False, showdelays=False,
            nicefont=False, text=None, keeponlygrid=False,
            jdrange=None, magrange=None, figsize=(12, 8), plotsize=(0.08, 0.96, 0.09, 0.95), showgrid=False,
            markersize=6, showerrorbars=True, showdatapoints=True, errorbarcolour="#BBBBBB", capsize=3, knotsize=0.015,
            legendloc="best", showspldp=False, colourprop=None, hidecolourbar=False, transparent=False,
            collapseref=False, jdmintickstep=100, magmintickstep=0.2, filename="screen", showinsert=None,
            insertname=None, verbose=False, ax=None):
    """
    Function that uses matplotlib to plot a **list** of lightcurves/splines/GPRs, either on screen or into a file.
    It uses lightcurve attributes such as ``lc.plotcolour`` and ``lc.showlabels``, and displays masked points
    as black circles. It's certainly a key function of pycs.
    You can also put tuples (lightcurve, listofseasons) in the lclist, and the seasons will be drawn.
    This function is intended both for interactive exploration and for producing publication plots.

    :param lclist: A list of lightcurve objects [lc1, lc2, ...] you want to plot.
    :type lclist: list

    :param splist: A list of spline or rslc (e.g., GPR) objects to display on top of the data points.
    :type splist: list

    :param title: Adds a title to the plot, center top of the axes, usually used for lens names.
        To nicely write a lens name, remember to use raw strings and LaTeX mathrm, e.g. :
        ``title = r"$\mathrm{RX\,J1131-1231}$"``
    :type title: string

    :param style: A shortcut to produce specific kinds of stylings for the plots.
        Available styles:

            * ``homepagepdf`` : for cosmograil homepage, ok also with long magnitude labels (like -13.2)

    :type style: string

    :param showlegend: Automatic legend (too technical/ugly for publication plots, uses str(lightcurve))
    :type showlegend: boolean

    :param showlogo: Adds an EPFL logo + www.cosmograil.org on the plot.
    :type showlogo: boolean

    :param logopos: Where to put it, "left" or "right"
    :type logopos: string

    :param showdates: If True, the upper x axis will show years and 12 month minor ticks.
    :type showdates: boolean

    :param showdelays: If True, the relative delays between the curves are written on the plot.
    :type showdelays: boolean

    :param nicefont: Sets default to serif fonts (terrible implementation, but works)
    :type nicefont: boolean

    :param text:
        Generic text that you want to display on top of the plot, in the form : [line1, line2, line3 ...]
        where line_i is (x, y, text, kwargs) where kwargs is e.g. {"fontsize":18} and x and y are relative positions (from 0 to 1).
    :type text: list

    :param jdrange: Range of jds to plot, e.g. (53000, 54000).
    :type jdrange: tuple

    :param magrange: Range of magnitudes to plot, like ``magrange = [-11, -13]``.
        If you give only a float, like ``magrange=1.0``, I'll plot +/- this number around the mean curve level
        Default is None -> automatic mag range.
    :type magrange: tuple

    :param figsize: Figure size (width, height) in inches, for interactive display or savefig.
    :type figsize: tuple

    :param plotsize: Position of the axes in the figure (left, right, bottom, top), default is (0.065, 0.96, 0.09, 0.95).
    :type plotsize: tuple

    :param showgrid: Show grid, that is vertical lines only, one for every year.
    :type showgrid: boolean

    :param markersize: Size of the data points, default is 6
    :type markersize: float

    :param showerrorbars: If False, the ploterrorbar settings of the lightcurves are disregared and no error bars are shown.
    :type showerrorbars: boolean

    :param showdatapoints: If False, no data points are shown. Useful if you want e.g. to plot only the microlensing
    :type showerrorbars: boolean

    :param keeponlygrid: If True, keeps the yearly grid from showdates but do not display the dates above the plot.
    :type keeponlygrid: boolean

    :param showinsert: If True, display the insertname image in the top-right corner of the main image
    :type showinsert: boolean

    :param insertname: path to the image you want to insert
    :type instername: string

    :param errorbarcolour: Color for the error bars
    :type errorbarcolour: string

    :param capsize: Size of the error bar "ticks"
    :type capsize: float

    :param knotsize: For splines, the length of the knot ticks, in magnitudes.
    :type knotsize: float

    :param legendloc: Position of the legend. It is passed to matplotlib legend. It can be useful to specify this if you want to make animations etc.

            * string	int
            * upper right	1
            * upper left	2
            * lower left	3
            * lower right	4
            * right	5


            * center left	6
            * center right	7
            * lower center	8
            * upper center	9
            * center	10

    :type legendloc: string or int

    :param showspldp: Show the acutal data points of spline objects (for debugging etc)
    :type showspldp: boolean

    :param colourprop: If this is set I will make a scatter plot with points coloured according to the given property, disregarding the lightcurve attribute plotcolour.
        Format : (property_name, display_name, minval, maxval), where display_name is a "nice" version of property_name, like "FWHM [arcsec]" instead of "seeing".
        Note that this property will be used in terms of a float. So you cannot use this for properties that are not floats.
    :type colourprop: tuple

    :param hidecolourbar: Set to True to hide the colourbar for the colourprop
    :type hidecolourbar: boolean

    :param transparent: Set transparency for the plot, if saved using filename
    :type transparent: boolean

    :param collapseref: Plot one single dashed line as the reference for the microlensing.
        Use this if you would otherwise get ugly overplotted dashed lines nearly at the same level ...
        This option is a bit ugly, as it does not correct the actual microlensing curves for the collapse.
    :type collapseref: boolean

    :param jdmintickstep: Minor tick step for jd axis
    :type jdmintickstep: float

    :param magmintickstep: Minor tick step for mag axis
    :type magmintickstep: float

    :param filename: If this is not "screen", I will save the plot to a file instead of displaying it. Try e.g. "test.png" or "test.pdf". Success depends on your matplotlib backend.
    :type filename: string

    :param ax: if not None, I will return what I plot in the given matplotlib axe you provide me with instead of plotting it.
    :type ax: matplotlib axes

    :param verbose: Set to True if you want me to print some details while I process the curves
    :type verbose: boolean


    """

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    import matplotlib.dates
    import matplotlib.lines

    if style == None:
        pass
    elif style == "homepagepdf":
        figsize = (10, 5)
        plotsize = (0.09, 0.97, 0.10, 0.95)
        showlogo = False
        nicefont = False
        showdelays = False
        showlegend = False
        showdates = True
        errorbarcolour = "#777777"
        markersize = 3.0
        capsize = 0
        jdmintickstep = 50
        magmintickstep = 0.2
        showgrid = True
        transparent = False

    elif style == "posterpdf":
        figsize = (10, 5.5)
        plotsize = (0.09, 0.97, 0.10, 0.95)
        showlogo = False
        nicefont = True
        showdelays = False
        showlegend = True
        showdates = True
        errorbarcolour = "#777777"
        markersize = 10.0
        capsize = 0
        jdmintickstep = 50
        magmintickstep = 0.2
        showgrid = False
        transparent = True
        title = None
        fontsize=18


    elif style == "internal":
        figsize = (10, 5.5)
        plotsize = (0.09, 0.97, 0.10, 0.95)
        showlogo = False
        nicefont = True
        showdelays = False
        showlegend = True
        showdates = True
        errorbarcolour = "#777777"
        markersize = 5.0
        capsize = 0
        jdmintickstep = 50
        magmintickstep = 0.2
        showgrid = True
        transparent = False



    else:
        raise RuntimeError("I do not know the style %s" % (style))

    if not (isinstance(lclist, list) or isinstance(lclist, tuple)):
        raise TypeError, "Hey, give me a LIST of lightcurves !"

    if colourprop != None:
        (colourpropname, colournicename, colourminval, colourmaxval) = colourprop

    labelfontsize = 18
    if nicefont:
        # mpl.rcParams['font.size'] = 20
        mpl.rcParams['font.family'] = 'serif'
    # labelfontsize = 20
    else:
        labelfontsize = 14

    if ax == None:
        fig = plt.figure(figsize=figsize)  # sets figure size
        fig.subplots_adjust(left=plotsize[0], right=plotsize[1], bottom=plotsize[2], top=plotsize[3])
        axes = plt.gca()

    else:
        ihaveax = True
        axes = ax

    if verbose: print "Plotting %i lightcurves and %i splines ..." % (len(lclist), len(splist))

    reflevels = []  # only used for collapseref

    # The lightcurves :
    for curve in lclist:
        if showdatapoints:
            if type(curve).__name__ == 'tuple':  # then we have both a lightcurve and a season to plot

                actualcurve = curve[0]
                curveseasons = curve[1]

                if not isinstance(curveseasons, list):
                    raise TypeError, "lc.display wants LISTs of seasons, not individual seasons !"
                for curveseason in curveseasons:
                    # the x lims :
                    (x1, x2) = curveseason.getjdlims(actualcurve)
                    # for y, lets take the median of that season
                    y = np.median(actualcurve.getmags()[curveseason.indices])

                    # we make this robust even with old versions of matplotlib, so no fancy arrows here.
                    axes.plot([x1, x2], [y, y], color=actualcurve.plotcolour, dashes=(1, 1))
                    axes.annotate(str(curveseason), ((x1 + x2) / 2.0, y), xytext=(-50, -15), textcoords='offset points',
                                  size=10, color=actualcurve.plotcolour)
                # plt.axvline(seasonjdlims[0], color = curve[0].plotcolour, dashes = (5,5))
                # plt.axvline(seasonjdlims[1], color = curve[0].plotcolour, dashes = (5,5))

                curve = curve[0]  # for the rest of this loop, curve is now only the lightcurve.

            if verbose: print "#   %s -> %s\n\t%s" % (curve, str(curve.plotcolour), "\n\t".join(curve.commentlist))
            # if verbose and (curve.ml != None):
            #	print curve.ml.longinfo()

            tmpjds = curve.getjds()
            tmpmags = curve.getmags()  # to avoid calculating the microlensing each time we need it

            if colourprop != None:
                scattervalues = np.array([float(propertydict[colourpropname]) for propertydict in curve.properties])
                axes.scatter(tmpjds, tmpmags, s=markersize, c=scattervalues, vmin=colourminval, vmax=colourmaxval,
                             edgecolors="None",alpha=.1)
            else:

                if curve.ploterrorbars and showerrorbars:
                    axes.errorbar(tmpjds, tmpmags, curve.magerrs, fmt=".", markersize=markersize,
                                  markeredgecolor=curve.plotcolour, color=curve.plotcolour, ecolor=errorbarcolour,
                                  capsize=capsize, label=str(curve), elinewidth=0.5)
                # plt.errorbar(tmpjds, tmpmags, curve.magerrs, linestyle="-", marker=".", color=curve.plotcolour, ecolor="#BBBBBB", label=str(curve))
                #
                else:
                    axes.plot(tmpjds, tmpmags, marker=".", markersize=markersize, linestyle="None",
                              markeredgecolor=curve.plotcolour, color=curve.plotcolour, label=str(curve))

            # We plot little round circles around masked points.
            axes.plot(tmpjds[curve.mask == False], tmpmags[curve.mask == False], linestyle="None", marker="o",
                      markersize=8., markeredgecolor="black", markerfacecolor="None", color="black")

        # And now we want to graphically display the microlensing in a nice way. This costs some cpu but anyway
        # for a display it's fine.
        # if curve.ml != None and curve.hideml == False:
        '''
        if curve.ml != None:
            if curve.ml.mltype in ["poly", "leg"]:
                for sfct in curve.ml.mllist:
                    smoothml = sfct.smooth(curve)
                    if not collapseref:
                        axes.plot(smoothml['jds'], smoothml['refmags'], color=curve.plotcolour,
                                  dashes=((3, 3)))  # the new ref
                    else:
                        reflevels.append(np.mean(smoothml['refmags']))
                    axes.plot(smoothml['jds'], smoothml['refmags'] + smoothml['ml'], color=curve.plotcolour)
            if curve.ml.mltype == "spline":
                smoothml = curve.ml.smooth(curve)

                if not collapseref:
                    axes.plot(smoothml['jds'], np.zeros(smoothml["n"]) + smoothml['refmag'], color=curve.plotcolour,
                              dashes=((3, 3)))  # the new ref
                else:
                    reflevels.append(smoothml['refmag'])

                if hasattr(curve, "hideml"):
                    if not curve.hideml:
                        axes.plot(smoothml['jds'], smoothml['refmag'] + smoothml['ml'], color=curve.plotcolour)
                else:
                    axes.plot(smoothml['jds'], smoothml['refmag'] + smoothml['ml'], color=curve.plotcolour)
                # We want to overplot the knots
                if hasattr(curve, "hideml"):
                    if not curve.hideml:
                        if getattr(curve.ml.spline, "showknots", True) == True:
                            axes.errorbar(smoothml['knotjds'], smoothml['knotmags'] + smoothml["refmag"],
                                          knotsize * np.ones(len(smoothml['knotjds'])), capsize=0,
                                          ecolor=curve.plotcolour, linestyle="none", marker="", elinewidth=1.5)
                else:
                    if getattr(curve.ml.spline, "showknots", True) == True:
                        axes.errorbar(smoothml['knotjds'], smoothml['knotmags'] + smoothml["refmag"],
                                      knotsize * np.ones(len(smoothml['knotjds'])), capsize=0, ecolor=curve.plotcolour,
                                      linestyle="none", marker="", elinewidth=1.5)
        '''
        # Labels if wanted :
        if curve.showlabels:
            for i, label in enumerate(curve.labels):
                if label != "":
                    # axes.annotate(label, (curve.jds[i], curve.mags[i]))
                    if len(label) > 4:  # Probably jd labels, we write vertically :
                        axes.annotate(label, (tmpjds[i], tmpmags[i]), xytext=(-3, -70), textcoords='offset points',
                                      size=12, color=curve.plotcolour, rotation=90)
                    else:  # horizontal writing
                        axes.annotate(label, (tmpjds[i], tmpmags[i]), xytext=(7, -6), textcoords='offset points',
                                      size=12, color=curve.plotcolour)

    if collapseref and len(reflevels) != 0:
        print "WARNING : collapsing the refs %s" % (reflevels)
        axes.axhline(np.mean(np.array(reflevels)), color="gray", dashes=((3, 3)))  # the new ref

    # The supplementary objects
    if len(splist) != 0:
        for stuff in splist:

            # We do some stupid type checking. But I like this as it does not require
            # to import spline and gpr etc.
            if hasattr(stuff, "knottype"):  # Then it's a spline
                spline = stuff
                if verbose: print "#   %s -> %s" % (str(spline), str(spline.plotcolour))

                npts = (spline.datapoints.jds[-1] - spline.datapoints.jds[0]) * 2.0
                xs = np.linspace(spline.datapoints.jds[0], spline.datapoints.jds[-1], npts)
                ys = spline.eval(jds=xs)
                axes.plot(xs[575:-600], ys[575:-600], "-", color=spline.plotcolour, zorder=+20, label=str(spline))
                # For the knots, we might not want to show them (by default we do show them) :
                if getattr(spline, "showknots", True) == True:
                    if ax != None:
                        ax.errorbar(spline.getinttex(), spline.eval(jds=spline.getinttex()),
                                    0.015 * np.ones(len(spline.getinttex())), capsize=0, ecolor=spline.plotcolour,
                                    linestyle="none", marker="", elinewidth=1.5, zorder=40, barsabove=True)
                        knotxs = spline.getinttex()
                        knotys = spline.eval(knotxs)
                        for (knotx, knoty) in zip(knotxs, knotys):
                            l = matplotlib.lines.Line2D([knotx, knotx], [knoty - knotsize, knoty + knotsize], zorder=30,
                                                        linewidth=1.5, color=spline.plotcolour)
                            ax.add_line(l)
                    else:
                        axes = plt.gca()
                        knotxs = spline.getinttex()
                        knotys = spline.eval(knotxs)
                        for (knotx, knoty) in zip(knotxs, knotys):
                            l = matplotlib.lines.Line2D([knotx, knotx], [knoty - knotsize, knoty + knotsize], zorder=30,
                                                        linewidth=1.5, color=spline.plotcolour)
                            axes.add_line(l)

                if showspldp:  # The datapoints of the spline (usually not shown)
                    axes.plot(spline.datapoints.jds, spline.datapoints.mags, marker=",", linestyle="None",
                              color=spline.plotcolour, zorder=-20)

                    #			if hasattr(stuff, "regfct"): # Then it's a GPR
                    #
                    # 				gpr = stuff
                    # 				npts = (gpr.jds[-1] - gpr.jds[0])*2.0
                    # 				xs = np.linspace(gpr.jds[0], gpr.jds[-1], npts)
                    # 				(ys, yerrs) = gpr.regfct(xs)
                    # 				#print "regfct evaluated"
                    # 				plt.plot(xs, ys, "-", color=gpr.plotcolour, zorder=+20, label=str(gpr))
                    # 				xf = np.concatenate((xs, xs[::-1]))
                    #         			yf = np.concatenate((ys+yerrs, (ys-yerrs)[::-1]))
                    #         			plt.fill(xf, yf, facecolor = gpr.plotcolour, alpha=0.1, edgecolor = (1,1,1))

            if hasattr(stuff, "pd"):  # Then it's a rslc

                rs = stuff
                # plt.plot(rs.getjds(), rs.mags, "-.", color=rs.plotcolour)
                axes.plot(rs.getjds(), rs.mags, "-", color=rs.plotcolour)
                xf = np.concatenate((rs.getjds(), rs.getjds()[::-1]))
                yf = np.concatenate((rs.mags + rs.magerrs, (rs.mags - rs.magerrs)[::-1]))
                plt.fill(xf, yf, facecolor=rs.plotcolour, alpha=0.2, edgecolor=(1, 1, 1), label=str(rs))

    # Astronomers like minor tick marks :
    minorxLocator = MultipleLocator(jdmintickstep)
    axes.xaxis.set_minor_locator(minorxLocator)
    # minorLocator = MultipleLocator(1) # so to have a tick every day
    # axes.xaxis.set_minor_locator(minorLocator)


    # Something for astronomers only : we invert the y axis direction !
    axes.set_ylim(axes.get_ylim()[::-1])

    if colourprop != None and hidecolourbar == False:
        cbar = plt.colorbar(orientation='vertical', shrink=1.0, fraction=0.065, pad=0.025)
        cbar.set_label(colournicename)

    # And we make custom title :

    if title == "None" or title == None or title == "none":
        # plt.title("Lightcurves", fontsize=18)
        pass
    else:
        # plt.title(title, fontsize=18)
        axes.annotate(title, xy=(0.5, 1.0), xycoords='axes fraction', xytext=(0, -4),
                      textcoords='offset points', ha='center', va='top', fontsize=25)

    if jdrange != None:
        axes.set_xlim(jdrange[0], jdrange[1])

    axes.set_xlabel("MJD [day]", fontsize=labelfontsize)
    axes.set_ylabel("AB Magnitude", fontsize=labelfontsize)

    if showdelays:
        txt = getnicetimedelays(lclist, separator="\n")
        axes.annotate(txt, xy=(0.0, 1.0), xycoords='axes fraction', xytext=(6, -6),
                      textcoords='offset points', ha='left', va='top')
        # plt.text(0.01, 0.99, txt,
        #	horizontalalignment='left', verticalalignment='top',
        #	transform = axes.transAxes)
        legendloc = 1
        if verbose:
            print "Delays between plotted curves :"
            print txt
    '''
    #micro
    if showlegend and (len(lclist) > 0 or len(splist) > 0):
        font = {'family': 'serif',
                'color': 'darkred',
                'weight': 'normal',
                'size': 18,
                }
        #plt.title("Multi-Knot Spline Fitting (Microlensing Included)",fontsize=20)
        plt.text(57295,24,r"Multi-Knot Spline Fitting",fontdict=font)
        plt.text(57301,24.5,r"(Microlensing Included)",fontdict=font)
        mpl.rcParams['legend.numpoints'] = 1
        mpl.rcParams['legend.loc'] = 'lower left'
        h,l=axes.get_legend_handles_labels()
        l=[x[8:16]+' (Delay='+str(format(abs(float(x[18:24])-delay),'.3f'))+', Shift='+str(format(float(x[25:30]),'.3f'))+')' for x in l[1:]]
        l=['Best Fit Spline']+l
        axes.legend(loc=3, prop=fm.FontProperties(size=12))
        axes.legend(h,l,title='SN Refsdal')

    #no micro
    if showlegend and (len(lclist) > 0 or len(splist) > 0):
        font = {'family': 'serif',
                'color': 'darkred',
                'weight': 'normal',
                'size': 18,
                }
        plt.text(57295,24.2,r"Simple Spline Fitting",fontdict=font)
        plt.text(57308,24.6,r"(No Microlensing)",fontdict=font)
        #plt.title("Simple Spline Fitting (No Microlensing)", fontsize=20)
        mpl.rcParams['legend.numpoints'] = 1
        mpl.rcParams['legend.loc'] = 'lower left'
        h, l = axes.get_legend_handles_labels()
        l = [x[8:16] + ' (Delay=' + str(format(abs(float(x[18:24]) - delay), '.3f')) + ', Shift=' + str(format(float(x[26:31]),'.3f')) + ')' for x in l[1:]]
        l = ['Simple Spline'] + l
        axes.legend(loc=3, prop=fm.FontProperties(size=12))
        axes.legend(h, l, title='SN Refsdal')
    '''
    #points
    if showlegend and (len(lclist) > 0 or len(splist) > 0):
        font = {'family': 'serif',
                'color': 'darkred',
                'weight': 'normal',
                'size': 18,
                }
        plt.text(57358, 24.7, r"SN Refsdal:", fontdict=font)
        plt.text(57305, 25.1, r"F160W Unshifted Data", fontdict=font)
        #plt.title("SN Refsdal: F160W Unshifted Data", fontsize=20)
        mpl.rcParams['legend.numpoints'] = 1
        mpl.rcParams['legend.loc'] = 'lower left'
        h, l = axes.get_legend_handles_labels()
        l = [x[8:16] for x in l]
        axes.legend(loc=3, prop=fm.FontProperties(size=12))
        axes.legend(h, l, title='SN Refsdal')

    if magrange != None:
        if type(magrange) == float or type(magrange) == int:
            # We find the mean mag of the stuff to plot :
            allmags = []
            for l in lclist:
                allmags.extend(l.getmags())
            meanlevel = np.mean(np.array(allmags))
            axes.set_ylim(meanlevel + magrange, meanlevel - magrange)
        else:
            axes.set_ylim(magrange[0], magrange[1])

    showdates=False
    if showdates:  # Be careful when you change something here, it could mess up the axes.
        # Especially watch out when you change the plot range.
        # This showdates stuff should come at the very end
        minjd = axes.get_xlim()[0]
        maxjd = axes.get_xlim()[1]
        # axes.set_xlim(minjd, maxjd)
        yearx = axes.twiny()
        yearxmin = util.datetimefromjd(minjd + 2400000.5)
        yearxmax = util.datetimefromjd(maxjd + 2400000.5)
        yearx.set_xlim(yearxmin, yearxmax)
        yearx.xaxis.set_minor_locator(matplotlib.dates.MonthLocator())
        yearx.xaxis.set_major_locator(matplotlib.dates.YearLocator())
        yearx.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
        yearx.xaxis.tick_top()
        yearx.tick_params(direction='in',pad=-20)
        if keeponlygrid:
            yearx.set_xticklabels([])
            # yearx.set_xlabel("Date")

    minoryLocator = MultipleLocator(magmintickstep)
    axes.yaxis.set_minor_locator(minoryLocator)

    if showgrid:
        plt.grid(zorder=20)

    if text != None:
        for line in text:
            axes.text(line[0], line[1], line[2], transform=axes.transAxes, **line[3])

    if showinsert:
        assert insertname != None
        from matplotlib._png import read_png
        from matplotlib.offsetbox import OffsetImage, AnnotationBbox
        im = read_png(insertname)
        imagebox = OffsetImage(im, zoom=0.5, interpolation="sinc", resample=True)
        ab = AnnotationBbox(imagebox, xy=(1.0, 1.0), xycoords='axes fraction', xybox=(-75, -75),
                            boxcoords="offset points",
                            pad=0.0, frameon=False
                            )
        axes.add_artist(ab)

    if showlogo:

        # The EPFL logo :
        from matplotlib._png import read_png
        from matplotlib.offsetbox import OffsetImage, AnnotationBbox
        logodir = os.path.dirname(__file__)
        im = read_png(os.path.join(logodir, "epfl.png"))
        imagebox = OffsetImage(im, zoom=0.13, interpolation="sinc", resample=True)

        if logopos == "left":
            ab = AnnotationBbox(imagebox, xy=(0.0, 0.0), xycoords='axes pixels', xybox=(52, 30),
                                boxcoords="offset points",
                                pad=0.0, frameon=False
                                )
            axes.add_artist(ab)
            axes.annotate("COSMOGRAIL.org", xy=(0.0, 0.0), xycoords='axes fraction', fontsize=16, xytext=(105, 7),
                          textcoords='offset points', ha='left', va='bottom', color="gray")

            ## name lightcurves:
            if 0:
                axes.annotate("A", xy=(0.0, 0.0), xycoords='axes fraction', fontsize=25, xytext=(20, 260),
                              textcoords='offset points', ha='center', va='bottom', color="red")
                axes.annotate("B", xy=(0.0, 0.0), xycoords='axes fraction', fontsize=25, xytext=(20, 150),
                              textcoords='offset points', ha='center', va='bottom', color="green")
                axes.annotate("C", xy=(0.0, 0.0), xycoords='axes fraction', fontsize=25, xytext=(20, 125),
                              textcoords='offset points', ha='center', va='bottom', color="blue")
                axes.annotate("D", xy=(0.0, 0.0), xycoords='axes fraction', fontsize=25, xytext=(20, 80),
                              textcoords='offset points', ha='center', va='bottom', color="purple")

        if logopos == "right":
            ab = AnnotationBbox(imagebox, xy=(1.0, 0.0), xycoords='axes fraction', xybox=(-200, 30),
                                boxcoords="offset points",
                                pad=0.0, frameon=False
                                )
            axes.add_artist(ab)
            axes.annotate("COSMOGRAIL.org", xy=(1.0, 0.0), xycoords='axes fraction', fontsize=16, xytext=(-10, 7),
                          textcoords='offset points', ha='right', va='bottom', color="gray")

        if logopos == "center":
            ab = AnnotationBbox(imagebox, xy=(0.55, 0.0), xycoords='axes fraction', xybox=(-80, 30),
                                boxcoords="offset points",
                                pad=0.0, frameon=False
                                )
            axes.add_artist(ab)
            axes.annotate("COSMOGRAIL.org", xy=(0.55, 0.0), xycoords='axes fraction', fontsize=16, xytext=(40, 7),
                          textcoords='offset points', ha='center', va='bottom', color="gray")

        # Alternative possibility (just to keep the idea) :
        """
        try:
            import Image
        except ImportError:
            print "Couldn't import PIL ! Therefore I won't be able to display the cosmograil logo."
        else:
            im = Image.open('epfl.png')
            height = im.size[1]
            print height
            im = np.array(im).astype(np.float) / 255
            fig.figimage(im, 0, fig.bbox.ymax - height)
            # With newer (1.0) versions of matplotlib, you can
            # use the "zorder" kwarg to make the image overlay
            # the plot, rather than hide behind it... (e.g. zorder=10)
            fig.figimage(im, 0, fig.bbox.ymax - height)
        """

    if ax != None:
        return

    if filename == "screen":
        plt.show()
    else:

        plt.savefig(filename, transparent=transparent,format='pdf')
        # if verbose:
        print "Plot written to %s" % filename
        #plt.show()  # this seems important so that the plot is not displayed when a next plt.show() is called.
        plt.clf()



