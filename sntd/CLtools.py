"""
Dan Coe (2009)
http://arxiv.org/abs/0906.4123
Fisher Matrices and Confidence Ellipses: A Quick Start Guide and Software

used in Coe & Moustakas (2009):
http://arxiv.org/abs/0906.4108
Cosmological Constraints from Gravitational Lens Time Delays
"""

from .coetools import *
from .coeplot import *
from os.path import join, exists


# Ellipse colors
# outer (lighter), inner (darker)
blues = [(0,1,1), (0,0,1)]
#reds = [(1,0,1), (1,0,0)]
reds = [(1,0.5,0.5), (1,0,0)]
purples = [(1,0,1), (0.5,0,0.5)]
yellows = [(1,1,0), (0.5,0.5,0)]
oranges = [(1,0.7,0.5), (1,0.5,0)]
darkoranges = [(1,0.5,0), (0.5,0.25,0)]
#greens = [(0,1,0), (0,0.8,0)]
#greens = [(0.5,1,0.5), (0,1,0)]
greens = [(0.25,1,0.25), (0,0.8,0)]
greys = [(0.5,0.5,0.5), (0,0,0)]
lightblues = [(0.7,1,1), (0.7,0.7,1)]
lightreds = [(1,0.7,1), (1,0.7,0.7)]
lightgreens = [(0.7,1,0.7), (0.5,0.8,0.5)]
lightyellows = [(1,1,0.7), (0.8,0.8,0.5)]
lightgreys = [(0.8,0.8,0.8), (0.5,0.5,0.5)]

# Axis labels
axlabs = {}
axlabs['om'] = '$\omega_m = \Omega_m h^2$'
axlabs['Om'] = '$\Omega_m$'
axlabs['OQ'] = '$\Omega_\Lambda$'
axlabs['OL'] = '$\Omega_\Lambda$'
axlabs['Ode'] = '$\Omega_{de}$'
axlabs['Ok'] = '$\Omega_k$'
axlabs['w']  = 'w'
axlabs['h']  = '$h$'
axlabs['w0']  = r'${\rm w}_0$'
axlabs['wa']  = r'${\rm w}_a$'

# List of variables for which axis labels should be 
# horizontal rather than vertical
axhor = 'h w w0 wa OQ OL Ode Om Ok'.split()

ndec = 2
def mapfmt(x):
    fmt = '%%.%df' % ndec
    return fmt % x

# ~/cosmo/DETFast1/plots/flat/plotell.py
def rightaxis(sh=1, ndeca=2, adj=0):
    """Plot right axis: Om = 1 - Ode"""
    # Right axis
    # e.g., ~/LensPerfect/A1689/analysis/NFWfitsplot.py
    #CLtools.ndec = ndec
    global ndec
    ndec = ndeca
    ax1 = gca()
    ax1.yaxis.tick_left()
    ax2 = gcf().add_axes(ax1.get_position(), frameon=False)
    ylim(ax1.get_ylim())
    ax2.set_xticks([])
    ax2.yaxis.tick_right()
    ytx = ax2.get_yticks()
    ytx = 1 - ytx
    ytxs = list(map(mapfmt, ytx))
    p = ax2.set_yticklabels(ytxs)
    ax2.set_ylabel('$\Omega_m$', rotation='horizontal')
    ax2.yaxis.set_label_position('right')
    if adj: subplots_adjust(left=0.1, right=0.9)  # default: 0.125, 0.9
    if sh: show()

def xvarlabel(xvar):
    """Mark x-axis label, given abbreviation xvar (e.g., 'w0')"""
    xlabel(axlabs[xvar])

def yvarlabel(yvar):
    """Mark x-axis label, given abbreviation xvar (e.g., 'wa')"""
    rot = ['vertical', 'horizontal'][yvar in axhor]
    ylabel(axlabs[yvar], rotation=rot)

def xylabels(xvar, yvar):
    """Mark x-axis and y-axis labels, given abbreviations (e.g., 'w0', 'wa')"""
    xlabel(axlabs[xvar])
    rot = ['vertical', 'horizontal'][yvar in axhor]
    ylabel(axlabs[yvar], rotation=rot)

def markinput(x0, y0, c='0.40', dc='w'):
    """Dotted gray lines and white circle (default)
    mark best fit input parameters x0, y0
    e.g., (w0, wa) = (-1, 0)"""
    xr = xlim()
    yr = ylim()
    axlines(x0, y0, c=c, ls=':', zorder=20)
    if dc != None: 
        dcf = dce = dc
        if type(dc) == str:
            if len(dc) > 1:
                dcf, dce = dc
        plot([x0], [y0], 'o', mfc=dcf, mec=dce, zorder=50)
    xlim(xr)
    ylim(yr)

def finishup(x0, y0, xvar, yvar, c='0.40', dc='w', addax=0, sh=1, ndec=2):
    """Mark best fit cosmo, add axis labels, etc."""
    markinput(x0, y0, c=c, dc=dc)
    xylabels(xvar, yvar)
    if addax: rightaxis(sh, ndec)
    elif sh: show()

def setnd(nd=2):
    """Set sigma for confidence contours
    nd = # dimensions (# parameters)
    Normally this should be 2
    unless there is a *total* degeneracy, unless it should be 1
    """
    global nsig, nsig2, nsig1, nsig0
    if nd == 2:
        nsig = sqrt(array([6.17, 2.3, 0]))
    elif nd == 1:
        nsig = array([2, 1, 0])
    nsig2, nsig1, nsig0 = nsig
    return nsig

setnd()  # DEFAULT 2-D
def chisqplot(chisq, xs, ys, cl=0, addax=0, colors=[(0,1,1), (0,0,1)], zorder=None, alpha=1, nlev=2, nD=2, sh=1, fill=1, lw=2):  # ls='-'
    """Plot confidence contours given chisq(xs,ys) (matrix[y,x])
    Sigma:
    (0, 1, 2) -- default
    (0, 2) -- nlev = 1
    (0, 1) -- nlev = -1
    To plot 3-sigma contours or more, 
    colors would need to be given additional color(s), and other tweaks?"""
    global CS, CS2
    colors = colors[::-1] # dark (inner), light (outer)
    chisq = chisq - min(chisq.flat)
    
    #nsig = 
    #dchisq = array([0, 2.3, 6.17])  # for ellipse, yields 95.4% & 68% CL (covariance.py)
    #if nD == 1: dchisq = array([0, 1, 2])  # 1-D
    #if nD == 1: dchisq = array([0, 1, 4])  # 1-D
    setnd(nD)
    dchisq = nsig**2
    if nlev == 1: dchisq = dchisq.take([0, 2])
    if nlev == -1: dchisq = dchisq.take([0, 1])
    
    if cl: clf()
    if fill:
        CS = contourf(xs, ys, chisq, dchisq, colors=colors[:nlev], zorder=zorder)
    else:
        CS = contour(xs, ys, chisq, dchisq, colors=colors[:nlev][::-1], zorder=zorder, linewidths=[lw])
    if (alpha < 1) or (zorder != None): 
        contour_set(CS, alpha=alpha, zorder=zorder, sh=sh)
    color = colors[0]
    if type(color) in (tuple, list): color = [color]
    dchisq = dchisq.take((0,-1))
    if fill:
        CS2 = contour(xs, ys, chisq, dchisq, colors=color, hold='on', linewidths=[lw], zorder=zorder, alpha=alpha)
    if (alpha < 1) or (zorder != None): 
        if zorder != None: zorder += 0.01  # Lines on top of solid regions
        # Default: lines = 2, solids = 1
        contour_set(CS2, alpha=alpha, zorder=zorder, color=color[0], sh=sh)
    
    if addax: rightaxis()

def Pplot(P, xs, ys, **other):
    """Plot confidence contours given probability P(xs,ys) (matrix[y,x])"""
    chisq = -2 * log(P)
    chisqplot(chisq, xs, ys, **other)

def setell(dx, dy, p, nD=2):
    """Return ellipse parameters A, B, ang[degrees] given dx, dy, p
    dx, dy = 1-sigma uncertainties in x, y
    p = correlation coefficient (0 = independent, 1 = completely correlated)
    Eqs. 2,3,4 from Fisher Quick Start Guide"""
    # Ellipse formulae
    # http://www.scribd.com/doc/2341150/Lecture-5-Bivariate-ND-Error-Ellipses
    dxy = p * dx * dy
    
    de1 = (dx**2 + dy**2) / 2
    de2 = (dx**2 - dy**2)**2 / 4
    de3 = sqrt(de2 + dxy**2)
    A = sqrt(de1 + de3)
    B = sqrt(de1 - de3)
    if nD == 1:
        A = A * sqrt(2)
        B = B * sqrt(2)
    
    sin2ang = 2 * dxy # / something
    cos2ang = dx**2 - dy**2 # / something
    ang = atanxy(cos2ang, sin2ang) / 2
    return A, B, ang

def plotells(xo, yo, A, B, ang, colors=blues, alpha=1, zorder=1, lw=1, fill=1, nD=2):
    """Plot 1-sigma & 2-sigma ellipses given A, B, ang[radians]"""
    setnd(nD)
    patch = ellpatch1(xo, yo, 2*nsig2*A, 2*nsig2*B, ang, 
                      fc=colors[0], ec=colors[1], alpha=alpha, zorder=zorder, fill=fill)
    patch.set_lw(lw + fill)
    gca().add_patch(patch)
    
    patch = ellpatch1(xo, yo, 2*nsig1*A, 2*nsig1*B, ang, 
                      fc=colors[1], ec=colors[1], alpha=alpha, zorder=zorder, fill=fill)
    patch.set_lw(lw + 1 - fill)
    gca().add_patch(patch)

def plotellsp(xo, yo, dx, dy, p, colors=blues, alpha=1, zorder=1, lw=1, fill=1, nD=2):
    """Plot 1-sigma & 2-sigma ellipses given dx, dy, p
    dx, dy = 1-sigma uncertainties in x, y
    p = correlation coefficient (0 = independent, 1 = completely correlated)
    """
    A, B, ang = setell(dx, dy, p)
    plotells(xo, yo, A, B, ang,
             colors=colors, alpha=alpha, zorder=zorder, lw=lw, fill=fill, nD=nD)

def plotell2p(xo, yo, dx, dy, p, color='b', fillcolor=None, alpha=1, zorder=1, lw=1, fill=False, nD=2):
    """Plot 2-sigma ellipse given dx, dy, p
    dx, dy = 1-sigma uncertainties in x, y
    p = correlation coefficient (0 = independent, 1 = completely correlated)
    """
    setnd(nD)
    if fill and fillcolor==None:
        fillcolor = color
    A, B, ang = setell(dx, dy, p)
    patch = ellpatch1(xo, yo, 2*nsig2*A, 2*nsig2*B, ang, 
                      fillcolor, color, alpha, zorder, fill)
    patch.set_lw(lw)
    gca().add_patch(patch)

def plotell1p(xo, yo, dx, dy, p, color='b', fillcolor=None, alpha=1, zorder=1, lw=1, fill=False, nD=2):
    """Plot 1-sigma ellipse given dx, dy, p
    dx, dy = 1-sigma uncertainties in x, y
    p = correlation coefficient (0 = independent, 1 = completely correlated)
    """
    setnd(nD)
    if fill and fillcolor==None:
        fillcolor = color
    A, B, ang = setell(dx, dy, p)
    patch = ellpatch1(xo, yo, 2*nsig1*A, 2*nsig1*B, ang, 
                      fillcolor, color, alpha, zorder, fill)
    patch.set_lw(lw)
    gca().add_patch(patch)


#################################
# ADD LEGEND

def initlegxy(xl0=0.75, yl0=0.95, sxl0=0.06, syl0=0.05):
    """Set legend coordinates"""
    global xl, yl, dxl, dyl, sxl, syl
    xr = xlim()
    yr = ylim()
    xl = interp([xl0], array([0, 1]), array(xr))[0]
    yl = interp([yl0], array([0, 1]), array(yr))[0]
    sxl = sxl0 * p2p(xr)
    syl = syl0 * p2p(yr)
    dxl = sxl * 0.8
    dyl = syl * 1.2
    
#xl, yl = 0.75, 0.535  # TOP-RIGHT (GOING DOWN)
#xl0, yl0 = 0.75, 0.95  # TOP-RIGHT

def addell(x, y, colors=blues, alpha=1, sx=0.05, sy=0.05, zorder=1, sh=0):
    """Add legend ellipse"""
    ax = gca()
    
    patch = Ellipse((x, y), sx, sy)
    ax.add_patch(patch)
    patch.set_fc(colors[0])
    patch.set_ec(colors[1])
    patch.set_alpha(alpha)
    patch.set_zorder(zorder)
    #patch.set_lw(2)
    
    #nsig = sqrt(array([6.17, 2.3, 0]))
    #fac = nsig1 / nsig2
    fac = sqrt(2.3 / 6.17)
    patch = Ellipse((x, y), fac*sx, fac*sy)
    ax.add_patch(patch)
    patch.set_fc(colors[1])
    patch.set_ec(colors[1])
    patch.set_alpha(alpha)
    patch.set_zorder(zorder)
    
    if sh: show()

xl = None

def addlab(lab, colors=blues, alpha=1, sx=0.05, sy=0.05, zorder=10, sh=0):
    """Add legend label"""
    global yl
    if xl == None: initlegxy()
    addell(xl, yl, colors, alpha, sxl, syl, zorder=zorder, sh=0)
    p = text(xl+dxl, yl, lab, va='center', zorder=1000)
    p.set_zorder(10000)
    #yl += dyl  # GOING UP
    yl -= dyl  # GOING DOWN
    if sh: show()
