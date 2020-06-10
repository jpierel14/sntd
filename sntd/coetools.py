## Automatically adapted for numpy Jun 08, 2006 by 

## Automatically adapted for numpy Jun 08, 2006 by 

# PYTHON TOOLS
# Dan Coe

import os
import sys

from numpy import *
sys.float_output_precision = 5  # PRINTING ARRAYS: # OF DECIMALS
from types import *  # TO TELL WHAT TYPE A VARIABLE IS

from time import *
import string  # LOAD THIS AFTER numpy, BECAUSE numpy HAS ITS OWN string

from numpy import *
from bisect import bisect
from scipy.integrate import quad
from scipy.special import erf
from numpy.random import *  # random
#from biggles import *
import string

# FOR MORE PLOTTING IDEAS, SEE ksbtools.py (USED biggles)

import matplotlib
#matplotlib.use('TkAgg')
from pylab import *
# I HAVE A FUNCTION close IN MLab_coe WHICH CONFLICTS WITH THE FIGURE CLOSER:

from pylab import close as closecurrentfig
from pylab import xlim as xlim1
from pylab import ylim as ylim1
import numpy.random as RandomArray
import math
#from MLab_coe import singlevalue
import os
from numpy.linalg import eig, svd
## Automatically adapted for numpy Jun 08, 2006
## By hand: 'float' -> float, Float -> float, Int -> int

# coeio.py
# INPUT / OUTPUT OF FILES

#import string

#import fitsio
try:
    import pyfits#, numarray
    pyfitsloaded = True
except:
    pyfitsloaded = False
    #pass # print "pyfits not installed, so not importing it"


from os.path import exists, join

"""
Dan Coe (2009)
http://arxiv.org/abs/0906.4123
Fisher Matrices and Confidence Ellipses: A Quick Start Guide and Software

used in Coe & Moustakas (2009):
http://arxiv.org/abs/0906.4108
Cosmological Constraints from Gravitational Lens Time Delays
"""

# ~/lp/mycolormaps.py

# ~/p/
# colormapnew.py
# colormapdata.py
# colormaps.py

# /Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/matplotlib/
# _cm.py
# pyplot.py
# colors.py

# TO SET ONE OF THESE COLORMAPS, JUST CALL, e.g.:
# >>> gray1070()
# OR:
# >>> colormap(x, y, z, cmap=cm.gray1070)

numerix = os.environ.get('NUMERIX', '')
LUTSIZE = mpl.rcParams['image.lut']

# cm.datad.keys()


#################################

lo, hi = 0.1, 0.7
x = 1  # DOESN'T MATTER??
_gray1070_data =  {
    'red':   ((0., x, lo), (1, hi, x)),
    'green': ((0., x, lo), (1, hi, x)),
    'blue':  ((0., x, lo), (1, hi, x))}

cm.gray1070 = matplotlib.colors.LinearSegmentedColormap('gray1070', _gray1070_data, LUTSIZE)

cm.datad['gray1070'] = _gray1070_data

def gray1070():
    rc('image', cmap='gray1070')
    im = gci()
    
    if im is not None:
        im.set_cmap(cm.gray1070)
    draw_if_interactive()

#################################

lo, hi = 0.2, 0.7
x = 1  # DOESN'T MATTER??
_gray2070_data =  {
    'red':   ((0., x, lo), (1, hi, x)),
    'green': ((0., x, lo), (1, hi, x)),
    'blue':  ((0., x, lo), (1, hi, x))}

cm.gray2070 = matplotlib.colors.LinearSegmentedColormap('gray2070', _gray2070_data, LUTSIZE)

cm.datad['gray2070'] = _gray2070_data

def gray2070():
    rc('image', cmap='gray2070')
    im = gci()
    
    if im is not None:
        im.set_cmap(cm.gray2070)
    draw_if_interactive()

#################################

lo, hi = 0.0, 1.0
x = 1  # DOESN'T MATTER??
_gray0010_data =  {
    'red':   ((0., x, lo), (1, hi, x)),
    'green': ((0., x, lo), (1, hi, x)),
    'blue':  ((0., x, lo), (1, hi, x))}

cm.gray0010 = matplotlib.colors.LinearSegmentedColormap('gray0010', _gray0010_data, LUTSIZE)

cm.datad['gray0010'] = _gray0010_data

def gray0010():
    rc('image', cmap='gray0010')
    im = gci()
    
    if im is not None:
        im.set_cmap(cm.gray0010)
    draw_if_interactive()

#################################

lo, hi = 1.0, 0.0
x = 1  # DOESN'T MATTER??
_grayinv_data =  {
    'red':   ((0., x, lo), (1, hi, x)),
    'green': ((0., x, lo), (1, hi, x)),
    'blue':  ((0., x, lo), (1, hi, x))}

cm.grayinv = matplotlib.colors.LinearSegmentedColormap('grayinv', _grayinv_data, LUTSIZE)

cm.datad['grayinv'] = _grayinv_data

def grayinv():
    rc('image', cmap='grayinv')
    im = gci()
    
    if im is not None:
        im.set_cmap(cm.grayinv)
    draw_if_interactive()

#################################

lo, hi = 1.0, 0.2
x = 1  # DOESN'T MATTER??
_grayinv20_data =  {
    'red':   ((0., x, lo), (1, hi, x)),
    'green': ((0., x, lo), (1, hi, x)),
    'blue':  ((0., x, lo), (1, hi, x))}

cm.grayinv20 = matplotlib.colors.LinearSegmentedColormap('grayinv20', _grayinv20_data, LUTSIZE)

cm.datad['grayinv20'] = _grayinv20_data

def grayinv20():
    rc('image', cmap='grayinv20')
    im = gci()
    
    if im is not None:
        im.set_cmap(cm.grayinv20)
    draw_if_interactive()

#################################

# Roman numerals
# I = 1
# V = 5
# X = 10
# L = 50
# C = 100
# D = 500
# M = 1000

def roman(x):
    s = ''
    n = x / 1000
    s = 'M' * n

    x = x % 1000
    n = x / 100
    if n == 9:
        s += 'CM'
    elif n == 4:
        s += 'CD'
    else:
        if n >= 5:
            s += 'D'
            n -= 5
        s += 'C' * n

    x = x % 100
    n = x / 10
    if n == 9:
        s += 'XC'
    elif n == 4:
        s += 'XL'
    else:
        if n >= 5:
            s += 'L'
            n -= 5
        s += 'X' * n

    x = x % 10
    n = x / 1
    if n == 9:
        s += 'IX'
    elif n == 4:
        s += 'IV'
    else:
        if n >= 5:
            s += 'V'
            n -= 5
        s += 'I' * n

    return s

#from numpy import compress as compress1b

# THE numpy1.0b compress SUCKS!
# I REIMPORT IT AS compress1b AND I USE IT HERE
def compress(c, m):
    #print 'compress2: NEW VERSION OF compress BECAUSE THE Numpy 1.0b VERSION IS SHWAG'
    m = array(m)
    if len(m.shape) == 1:
        mc = compress(c, m)
    elif len(m.shape) == 2:
        nm, nc = m.shape
        if type(c) != list:
            c = c.tolist()
        c = c * nm
        m = ravel(m)  # REQUIRED ON SOME MACHINES
        mc = compress(c, m)
        mc = reshape(mc, (nm, len(mc)/nm))
    else:
        print('MORE THAN 2 AXES NOT SUPPORTED BY compress2')
    return mc

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


def strspl(s):
    if type(s) == str:
        if s.find(' ') > -1:
            s = s.split()
    return s

def pint(A, n=0):
    """Makes it easier to view float arrays:
    prints A.astype(int)"""
    if type(A) in [list, tuple]:
        A = array(A)
    if n != 0:
        A = A * 10**n
    print(A.astype(int))

def pintup(A, n=0):
    """Makes it easier to view float arrays:
    prints A.astype(int)"""
    pint(flipud(A), n)

# UNLESS $NUMERIX IS SET TO numpy, pyfits(v1.1b) USES NumArray
if numerix != 'numpy':
    print('You probably should have done this first: setenv NUMERIX numpy')

#pyfitsusesnumpy = (string.atof(pyfits.__version__[:3]) >= 1.1) and (numerix == 'numpy')
#if not pyfitsusesnumpy and pyfitsloaded:
#    import numarray

def recapfile(name, ext):
    """CHANGE FILENAME EXTENSION"""
    if ext[0] != '.':
        ext = '.' + ext
    i = name.rfind(".")
    if i == -1:
        outname = name + ext
    else:
        outname = name[:i] + ext
    return outname

def capfile(name, ext):
    """ADD EXTENSION TO FILENAME IF NECESSARY"""
    if ext[0] != '.':
        ext = '.' + ext
    n = len(ext)
    if name[-n:] != ext:
        name += ext
    return name

def decapfile(name, ext=''):
    """REMOVE EXTENSION FROM FILENAME IF PRESENT
    IF ext LEFT BLANK, THEN ANY EXTENSION WILL BE REMOVED"""
    if ext:
        if ext[0] != '.':
            ext = '.' + ext
        n = len(ext)
        if name[-n:] == ext:
            name = name[:-n]
    else:
        i = name.rfind('.')
        if i > -1:
            name = name[:i]
    return name

uncapfile = decapfile


def params_cl():
    """RETURNS PARAMETERS FROM COMMAND LINE ('cl') AS DICTIONARY:
    KEYS ARE OPTIONS BEGINNING WITH '-'
    VALUES ARE WHATEVER FOLLOWS KEYS: EITHER NOTHING (''), A VALUE, OR A LIST OF VALUES
    ALL VALUES ARE CONVERTED TO INT / FLOAT WHEN APPROPRIATE"""
    list = sys.argv[:]
    i = 0
    dict = {}
    oldkey = ""
    key = ""
    list.append('')  # EXTRA ELEMENT SO WE COME BACK AND ASSIGN THE LAST VALUE
    while i < len(list):
        if striskey(list[i]) or not list[i]:  # (or LAST VALUE)
            if key:  # ASSIGN VALUES TO OLD KEY
                if value:
                    if len(value) == 1:  # LIST OF 1 ELEMENT
                        value = value[0]  # JUST ELEMENT
                dict[key] = value
            if list[i]:
                key = list[i][1:] # REMOVE LEADING '-'
                value = None
                dict[key] = value  # IN CASE THERE IS NO VALUE!
        else: # VALUE (OR HAVEN'T GOTTEN TO KEYS)
            if key: # (HAVE GOTTEN TO KEYS)
                if value:
                    value.append(str2num(list[i]))
                else:
                    value = [str2num(list[i])]
        i += 1

    return dict


def delfile(file, silent=0):
    if os.path.exists(file) or os.path.islink(file): # COULD BE BROKEN LINK!
        if not silent:
            print('REMOVING ', file, '...')
        os.remove(file)
    else:
        if not silent:
            print("CAN'T REMOVE", file, "DOES NOT EXIST.")

rmfile = delfile

def dirfile(filename, dir=""):
    """RETURN CLEAN FILENAME COMPLETE WITH PATH
    JOINS filename & dir, CHANGES ~/ TO home"""
    if filename[0:2] == '~/':
        filename = os.path.join(home, filename[2:])
    else:
        if dir[0:2] == '~/':
            dir = os.path.join(home, dir[2:])
        filename = os.path.join(dir, filename)
    return filename
    

def loadfile(filename, dir="", silent=0, keepnewlines=0):
    infile = dirfile(filename, dir)
    if not silent:
        print("Loading ", infile, "...\n")
    fin = open(infile, 'r')
    sin = fin.readlines()
    fin.close()
    if not keepnewlines:
        for i in range(len(sin)):
            sin[i] = sin[i][:-1]
    return sin

def loadheader(filename, dir="", silent=0, keepnewlines=0):
    infile = dirfile(filename, dir)
    if not silent:
        print("Loading ", infile, "...\n")
    fin = open(infile, 'r')
    line = '#'
    sin = []
    while line:
        line = fin.readline()
        if line[0] != '#':
            break
        else:
            sin.append(line)
    fin.close()
    if not keepnewlines:
        for i in range(len(sin)):
            sin[i] = sin[i][:-1]
    return sin

def fileempty(filename, dir="", silent=0, delifempty=0):
    """CHECK IF A FILE ACTUALLY HAS ANYTHING IN IT
    OR IF IT'S JUST CONTAINS BLANK / COMMENTED LINES"""
    filename = dirfile(filename, dir)
    gotdata = 0
    if os.path.exists(filename):
        fin = open(filename, 'r')
        line = 'x'
        while line and not gotdata:
            line = fin.readline()
            if line:
                if line[0] != '#':
                    gotdata = 1
        if delifempty:
            if not gotdata:
                os.remove(filename)
        fin.close()
    return (gotdata == 0)

def delfileifempty(filename, dir="", silent=0):
    fileempty(filename, dir, silent, 1)

def assigndict(keys, values):
    n = len(keys)
    if n != len(values):
        print("keys & values DON'T HAVE SAME LENGTH IN coeio.assigndict!")
    else:
        d = {}
        for i in range(n):
            d[keys[i]] = values[i]
        return d

def loaddict1(filename, dir="", silent=0):
    lines = loadfile(filename, dir, silent)
    dict = {}
    for line in lines:
        if line[0] != '#':
            words = line.split()
            key = str2num(words[0])
            val = ''  # if nothing there
            if len(words) == 2:
                val = str2num(words[1])
            elif len(words) > 2:
                val = []
                for word in words[1:]:
                    val.append(str2num(word))
                
            dict[key] = val
    return dict


def loaddict(filename, dir="", silent=0):
    lines = loadfile(filename, dir, silent)
    dict = {}
    for line in lines:
        if line[0] != '#':
            words = line.split()
            key = str2num(words[0])
            val = ''  # if nothing there
            valstr = ' '.join(words[1:])
            valtuple = False
            if valstr[0] in '[(' and valstr[-1] in '])':  # LIST / TUPLE!
                valtuple = valstr[0] == '('
                valstr = valstr[1:-1].replace(',', '')
                words[1:] = valstr.split()
            if len(words) == 2:
                val = str2num(words[1])
            elif len(words) > 2:
                val = []
                for word in words[1:]:
                    val.append(str2num(word))
                if valtuple:
                    val = tuple(val)
                
            dict[key] = val
    return dict


# THE LONG AWAITED MIXED FORMAT LOADER!
def loadcols(infile, format='', pl=0):
    """LOADS A DATA FILE CONTAINING COLUMNS OF DIFFERENT TYPES (STRING, FLOAT, & INT
    RETURNS A LIST OF LISTS
    format (OPTIONAL) INPUT AS A STRING, ONE LETTER (s, d, or f) FOR EACH COLUMN
    ARRAY OUTPUT FOR NUMBERS: ADD AN 'A' TO THE BEGINNING OF format 
    USAGE: labels, x, y = loadcols('~/A1689/arcsnew_lab.txt', format='sdd')"""
    txt = loadfile(infile)
##     line = txt[0]
##     words = string.split(line)
    while txt[0][0] == '#':
        txt = txt[1:]
    line = txt[0]
    words = line.split()
    ncols = len(words)
    data = [[]]
    for icol in range(ncols-1):
        data.append([])

    arrayout = 0

    if format:
        if format[0] == 'A':
            format = format[1:]
            arrayout = 1

    if not format:  # FIGURE IT OUT BASED ON FIRST LINE ONLY
        for word in words:
            try:
                datum = int(word)
                format += 'd'
            except:
                try:
                    datum = float(word)
                    format += 'f'
                except:
                    format += 's'

    #print format
    roundcols = []
    for line in txt:
        if line:
            if line[0] != '#':
                words = line.split()
                if pl:
                    print(line)
                for iword in range(len(words)):
                    if iword > len(format)-1:
                        print('EXTRA CONTENT IN LINE: ', end=' ')
                        print(words[iword:].join())
                        break
                    #print iword
                    word = words[iword]
                    formatum = format[iword]
                    if formatum == 'f':
                        datum = float(word)
                    elif formatum == 'd':
                        try:
                            datum = int(word)
                        except:
                            #datum = int(round(string.atof(word)))
                            datum = float(word)
                            try:
                                datum = roundint(datum)
                                if not (iword+1) in roundcols:
                                    roundcols.append(iword+1)
                            except:
                                pass
                    else:
                        datum = word
                    data[iword].append(datum)
    
    if roundcols:
        if len(roundcols) > 1:
            print('WARNING, THE FOLLOWING COLUMNS WERE ROUNDED FROM FLOAT TO INT: ', roundcols)
        else:
            print('WARNING, THE FOLLOWING COLUMN WAS ROUNDED FROM FLOAT TO INT: ', roundcols)

    if arrayout:
        for icol in range(ncols):
            if format[icol] in 'df':
                data[icol] = array(data[icol])

    return data

# CRUDE
def savecols(data, filename, format=''):
    ncols = len(data)
    nrows = len(data[0])
    if not format:
        for icol in range(ncols):
            datum = data[icol][0]
            if type(datum) == int:
                format += 'd'
            elif type(datum) == float:
                format += 'f'
            else:
                format += 's'

    # CHANGE format from 'sdd' TO ' %s %d %d\n'
    ff = ' '
    for f in format:
        if f == 'f':
            ff += '%.3f '
        else:
            ff += '%' + f + '  '
    format = ff[:-1]
    format += '\n'
    fout = open(filename, 'w')
    for irow in range(nrows):
        dataline = []
        for icol in range(ncols):
            dataline.append(data[icol][irow])
        fout.write(format % tuple(dataline))
    
    fout.close()


def savedata(data, filename, dir="", header="", separator="  ", format='', labels='', descriptions='', units='', notes=[], pf=0, maxy=300, machine=0, silent=0):
    """Saves an array as an ascii data file into an array."""
    # AUTO FORMATTING (IF YOU WANT, ALSO OUTPUTS FORMAT SO YOU CAN USE IT NEXT TIME w/o HAVING TO CALCULATE IT)
    # maxy: ONLY CHECK THIS MANY ROWS FOR FORMATTING
    # LABELS MAY BE PLACED ABOVE EACH COLUMN
    # IMPROVED SPACING

    dow = filename[-1] == '-'  # DON'T OVERWRITE
    if dow:
        filename = filename[:-1]

    tr = filename[-1] == '+'
    if tr:
        data = transpose(data)
        filename = filename[:-1]

    if machine:
        filename = 'datafile%d.txt' % machine # doubles as table number
    outfile = dirfile(filename, dir)

    if dow and os.path.exists(outfile):
        print(outfile, " ALREADY EXISTS")
    else:
        skycat = strend(filename, '.scat')
        if skycat:
            separator = '\t'
        if len(data.shape) == 1:
            data = reshape(data, (len(data), 1))
            #data = data[:,NewAxis]
        [ny,nx] = data.shape
        colneg = [0] * nx  # WHETHER THE COLUMN HAS ANY NEGATIVE NUMBERS: 1=YES, 0=NO
        collens = []
        if format:
            if type(format) == dict:  # CONVERT DICTIONARY FORM TO LIST
                dd = ' '
                for label in labels:
                    if label in list(format.keys()):
                        dd += format[label]
                    else:
                        print("WARNING: YOU DIDN'T SUPPLY A FORMAT FOR", label + ".  USING %.3f AS DEFAULT")
                        dd += '%.3f'
                    dd += '  '
                dd = dd[:-2] + '\n'  # REMOVE LAST WHITESPACE, ADD NEWLINE
                format = dd
                #print format
        else:
            if not silent:
                print("Formatting... ")
            coldec = [0] * nx  # OF DECIMAL PLACES
            colint = [0] * nx  # LENGTH BEFORE DECIMAL PLACE (INCLUDING AN EXTRA ONE IF IT'S NEGATIVE)
            #colneg = [0] * nx  # WHETHER THE COLUMN HAS ANY NEGATIVE NUMBERS: 1=YES, 0=NO
            colexp = [0] * nx  # WHETHER THE COLUMN HAS ANY REALLY BIG NUMBERS THAT NEED exp FORMAT : 1=YES, 0=NO

            if machine:
                maxy = 0
            if (ny <= maxy) or not maxy:
                yyy = list(range(ny))
            else:
                yyy = arange(maxy) * ((ny - 1.) / (maxy - 1.))
                yyy = yyy.astype(int)
            for iy in yyy:
                for ix in range(nx):
                    datum = data[iy,ix]
                    if isNaN(datum):
                        ni, nd = 1, 1
                    else:
                        if (abs(datum) > 1.e9) or (0 < abs(datum) < 1.e-5): # IF TOO BIG OR TOO SMALL, NEED exp FORMAT
                            ni, nd = 1, 3
                            colexp[ix] = 1
                        else:
                            ni = len("% d" % datum) - 1
                            if ni <= 3:
                                nd = ndec(datum, max=4)
                            else:
                                nd = ndec(datum, max=7-ni)
                            # Float32: ABOUT 7 DIGITS ARE ACCURATE (?)

                    if ni > colint[ix]:  # IF BIGGEST, YOU GET TO DECIDE NEG SPACE OR NO
                        colneg[ix] = (datum < 0)
                        #print '>', ix, colneg[ix], nd, coldec[ix]
                    elif ni == colint[ix]:  # IF MATCH BIGGEST, YOU CAN SET NEG SPACE ON (NOT OFF)
                        colneg[ix] = (datum < 0) or colneg[ix]
                        #print '=', ix, colneg[ix], nd, coldec[ix]
                    coldec[ix] = max([ nd, coldec[ix] ])
                    colint[ix] = max([ ni, colint[ix] ])

            #print colneg
            #print colint
            #print coldec

            collens = []
            for ix in range(nx):
                if colexp[ix]:
                    collen = 9 + colneg[ix]
                else:
                    collen = colint[ix] + coldec[ix] + (coldec[ix] > 0) + (colneg[ix] > 0)  # EXTRA ONES FOR DECIMAL POINT / - SIGN
                if labels and not machine:
                    collen = max((collen, len(labels[ix])))  # MAKE COLUMN BIG ENOUGH TO ACCOMODATE LABEL
                collens.append(collen)

            format = ' '
            for ix in range(nx):
                collen = collens[ix]
                format += '%'
                if colneg[ix]:  # NEGATIVE
                    format += ' '
                if colexp[ix]:  # REALLY BIG (EXP FORMAT)
                    format += '.3e'
                else:
                    if coldec[ix]:  # FLOAT
                        format += "%d.%df" % (collen, coldec[ix])
                    else:  # DECIMAL
                        format += "%dd" % collen
                if ix < nx - 1:
                    format += separator
                else:
                    format += "\n"
            if pf:
                print("format='%s\\n'" % format[:-1])
                
        # NEED TO BE ABLE TO ALTER INPUT FORMAT
        if machine:  # machine readable
            collens = [] # REDO collens (IN CASE format WAS INPUT)
            mformat = ''
            separator = ' '
            colformats = format.split('%')[1:]
            format = ''  # redoing format, too
            for ix in range(nx):
                #print ix, colformats
                cf = colformats[ix]
                format += '%'
                if cf[0] == ' ':
                    format += ' '
                cf = cf.strip()
                format += cf
                mformat += {'d':'I', 'f':'F', 'e':'E'}[cf[-1]]
                mformat += cf[:-1]
                if ix < nx - 1:
                    format += separator
                    mformat += separator
                else:
                    format += "\n"
                    mformat += "\n"
                # REDO collens (IN CASE format WAS INPUT)
                colneg[ix] = cf.find(' ') == -1
                if cf.find('e') > -1:
                    collen = 9 + colneg[ix]
                else:
                    cf = cf.split('.')[0]  # FLOAT: Number before '.'
                    cf = cf.split('d')[0]  # INT:   Number before 'd'
                    collen = int(cf)
                collens.append(collen)
        else:
            if not collens:
                collens = [] # REDO collens (IN CASE format WAS INPUT)
                colformats = format.split('%')[1:]
                for ix in range(nx):
                    cf = colformats[ix]
                    colneg[ix] = cf.find(' ') == -1
                    if cf.find('e') > -1:
                        collen = 9 + colneg[ix]
                    else:
                        cf = cf.split('.')[0]  # FLOAT: Number before '.'
                        cf = cf.split('d')[0]  # INT:   Number before 'd'
                        collen = int(cf)
                    if labels:
                        collen = max((collen, len(labels[ix])))  # MAKE COLUMN BIG ENOUGH TO ACCOMODATE LABEL
                    collens.append(collen)
            

##             if machine:  # machine readable
##                 mformat = ''
##                 for ix in range(nx):
##                     collen = collens[ix]
##                     if colexp[ix]: # EXP
##                         mformat += 'E5.3'
##                     elif coldec[ix]: # FLOAT
##                         mformat += 'F%d.%d' % (collen, coldec[ix])
##                     else: # DECIMAL
##                         mformat += 'I%d' % collen
##                     if ix < nx - 1:
##                         mformat += separator
##                     else:
##                         mformat += "\n"
                
        if descriptions:
            if type(descriptions) == dict:  # CONVERT DICTIONARY FORM TO LIST
                dd = []
                for label in labels:
                    dd.append(descriptions.get(label, ''))
                descriptions = dd

        if units:
            if type(units) == dict:  # CONVERT DICTIONARY FORM TO LIST
                dd = []
                for label in labels:
                    dd.append(units.get(label, ''))
                units = dd

        if not machine:
##             if not descriptions:
##                 descriptions = labels
            if labels:
                headline = ''
                maxcollen = 1
                for label in labels:
                    maxcollen = max([maxcollen, len(label)])
                for ix in range(nx):
                    label = labels[ix].ljust(maxcollen)
                    headline += '# %2d %s' % (ix+1, label)
                    if descriptions:
                        if descriptions[ix]:
                            headline += '  %s' % descriptions[ix]
                    headline += '\n'
                    ##headline += '# %2d %s  %s\n' % (ix+1, label, descriptions[ix])
                    #ff = '# %%2d %%%ds  %%s\n' % maxcollen  # '# %2d %10s %s\n' 
                    #headline += ff % (ix+1, labels[ix], descriptions[ix])
                    #headline += '# %2d %s\n' % (ix+1, descriptions[ix])
                headline += '#\n'
                headline += '#'
                colformats = format.split('%')[1:]
                if not silent:
                    print()
                for ix in range(nx):
                    cf = colformats[ix]
                    collen = collens[ix]
                    #label = labels[ix][:collen]  # TRUNCATE LABEL TO FIT COLUMN DATA
                    label = labels[ix]
                    label = label.center(collen)
                    headline += label + separator
                headline += '\n'
                if not header:
                    header = [headline]
                else:
                    if header[-1] != '.':  # SPECIAL CODE TO REFRAIN FROM ADDING TO HEADER
                        header.append(headline)

                if skycat:
                    headline1 = ''
                    headline2 = ''
                    for label in labels:
                        headline1 += label + '\t'
                        headline2 += '-' * len(label) + '\t'
                    headline1 = headline1[:-1] + '\n'
                    headline2 = headline2[:-1] + '\n'
                    header.append(headline1)
                    header.append(headline2)

        elif machine:  # Machine readable Table!
            maxlabellen = 0
            for ix in range(nx):
                flaglabel = 0
                if len(labels[ix]) >= 2:
                    flaglabel = (labels[ix][1] == '_')
                if not flaglabel:
                    labels[ix] = '  ' + labels[ix]
                if len(labels[ix]) > maxlabellen:
                    maxlabellen = len(labels[ix])
            # labelformat = '%%%ds' % maxlabellen
            if not header:
                header = []
                header.append('Title:\n')
                header.append('Authors:\n')
                header.append('Table:\n')
            header.append('='*80+'\n')
            header.append('Byte-by-byte Description of file: %s\n' % filename)
            header.append('-'*80+'\n')
            #header.append('   Bytes Format Units   Label    Explanations\n')
            headline = '   Bytes Format Units   '
            headline += 'Label'.ljust(maxlabellen-2)
            headline += '  Explanations\n'
            header.append(headline)
            header.append('-'*80+'\n')
            colformats = mformat.split()
            byte = 1
            for ix in range(nx):
                collen = collens[ix]
                headline = ' %3d-%3d' % (byte, byte + collen - 1) # bytes
                headline += ' '
                # format:
                cf = colformats[ix]
                headline += cf
                headline += ' ' * (7 - len(cf))
                # units:
                cu = ''
                if units:
                    cu = units[ix]
                if not cu:
                    cu = '---'
                headline += cu
                headline += '   '
                # label:
                label = labels[ix]
                headline += labels[ix].ljust(maxlabellen)
                # descriptions:
                if descriptions:
                    headline += '  '
                    headline += descriptions[ix]
                headline += '\n'
                header.append(headline)
                byte += collen + 1
            header.append('-'*80+'\n')
            if notes:
                for inote in range(len(notes)):
                    headline = 'Note (%d): ' % (inote+1)
                    note = notes[inote].split('\n')
                    headline += note[0]
                    if headline[-1] != '\n':
                        headline += '\n'
                    header.append(headline)
                    if len(note) > 1:
                        for iline in range(1, len(note)):
                            if note[iline]: # make sure it's not blank (e.g., after \n)
                                headline = ' ' * 10
                                headline += note[iline]
                                if headline[-1] != '\n':
                                    headline += '\n'
                                header.append(headline)
            header.append('-'*80+'\n')

        if not silent:
            print("Saving ", outfile, "...\n")

        fout = open(outfile, 'w')

        # SPECIAL CODE TO REFRAIN FROM ADDING TO HEADER:
        # LAST ELEMENT IS A PERIOD
        if header:
            if header[-1] == '.':
                header = header[:-1]
        
        for headline in header:
            fout.write(headline)
            if not (headline[-1] == '\n'):
                fout.write('\n')

        for iy in range(ny):
            fout.write(format % tuple(data[iy].tolist()))

        fout.close()


def loaddata(filename, dir="", silent=0, headlines=0):
    """Loads an ascii data file (OR TEXT BLOCK) into an array.
    Skips header (lines that begin with #), but saves it in the variable 'header', which can be accessed by:
    from coeio import header"""

    global header

    tr = 0
    if filename[-1] == '+':  # TRANSPOSE OUTPUT
        tr = 1
        filename = filename[:-1]

    if len(filename[0]) > 1:
        sin = filename
    else:
        sin = loadfile(filename, dir, silent)

    header = sin[0:headlines]
    sin = sin[headlines:]

    headlines = 0
    while (headlines < len(sin)) and (sin[headlines][0] == '#'):
        headlines = headlines + 1
    header[len(header):] = sin[0:headlines]

    ny = len(sin) - headlines
    if ny == 0:
        if headlines:
            ss = sin[headlines-1].split()[1:]
        else:
            ss = []
    else:
        ss = sin[headlines].split()

    nx = len(ss)
    #size = [nx,ny]
    data = FltArr(ny,nx)

    sin = sin[headlines:ny+headlines]

    for iy in range(ny):
        ss = sin[iy].split()
        for ix in range(nx):
            try:
                data[iy,ix] = float(ss[ix])
            except:
                print(ss)
                print(ss[ix])
                data[iy,ix] = float(ss[ix])
    
    if tr:
        data = transpose(data)
    
    if data.shape[0] == 1:  # ONE ROW
        return ravel(data)
    else:
        return data


def loadlist(filename, dir="./"):
    """Loads an ascii data file into a list.
    The file has one number on each line.
    Skips header (lines that begin with #), but saves it in the variable 'header'."""

    global header

    #    os.chdir("/home/coe/imcat/ksb/A1689txitxo/R/02/")

    infile = dirfile(filename, dir)
    print("Loading ", infile, "...\n")

    fin = open(infile, 'r')
    sin = fin.readlines()
    fin.close

    headlines = 0
    while sin[headlines][0] == '#':
        headlines = headlines + 1
    header = sin[0:headlines-1]

    n = len(sin) - headlines

    sin = sin[headlines:n+headlines]

    list = []
    for i in range(n):
        list.append(float(sin[i]))

    return list


def machinereadable(filename, dir=''):
    if filename[-1] == '+':
        filename = filename[:-1]
    filename = dirfile(filename, dir)
    fin = open(filename, 'r')
    line = fin.readline()
    return line[0] == 'T'  # BEGINS WITH Title:
    

def loadmachine(filename, dir="", silent=0):
    """Loads machine-readable ascii data file into a VarsClass()
    FORMAT...
    Title:
    Authors:
    Table:
    ================================================================================
    Byte-by-byte Description of file: datafile1.txt
    --------------------------------------------------------------------------------
       Bytes Format Units   Label    Explanations 
    --------------------------------------------------------------------------------
    (columns & descriptions)
    --------------------------------------------------------------------------------
    (notes)
    --------------------------------------------------------------------------------
    (data)
    """

    cat = VarsClass('')
    filename = dirfile(filename, dir)
    fin = open(filename, 'r')

    # SKIP HEADER
    line = ' '
    while line.find('Bytes') == -1:
        line = fin.readline()
    fin.readline()

    # COLUMNS & DESCRIPTIONS
    cols = []
    cat.labels = []
    line = fin.readline()
    while line[0] != '-':
        xx = []
        xx.append(int(line[1:4]))
        xx.append(int(line[5:8]))
        cols.append(xx)
        cat.labels.append(line[9:].split()[2])
        line = fin.readline()
    
    nx = len(cat.labels)

    # NOW SKIP NOTES:
    line = fin.readline()
    while line[0] != '-':
        line = fin.readline()

    # INITIALIZE DATA
    for ix in range(nx):
        exec('cat.%s = []' % cat.labels[ix])
    
    # LOAD DATA
    while line:
        line = fin.readline()
        if line:
            for ix in range(nx):
                s = line[cols[ix][0]-1:cols[ix][1]]
                #print cols[ix][0], cols[ix][1], s
                val = float(s)
                exec('cat.%s.append(val)' % cat.labels[ix])
    
    # FINALIZE DATA
    for ix in range(nx):
        exec('cat.%s = array(cat.%s)' % (cat.labels[ix], cat.labels[ix]))
    
    return cat


def loadpymc(filename, dir="", silent=0):
    filename = dirfile(filename, dir)
    ind = loaddict(filename+'.ind')
    i, data = loaddata(filename+'.out+')
    
    cat = VarsClass()
    for label in list(ind.keys()):
        ilo, ihi = ind[label]
        chunk = data[ilo-1:ihi]
        cat.add(label, chunk)
    
    return cat

class Cat2D:
    def __init__(self, filename='', dir="", silent=0, labels='x y z'.split()):
        if len(labels) == 2:
            labels.append('z')
        self.labels = labels
        if filename:
            if filename[-1] != '+':
                filename += '+'
            self.data = loaddata(filename, dir)
            self.assigndata()
    def assigndata(self):
        exec('self.%s = self.x = self.data[1:,0]' % self.labels[0])
        exec('self.%s = self.y = self.data[0,1:]' % self.labels[1])
        exec('self.%s = self.z = self.data[1:,1:]' % self.labels[2])
    def get(self, x, y, dointerp=0):
        ix = interp(x, self.x, arange(len(self.x)))
        iy = interp(y, self.y, arange(len(self.y)))
        if not dointerp:  # JUST GET NEAREST
            #ix = searchsorted(self.x, x)
            #iy = searchsorted(self.y, y)
            ix = roundint(ix)
            iy = roundint(iy)
            z = self.z[ix,iy]
        else:
            z = bilin2(iy, ix, self.z)
        return z


def loadcat2d(filename, dir="", silent=0, labels='x y z'):
    """INPUT: ARRAY w/ SORTED NUMERIC HEADERS (1ST COLUMN & 1ST ROW)
    OUTPUT: A CLASS WITH RECORDS"""
    if type(labels) == str:
        labels = labels.split()
    outclass = Cat2D(filename, dir, silent, labels)
    return outclass

def savecat2d(data, x, y, filename, dir="", silent=0):
    """OUTPUT: FILE WITH data IN BODY AND x & y ALONG LEFT AND TOP"""
    #y = y[NewAxis, :]
    y = reshape(y, (1, len(y)))
    data = concatenate([y, data])
    x = concatenate([[0], x])
    #x = x[:, NewAxis]
    x = reshape(x, (len(x), 1))
    data = concatenate([x, data], 1)
    if filename[-1] != '+':
        filename += '+'
    #savedata(data, filename)
    savedata1(data, filename)

def savedata1d(data, filename, dir="./", format='%6.5e ', header=""):
    fout = open(filename, 'w')
    for datum in data:
        fout.write('%d\n' % datum)
    fout.close()
                                                                               
def loadvars(filename, dir="", silent=0):
    """INPUT: CATALOG w/ LABELED COLUMNS
    OUTPUT: A STRING WHICH WHEN EXECUTED WILL LOAD DATA INTO VARIABLES WITH NAMES THE SAME AS COLUMNS
    >>> exec(loadvars('file.cat'))
    NOTE: DATA IS ALSO SAVED IN ARRAY data"""
    global data, labels, labelstr
    if filename[-1] != '+':
        filename += '+'
    data = loaddata(filename, dir, silent)
    labels = header[-1][1:].split()
    labelstr = labels.join(',')
    print(labelstr + ' = data')
    return 'from coeio import data,labels,labelstr\n' + labelstr + ' = data'  # STRING TO BE EXECUTED AFTER EXIT
    #return 'from coeio import data\n' + labelstr + ' = data'  # STRING TO BE EXECUTED AFTER EXIT


class VarsClass:
    def __init__(self, filename='', dir="", silent=0, labels='', labelheader='', headlines=0, loadheader=0):
        self.header = ''
        if filename:
            if strend(filename, '.fits'):  # FITS TABLE
                self.name = filename
                filename = dirfile(filename, dir)
                hdulist = pyfits.open(filename)
                self.labels = hdulist[1].columns.names
                tbdata = hdulist[1].data
                self.labels = labels or self.labels
                for label in self.labels:
                    print(label)
                    print(tbdata.field(label)[:5])
                    exec("self.%s = array(tbdata.field('%s'), 'f')" % (label, label))
                    print(self.get(label)[:5])
                self.updatedata()
            elif machinereadable(filename, dir):
                #self = loadmachine(filename, dir, silent)
                self2 = loadmachine(filename, dir, silent)
                self.labels = self2.labels[:]
                for label in self.labels:
                    exec('self.%s = self2.%s[:]' % (label, label))
                self.name = filename
                #self.header = '' # for now...
            else:
                if filename[-1] != '+':
                    filename += '+'
                self.name = filename[:-1]
                self.data = loaddata(filename, dir, silent, headlines)
                # NOTE header IS global, GETS READ IN BY loaddata
                if loadheader:
                    self.header = header or labelheader
                if header:
                    labelheader = labelheader or header[-1][1:]
                self.labels = labels or labelheader.split()
                # self.labels = string.split(header[-1][1:])
                # self.labelstr = string.join(self.labels, ', ')[:-2]
                self.assigndata()
        else:
            self.name = ''
            self.labels = []
        self.descriptions = {}
        self.units = {}
        self.notes = []
    def assigndata(self):
        for iii in range(len(self.labels)):
            label = self.labels[iii]
            try:
                exec('self.%s = self.data[iii]' % label)
            except:
                print('BAD LABEL NAME:', label)
                xxx[9] = 3
    def copy(self):
        #return copy.deepcopy(self)
        selfcopy = VarsClass()
        selfcopy.labels = self.labels[:]
        selfcopy.data = self.updateddata()
        selfcopy.assigndata()
        selfcopy.descriptions = self.descriptions
        selfcopy.units = self.units
        selfcopy.notes = self.notes
        selfcopy.header = self.header
        return selfcopy
    def updateddata(self):
        selflabelstr = ''
        for label in self.labels:
            #if label <> 'flags':
            selflabelstr += 'self.' + label + ', '
            #print label
            #exec('print self.%s') % label
            #exec('print type(self.%s)') % label
            #exec('print self.%s.shape') % label
            #print
        selflabelstr = selflabelstr[:-2]
        #print 'x', selflabelstr, 'x'
        #data1 = array([self.id, self.area])
        #print array([self.id, self.area])
        #s = 'data3 = array([%s])' % selflabelstr
        #print s
        #exec(s)
        #print data1
        #print 'data3 = array([%s])' % selflabelstr
        exec('data3 = array([%s])' % selflabelstr)
        return data3
    def updatedata(self):
        self.data = self.updateddata()
    def len(self):
        if self.labels:
            x = self.get(self.labels[0])
            l = len(x)
            try:
                l = l[:-1]
            except:
                pass
        else:
            l = 0
        return l
    def subset(self, good):
        #selfcopy = self.copy()
##         if len(self.id) <> len(good):
##             print "VarsClass: SUBSET CANNOT BE CREATED: good LENGTH = %d, data LENGTH = %d" % (len(self.id), len(good))
        if self.len() != len(good):
            print("VarsClass: SUBSET CANNOT BE CREATED: good LENGTH = %d, data LENGTH = %d" % (self.len(), len(good)))
        else:
            selfcopy = self.copy()
            data = self.updateddata()
            selfcopy.data = compress(good, data)
            selfcopy.assigndata()
            selfcopy.data = compress(good, self.data) # PRESERVE UN-UPDATED DATA ARRAY
            selfcopy.taken = compress(good, arange(self.len()))
            selfcopy.good = good
            return selfcopy
    def between(self, lo, labels, hi):
        """labels = list of labels or just one label"""
        if type(labels) == list:
            exec('good = between(lo, self.%s, hi)' % labels[0])
            for label in labels[1:]:
                exec('good = good * between(lo, self.%s, hi)' % label)
        else:
            exec('good = between(lo, self.%s, hi)' % labels)
        self.good = good
        return self.subset(good)
    def take(self, indices):
        indices = indices.astype(int)
        sub = VarsClass()
        sub.labels = self.labels[:]
        sub.taken = sub.takeind = indices
        sub.data = take(self.updateddata(), indices, 1)
        sh = sub.data.shape
        if len(sh) == 3:
            sub.data = reshape(sub.data, sh[:2])
        sub.assigndata()
        return sub
    def put(self, label, indices, values):
        #exec('put(self.%s, indices, values)' % label)
        exec('x = self.%s.copy()' % label)
        put(x, indices, values)
        exec('self.%s = x' % label)
    def takeid(self, id, idlabel='id'):
        selfid = self.get(idlabel).astype(int) # [6 4 5]
        i = argmin(abs(selfid - id))
        if selfid[i] != id:
            print("PROBLEM! ID %d NOT FOUND IN takeid" % id)
            return None
        else:
            return self.take(array([i]))
    def putid(self, label, id, value, idlabel='id'):
        #print "putid UNTESTED!!"  -- STILL TRUE
        #print "(Hit Enter to continue)"
        #pause()
        selfid = self.get(idlabel).astype(int) # [6 4 5]
        i = argmin(abs(selfid - id))
        if selfid[i] != id:
            print("PROBLEM! ID %d NOT FOUND IN putid" % id)
            return None
        else:
            exec('x = self.%s.copy()' % label)
            put(x, i, value)
            exec('self.%s = x' % label)
        #print self.takeid(id).get(label)
    def takeids(self, ids, idlabel='id'):
        #selfid = self.id.astype(int) # [6 4 5]
        selfid = self.get(idlabel).astype(int) # [6 4 5]
        indexlist = zeros(max(selfid)+1, int) - 1
        put(indexlist, selfid, arange(len(selfid)))  # [- - - - 1 2 0]
        #self.good = greater(selfid, -1)  # TOTALLY WRONG!  USED IN bpzphist
        indices = take(indexlist, array(ids).astype(int))  # ids = [4 6]  ->  indices = [1 0]
        #print type(indices[0])
        goodindices = compress(greater(indices, -1), indices)
        good = zeros(self.len(), int)
        #print 'takeids'
        good = good.astype(int)
        goodindices = goodindices.astype(int)
        #print type(good[0]) #good.type()
        #print type(goodindices[0])  #goodindices.type()
        #pause()
        put(good, goodindices, 1)
        self.good = good
        if -1 in indices:
            print("PROBLEM! NOT ALL IDS FOUND IN takeids!")
            print(compress(less(indices, 0), ids))
        return self.take(indices)
    def putids(self, label, ids, values, idlabel='id'):
        #print "putids UNTESTED!!"
        #print "putids not fully tested"
        #print "(Hit Enter to continue)"
        #pause()
        #selfid = self.id.astype(int) # [6 4 5]
        selfid = self.get(idlabel).astype(int) # [6 4 5]
        indexlist = zeros(max(selfid)+1, int) - 1
        put(indexlist, selfid, arange(len(selfid)))  # [- - - - 1 2 0]
        indices = take(indexlist, array(ids).astype(int))  # ids = [4 6]  ->  indices = [1 0]
        if -1 in indices:
            print("PROBLEM! NOT ALL IDS FOUND IN putids!")
            print(compress(less(indices, 0), ids))
        exec('x = self.%s.copy()' % label)
        if singlevalue(values):
            values = zeros(self.len(), float) + values
        put(x, indices, values)
        exec('self.%s = x' % label)
    def takecids(self, ids, idlabel='id'):  # only take common ids
        #selfid = self.id.astype(int) # [6 4 5]
        selfid = self.get(idlabel).astype(int) # [6 4 5]
        n = max((max(selfid), max(ids)))
        indexlist = zeros(n+1, int)
        #indexlist = zeros(max(selfid)+1)
        put(indexlist, selfid, arange(len(selfid))+1)  # [- - - - 1 2 0]
        indices = take(indexlist, array(ids).astype(int))  # ids = [4 6]  ->  indices = [1 0]
        indices = compress(indices, indices-1)
        goodindices = compress(greater(indices, -1), indices)
        good = zeros(self.len(), int)
        put(good, goodindices, 1)
        self.good = good
        return self.take(indices)
    def removeids(self, ids, idlabel='id'):
        selfid = self.get(idlabel).astype(int) # [6 4 5]
        if singlevalue(ids):
            ids = [ids]
        #newids = set(selfid) - set(ids)  # SCREWS UP ORDER!
        #newids = list(newids)
        newids = invertselection(ids, selfid)
        return self.takeids(newids)
    def get(self, label, orelse=None):
        if label in self.labels:
            exec('out = self.%s' % label)
        else:
            out = orelse
        return out
    def set(self, label, data):
        if singlevalue(data):
            data = zeros(self.len(), float) + data
        exec('self.%s = data' % label)
    def add(self, label, data):
        if 1:  #self.labels:
            if singlevalue(data):
                if self.len():
                    data = zeros(self.len(), float) + data
                else:
                    data = array([float(data)])
            elif self.len() and (len(data) != self.len()):
                print('WARNING!! in loadvarswithclass.add:')
                print('len(%s) = %d BUT len(id) = %d' % (label, len(data), self.len()))
                print()
        self.labels.append(label)
        exec('self.%s = data.astype(float)' % label)
    def assign(self, label, data):
        if label in self.labels:
            self.set(label, data)
        else:
            self.add(label, data)
    def append(self, self2, over=0):
        # APPENDS THE CATALOG self2 TO self
        labels = self.labels[:]
        labels.sort()
        labels2 = self2.labels[:]
        labels2.sort()
        if labels != labels2:
            print("ERROR in loadvarswithclass.append: labels don't match")
            xxx[9] = 3
        else:
            if over:  # OVERWRITE OLD OBJECTS WITH NEW WHERE IDS ARE THE SAME
                commonids = common(self.id, self2.id)
                if commonids:
                    selfuniqueids = invertselection(commonids, self.id)
                    self = self.takeids(selfuniqueids)
            for label in self.labels:
                exec('self.%s = concatenate((self.get(label), self2.get(label)))' % label)
            self.updatedata()
        return self
    def merge(self, self2, labels=None, replace=0):
        # self2 HAS NEW INFO (LABELS) TO ADD TO self
        if 'id' in self.labels:
            if 'id' in self2.labels:
                self2 = self2.takeids(self.id)
        labels = labels or self2.labels
        for label in labels:
            if label not in self.labels:
                self.add(label, self2.get(label))
            elif replace:
                exec('self.%s = self2.%s' % (label, label))
    def sort(self, label): # label could also be an array
        if type(label) == str:
            if (label == 'random') and ('random' not in self.labels):
                SI = argsort(random(self.len()))
            else:
                if label[0] == '-':  # Ex.: -odds
                    label = label[1:]
                    reverse = 1
                else:
                    reverse = 0
                exec('SI = argsort(self.%s)' % label)
                if reverse:
                    SI = SI[::-1]
        else:
            SI = argsort(label)  # label contains an array
        self.updatedata()
        self.data = take(self.data, SI, 1)
        self.assigndata()
    def findmatches(self, searchcat1, dtol=4):
        """Finds matches for self in searchcat1
        match distances within dtol, but may not always be closest
        see also findmatches2"""
        matchids = []
        dists = []
        searchcat = searchcat1.copy()
        #searchcat.sort('x')
        if 'dtol' in self.labels:
            dtol = self.dtol
        else:
            dtol = dtol * ones(self.len())
        for i in range(self.len()):
            if not (i % 100):
                print("%d / %d" % (i, self.len()))
            #matchid, dist = findmatch(searchcat.x, searchcat.y, self.x[i], self.y[i], dtol=dtol[i], silent=1, returndist=1, xsorted=1)  # silent=2*(i<>38)-1
            matchid, dist = findmatch(searchcat.x, searchcat.y, self.x[i], self.y[i], dtol=dtol[i], silent=1, returndist=1, xsorted=0)  # silent=2*(i<>38)-1
##             print self.x[i], self.y[i], matchid,
##             if matchid < self.len():
##                 print searchcat.id[matchid], searchcat.x[matchid], searchcat.y[matchid]
##             else:
##                 print
##             pause()
            matchids.append(matchid)
            dists.append(dist)
        matchids = array(matchids)
        dists = array(dists)
        matchids = where(equal(matchids, searchcat.len()), -1, matchids)
        self.assign('matchid', matchids)
        self.assign('dist', dists)
    def findmatches2(self, searchcat, dtol=0):
        """Finds closest matches for self within searchcat"""
        i, d = findmatches2(searchcat.x, searchcat.y, self.x, self.y)
        if dtol:
            i = where(less(d, dtol), i, -1)
        self.assign('matchi', i)
        self.assign('dist', d)
    def rename(self, oldlabel, newlabel):
        self.set(newlabel, self.get(oldlabel))
        i = self.labels.index(oldlabel)
        self.labels[i] = newlabel
        if self.descriptions:
            if oldlabel in list(self.descriptions.keys()):
                self.descriptions[newlabel] = self.descriptions[oldlabel]
        if self.units:
            if oldlabel in list(self.units.keys()):
                self.units[newlabel] = self.units[oldlabel]
    def save(self, name='', header='', format='', labels=1, pf=0, maxy=300, machine=0, silent=0):
        if type(labels) == list:
            self.labels = labels
        labels = labels and self.labels  # if labels then self.labels, else 0
        name = name or self.name  # if name then name, else self.name
        header = header or self.header  # if header then header, else self.header
        savedata(self.updateddata(), name+'+', labels=labels, header=header, format=format, pf=pf, maxy=maxy, machine=machine, descriptions=self.descriptions, units=self.units, notes=self.notes, silent=silent)
    def savefitstable(self, name='', header='', format={}, labels=1, overwrite=1):  # FITS TABLE
        name = name or recapfile(self.name, 'fits')  # if name then name, else self.name
        name = capfile(name, 'fits')  # IF WAS name (PASSED IN) NEED TO capfile
        if (not overwrite) and os.path.exists(name):
            print(name, 'ALREADY EXISTS, AND YOU TOLD ME NOT TO OVERWRITE IT')
        else:
            units = self.units
            header = header or self.header  # if header then header, else self.header
            if type(labels) == list:
                self.labels = labels
            labels = labels and self.labels  # if labels then self.labels, else 0
            collist = []
            for label in self.labels:
                a = self.get(label)
                if not pyfitsusesnumpy:
                    a = numarray.array(a)
                if label in list(units.keys()):
                    col = pyfits.Column(name=label, format=format.get(label, 'E'), unit=units[label], array=a)
                else:
                    col = pyfits.Column(name=label, format=format.get(label, 'E'), array=a)
                collist.append(col)
            cols = pyfits.ColDefs(collist)
            tbhdu = pyfits.new_table(cols)
            if not self.descriptions:
                delfile(name)
                tbhdu.writeto(name)
            else:
                hdu = pyfits.PrimaryHDU()
                hdulist = pyfits.HDUList(hdu)
                hdulist.append(tbhdu)
                prihdr = hdulist[1].header
                descriptions = self.descriptions
                for ilabel in range(len(labels)):
                    label = labels[ilabel]
                    if label in list(descriptions.keys()):
                        description = descriptions[label]
                        if len(description) < 48:
                            description1 = description
                            description2 = ''
                        else:
                            i = description[:45].rfind(' ')
                            description1 = description[:i]+'...'
                            description2 = '...'+description[i+1:]
                        prihdr.update('TTYPE%d'%(ilabel+1), label, description1)
                        if description2:
                            prihdr.update('TFORM%d'%(ilabel+1), format.get(label, 'E'), description2)
                for inote in range(len(self.notes)):
                    words = self.notes[inote].split('\n')
                    for iword in range(len(words)):
                        word = words[iword]
                        if word:
                            if iword == 0:
                                prihdr.add_comment('(%d) %s' % (inote+1, word))
                            else:
                                prihdr.add_comment('    %s' % word)
                                #prihdr.add_blank(word)
                headlines = header.split('\n')
                for headline in headlines:
                    if headline:
                        key, value = headline.split('\t')
                        prihdr.update(key, value)
                hdulist.writeto(name)
    def pr(self):  # print
        self.save('tmp.cat')
        os.system('cat tmp.cat | more')
        os.remove('tmp.cat')
##     def takecids(self, ids):
##         selfid = self.id.astype(int)
##         ids = ids.astype(int)
##         n = max((max(selfid), max(ids)))
##         takeme1 = zeros(n+1)
##         takeme2 = zeros(n+1)
##         put(takeme1, selfid, 1)
##         put(takeme2, ids, 1)
##         takeme = takeme1 * takeme2
##         takeme = take(takeme, selfid)
##         return self.subset(takeme)
##     def labelstr(self):
##         return string.join(labels, ', ')[:-2]
        # FORGET THIS!  JUST USE copy.deepcopy
##         selfcopy = VarsClass()
##         selfcopy.data = self.data[:]
##         selfcopy.labels = self.labels.copy()

def loadvarswithclass(filename, dir="", silent=0, labels='', header='', headlines=0):
    """INPUT: CATALOG w/ LABELED COLUMNS
    OUTPUT: A CLASS WITH RECORDS NAMED AFTER EACH COLUMN
    >>> mybpz = loadvars('my.bpz')
    >>> mybpz.id
    >>> mybpz.data -- ARRAY WITH ALL DATA"""
    outclass = VarsClass(filename, dir, silent, labels, header, headlines)
    #outclass.assigndata()
    return outclass

loadcat = loadvarswithclass

#def loadcat(filename, dir="", silent=0):
def loadimcat(filename, dir="", silent=0):
    """LOADS A CATALOG CREATED BY IMCAT
    STORES VARIABLES IN A DICTIONARY OF ARRAYS!"""
    
    infile = dirfile(filename, dir)
    if not silent:
        print("Loading ", infile, "...\n")

    fin = open(infile, 'r')
    sin = fin.readlines()
    fin.close

    headlines = 0
    while sin[headlines][0] == '#':
        headlines = headlines + 1
    names = sin[headlines-1][1:].split()

    sin = sin[headlines:]  # REMOVE HEADLINES
    nx = len(names)
    ny = len(sin)
    data = FltArr(ny,nx)
    
    for iy in range(ny):
        ss = sin[iy].split()
        for ix in range(nx):
            data[iy,ix] = float(ss[ix])

    cat = {}
    for i in range(nx):
        cat[names[i]] = data[:,i]

    return cat

def savedict(dict, filename, dir="", silent=0):
    """SAVES A DICTIONARY OF STRINGS"""
    outfile = dirfile(filename, dir)
    fout = open(outfile, 'w')
    for key in list(dict.keys()):
        fout.write('%s %s\n' % (key, dict[key]))
    fout.close()

def savefile(lines, filename, dir="", silent=0):
    """SAVES lines TO filename"""
    outfile = dirfile(filename, dir)
    fout = open(outfile, 'w')
    for line in lines:
        if line[-1] != '\n':
            line += '\n'
        fout.write(line)
    fout.close()


#def savecat(cat, filename, dir="./", silent=0):
def saveimcat(cat, filename, dir="./", silent=0):
    # DOESN'T WORK RIGHT YET!!!  HEADER INCOMPLETE.
    """SAVES A DICTIONARY OF 1-D ARRAYS AS AN IMCAT CATALOGUE"""
    outfile = dirfile(filename, dir)
    fout = open(outfile, 'w')
    fout.write("# IMCAT format catalogue file -- edit with 'lc' or my Python routines\n")
    
    # COLUMN HEADERS
    fout.write("#")
    for key in list(cat.keys()):
        fout.write(key.rjust(15))
    fout.write("\n")

    n = len(cat[key])
    for i in range(n):
        fout.write(" ")
        keys = list(cat.keys())
        keys.sort()
        for key in keys:
            x = cat[key][i]
            if (x - int(x)):
                fout.write("%15.5f" % x)
            else:
                fout.write("%15d" % x)      
        fout.write("\n")

    fout.close()

def prunecols(infile, cols, outfile, separator=" "):
    """TAKES CERTAIN COLUMNS FROM infile AND OUTPUTS THEM TO OUTFILE
    COLUMN NUMBERING STARTS AT 1!
    ALSO AVAILABLE AS STANDALONE PROGRAM prunecols.py"""
    fin = open(infile, 'r')
    sin = fin.readlines()
    fin.close()

    fout = open(outfile, 'w')
    for line in sin:
        print(line)
        line = line.strip()
        words = line.split(separator)
        print(words)
        for col in cols:
            fout.write(words[col-1] + separator)
        fout.write("\n")
    fout.close()


#################################
# SExtractor/SExSeg CATALOGS / CONFIGURATION FILES

class SExSegParamsClass:
    def __init__(self, filename='', dir="", silent=0, headlines=0):
        # CONFIGURATION
        #   configkeys -- PARAMETERS IN ORDER
        #   config[key] -- VALUE
        #   comments[key] -- COMMENTS (IF ANY)
        # PARAMETERS
        #   params -- PARAMETERS IN ORDER
        #   comments[key] -- COMMENTS (IF ANY)
        self.name = filename
        self.configkeys = []
        self.config = {}
        self.comments = {}
        self.params = []
        txt = loadfile(filename, dir, silent)
        for line in txt:
            if line.strip() and (line[:1] != '#'):
                # READ FIRST WORD AND DISCARD IT FROM line
                key = line.split()[0]
                line = line[len(key):]
                # READ COMMENT AND DISCARD IT FROM line
                i = line.find('#')
                if i > -1:
                    self.comments[key] = line[i:]
                    line = line[:i]
                else:
                    self.comments[key] = ''
                # IF ANYTHING IS LEFT, IT'S THE VALUE, AND YOU'VE BEEN READING FROM THE CONFIGURATION SECTION
                # OTHERWISE IT WAS A PARAMETER (TO BE INCLUDED IN THE SEXTRACTOR CATALOG)
                line = line.strip()
                if line.strip():  # CONFIGURATION
                    self.configkeys.append(key)
                    self.config[key] = line
                else:  # PARAMETERS
                    self.params.append(key)

    def save(self, name='', header=''):
        name = name or self.name  # if name then name, else self.name
        # QUICK CHECK: IF ANY CONFIG PARAMS WERE ADDED TO THE DICT, BUT NOT TO THE LIST:
        for key in list(self.config.keys()):
            if key not in self.configkeys:
                self.configkeys.append(key)
        # OKAY...
        fout = open(name, 'w')
        #fout.write('#######################################\n')
        #fout.write('# CONFIGURATION\n')
        #fout.write('\n')
        fout.write('# ----- CONFIGURATION -----\n')
        for key in self.configkeys:
            line = ''
            line += key.ljust(20) + ' '
            value = self.config[key]
            comment = self.comments[key]
            if not comment:
                line += value
            else:
                line += value.ljust(20) + ' '
                line += comment
            line += '\n'
            fout.write(line)
        fout.write('\n')
        #fout.write('#######################################\n')
        #fout.write('# PARAMETERS\n')
        #fout.write('\n')
        fout.write('# ----- PARAMETERS -----\n')
        for param in self.params:
            line = ''
            comment = self.comments[param]
            if not comment:
                line += param
            else:
                line += param.ljust(20) + ' '
                line += comment
            line += '\n'
            fout.write(line)
        fout.close()

    def merge(self, filename='', dir="", silent=0, headlines=0):
        self2 = loadsexsegconfig(filename, dir, silent, headlines)
        for key in self2.configkeys:
            self.config[key] = self2.config[key]
        if self2.params:
            self.params = self2.params
        for key in list(self2.comments.keys()):
            if self2.comments[key]:
                self.comments[key] = self2.comments[key]
        

def loadsexsegconfig(filename='', dir="", silent=0, headlines=0):
    return SExSegParamsClass(filename, dir, silent, headlines)


def loadsexcat(infile, purge=1, maxflags=8, minfwhm=1, minrf=0, maxmag=99, magname="MAG_AUTO", ma1name='APER', silent=0, dir=''):
    """>>> exec(loadsexcat('sexfile.cat'<, ...>))
    LOADS A SEXTRACTOR CATALOG DIRECTLY INTO VARIABLES x, y, fwhm, etc.
    PURGES (DEFAULT) ACCORDING TO flags, fwhm, mag (AS SET)
    NOW MODELED AFTER loadvars -- OUTPUT STRING NEEDS TO BE EXECUTED, THEN ALL VARIABLES ARE LOADED
    NOW TAKES ON *ALL* VARIABLES, AND ADJUSTS NAMES ACCORDINGLY"""
    # outdata is a list of arrays (most are 1-D, but some (mag_aper) are 2-D)
    #global data, labels, labelstr, params, paramstr, outdata
    global params, paramstr, data, fullparamnames
    #global id, x, y, fwhm, mag, magerr, magauto, magerrauto, magautoerr, flags, a, b, theta, stellarity, rf, ell, rk, assoc, magaper, magapererr, magerraper, ma1, ema1, mb, emb, cl, flag, xpeak, ypeak, area

    # Txitxo variable translation:
    # cl = stellarity
    # flag = flags
    # ma1 = mag_aper   ema1 = error for mag_aper
    # mb  = mag_auto   emb  = error for mag_auto

    infile = join(dir, infile)
    if not silent:
        print("LOADING SExtractor catalog " + infile, end=' ') 
    
    #req = {'fwhm': 1, 'mag': 99, 'flags': 4}  # REQUIREMENTS FOR DATA TO BE INCLUDED (NOT PURGED)
    req = {}
    req['FLAGS'] = maxflags
    req['FWHM'] = minfwhm
    req['MAG'] = maxmag
    req['RF'] = minrf

    if magname:
        magname = magname.upper()
        if magname[:4] != 'MAG_':
            magname = 'MAG_' + magname
        #magerrname = 'MAGERR_' + magname[-4:]
        magerrname = 'MAGERR_' + magname[4:]
    else:
        magerrname = ''
        ma1name = ''

    sin = loadfile(infile, silent=1)

    # REMOVE HEADLINES FROM sin, CREATE header
    header = []
    while sin[0][0] == "#":
        if sin[0][1] != '#':
            header.append(sin[0])  # Only add lines beginning with single #
        sin = sin[1:]

    nx = len(sin[0].split())
    ny = len(sin)
    data = FltArr(ny,nx)

    for iy in range(ny):
        ss = sin[iy].split()
        for ix in range(nx):
            try:
                data[iy,ix] = float(ss[ix])
            except:
                print(iy, ix, nx)
                print(ss)
                die()

    data = transpose(data)
    paramcol = {}
    params = []
    
    flags = None
    fwhm = None
    mag = None
    rf = None

    #print 'TRYING NEW MAG ASSIGNMENT...'
    lastcol = 0  # COLUMN OF PREVIOUS PARAM
    lastparam = ''
    params = []
    fullparamnames = []
    for headline in header:
        ss = headline.split()  # ['#', '15', 'X_IMAGE', 'Object position along x', '[pixel]']
        if len(ss) == 1:
            break
        col = int(ss[1])  # 15  -- DON'T SUBTRACT 1 FROM col!  DON'T WANT A 0 COLUMN!  FACILITATES DATA DISTRIBUTION
        ncols = col - lastcol
        param = ss[2]    # "X_IMAGE"
        fullparamnames.append(param)
        if param[-1] == ']':
            param = param.split('[')[0]
        if param[:4] == "MAG_":  # MAG_AUTO or MAG_APER but not MAGERR_AUTO
            #if (param == magname) or not magname or 'MAG' not in paramcol.keys():  # magname IF YOU ONLY WANT MAG_AUTO (DEFAULT)
            if (param == magname) or not magname:  # magname IF YOU ONLY WANT MAG_AUTO (DEFAULT)
                magname = param
                param = "MAG"
        if param[:7] == "MAGERR_":  # MAGERR_AUTO or MAGERR_APER
            #if (param == magerrname) or not magerrname or 'MAG' not in paramcol.keys():  # magname IF YOU ONLY WANT MAG_AUTO (DEFAULT)
            if (param == magerrname) or not magerrname:  # magname IF YOU ONLY WANT MAG_AUTO (DEFAULT)
                magerrname = param
                param = "MAGERR"
        if param[-6:] == "_IMAGE":  # TRUNCATE "_IMAGE"
            param = param[:-6]
        if param in ["FLAGS", "IMAFLAGS_ISO"]:
            if not flags:
                flags = ravel(data[col-1]).astype(int)
                param = "FLAGS"
            else:
                flags = bitwise_or(flags, ravel(data[col-1]).astype(int))  # "FLAGS" OR "IMAFLAGS_ISO"
                param = ''
                lastcol += 1
##      if (param == "FLAGS") and paramcol.has_key("FLAGS"):
##          param = "SHWAG"  # "IMAFLAGS_ISO" (THE REAL FLAGS) HAVE ALREADY BEEN FOUND
##      if param == "IMAFLAGS_ISO":  # FLAGS OR-ED WITH FLAG MAP IMAGE
##          param = "FLAGS"
        #paramcol[param] = col
        #if vector > 1
        # ASSIGN COLUMN(S), NOW THAT WE KNOW HOW MANY THERE ARE
        if param != lastparam and param:
            if lastcol:
                paramcol[lastparam] = arange(ncols) + lastcol
            lastcol = col
            lastparam = param
            params.append(param)
        #print params

##     # IN CASE WE ENDED ON AN ARRAY (MAG_APER[4]) -- TAKEN CARE OF BELOW?
##     if param == lastparam:
##         if lastcol:
##             paramcol[lastparam] = arange(ncols) + lastcol
##         lastcol = col
##         lastparam = param
##         params.append(param)
   
    #print len(params)


    bigparamnames = params[:]
    paramstr = params.join(',')
    # ASSIGN LAST COLUMN(S)
    ncols = nx - lastcol + 1
    paramcol[lastparam] = arange(ncols) + lastcol

    col = paramcol.get("FWHM")
    #fwhm = col and ravel(data[col-1])
    if col.any(): fwhm = ravel(data[col-1])
    col = paramcol.get("FLUX_RADIUS")
    #rf = col and ravel(data[col-1])
    if col.any(): rf = ravel(data[col-1])
    col = paramcol.get("MAG")
    #mag = col and ravel(data[col-1])
    if col.any(): mag = ravel(data[col-1])

    good = ones(ny)
    if not silent:
        print(sum(good), end=' ')
    if purge:
        if "FLAGS" in req and (flags != None):
            good = less(flags, maxflags)
        if "FWHM" in req and (fwhm != None):
            good = good * greater(fwhm, minfwhm)
        if "RF" in req and (rf != None):
            good = good * greater(rf, minrf)
        if "MAG" in req and (mag != None):
            good = good * less(mag, maxmag)
            
    if not silent:
        print(sum(good))

    if purge and not alltrue(good):
        data = compress(good, data)
        if (flags != None): flags = compress(good, flags)
        if (mag != None): mag = compress(good, mag)
        if (fwhm != None): fwhm = compress(good, fwhm)
        if (rf != None): rf = compress(good, rf)

    outdata = []
    #params = paramcol.keys()
    # RENAME params
    paramtranslate = {'NUMBER':'id', 'CLASS_STAR':'stellarity', 'KRON_RADIUS':'rk', 'FLUX_RADIUS':'rf', 'ISOAREA':'area'}
    for ii in range(len(params)):
        param = params[ii]
        param = paramtranslate.get(param, param)  # CHANGE IT IF IN DICTIONARY, OTHERWISE LEAVE IT ALONE
        param = param.replace('_IMAGE', '')
        param = param.replace('PROFILE', 'PROF')
        param = param.lower()
        param = param.replace('_', '')
        #param = string.replace(param, 'magerr', 'dmag')
        #if param in ['a', 'b']:
        #    param = string.upper(param)
        params[ii] = param
        
    #print params
    #for kk in bigparamnames: #paramcol.keys():
    for ii in range(len(bigparamnames)): #paramcol.keys():
        pp = params[ii]
        #print
        #print pp
        if pp in ['flags', 'fwhm', 'rf', 'mag']:
            #exec('print type(%s)' % pp)
            #exec('print '+pp)
            exec('outdata.append(%s)' % pp)
            #outdata.append(flags)
        else:
            kk = bigparamnames[ii]
            col = paramcol[kk]
            if len(col) == 1:
                #exec(kk + '= data[col-1]')
                #print data[col-1]
                #print shape(data[col-1])
                #print type(data[col-1])
                outdata.append(ravel(data[col-1]))
            else:
                #exec(kk + '= take(data, col-1)')
                outdata.append(ravel(take(data, col-1)))

    paramstr = params.join(',')
    #exec(paramstr + ' = outdata')

    # CALCULATE ell (IF NOT CALCULATED ALREADY)
    #print params
    #print params.index('a')
    #print len(outdata)
    if 'ell' not in params:
        if 'a' in params and 'b' in params:
            a = outdata[params.index('a')]
            b = outdata[params.index('b')]
            try:
                ell = 1 - b / a
            except:
                ell = a * 0.
                for ii in range(len(a)):
                    if a[ii]:
                        ell[ii] = 1 - b[ii] / a[ii]
                    else:
                        ell[ii] = 99
            params.append('ell')
            paramstr += ', ell'
            outdata.append(ell)
            fullparamnames.append('ELLIPTICITY')

    # COMBINE flags & imaflagsiso
    if 'imaflagsiso' in params:
        flags = outdata[params.index('flags')]
        imaflagsiso = outdata[params.index('imaflagsiso')]
        flags = bitwise_or(flags.astype(int), imaflagsiso.astype(int))  # "FLAGS" OR "IMAFLAGS_ISO"
        outdata[params.index('flags')] = flags.astype(float)

    
    # FOR Txitxo's photometry.py
    #print 'COMMENTED OUT photometry.py LINES...'
    photcom = '\n'
    if 'stellarity' in params:
        photcom += 'cl = stellarity\n'
    if 'flags' in params:
        photcom += 'flag = flags\n'
    if 'mag' in params:
        photcom += 'mb = mag\n'
    if 'magerr' in params:
        photcom += 'emb = ravel(magerr)\n'

    if ma1name:
        #magtype = string.lower(string.split(ma1name, '_')[-1])
        magtype = ma1name.split('[')[0].lower()
        magtype = {'profile':'prof', 'isophotal':'iso'}.get(magtype, magtype)
        pos = ma1name.find('[')
        if pos > -1:
            ma1i = int(ma1name[pos+1:-1])
        else:
            ma1i = 0
        if magtype == 'aper':
            photcom += 'ma1 = ravel(magaper[%d])\n' % ma1i
            photcom += 'ema1 = ravel(magerraper[%d])\n' % ma1i
        else:
            photcom += 'ma1 = ravel(mag%s)\n' % magtype
            photcom += 'ema1 = ravel(magerr%s)\n' % magtype
        
    data = outdata
    #print 'params:', params
    #print photcom
    #outstr = 'from coeio import params,paramstr,data,fullparamnames\n' + paramstr + ' = data' + photcom
    #print outstr
    #return outstr
    #return 'from coeio import data,labels,labelstr,params,outdata\n' + paramstr + ' = outdata' + photcom  # STRING TO BE EXECUTED AFTER EXIT
    return 'from coeio import params,paramstr,data,fullparamnames\n' + paramstr + ' = data' + photcom  # STRING TO BE EXECUTED AFTER EXIT

def loadsexcat2(infile, purge=1, maxflags=8, minfwhm=1, minrf=0, maxmag=99, magname="MAG_AUTO", ma1name='APER', silent=0, dir=''):
    """RETURNS A VarsClass() VERSION OF THE CATALOG"""
    loadsexcat(infile, purge=purge, maxflags=maxflags, minfwhm=minfwhm, minrf=minrf, maxmag=maxmag, magname=magname, ma1name=ma1name, silent=silent, dir=dir)
    # LOADS infile INTO data, params...
    cat = VarsClass()
    cat.name = infile
    cat.data = data
    cat.labels = params
    cat.assigndata()
    for label in cat.labels:
        exec('cat.%s = ravel(cat.%s).astype(float)' % (label, label))
    #cat.flags = cat.flags[NewAxis,:]
    #cat.flags = cat.flags.astype(float)
    return cat

def loadsexdict(sexfile):
    """LOADS A SEXTRACTOR CONFIGURATION (.sex) FILE INTO A DICTIONARY
       COMMENTS NOT LOADED"""
    sextext = loadfile(sexfile)
    sexdict = {}
    for line in sextext:
        if line:
            if line[0] != '#':
                words = line.split()
                if len(words) > 1:
                    key = words[0]
                    if key[0] != '#':
                        # sexdict[words[0]] = str2num(words[1])
                        restofline = words[1:].join()
                        value = restofline.split('#')[0]
                        if value[0] == '$':
                            i = value.find('/')
                            value = os.getenv(value[1:i]) + value[i:]
                        sexdict[key] = str2num(value)
    return sexdict

def savesexdict(sexdict, sexfile):
    """SAVES A SEXTRACTOR CONFIGURATION (.sex) FILE
    BASED ON THE sexdict DICTIONARY"""
    fout = open(sexfile, 'w')
    keys = list(sexdict.keys())
    keys.sort()
    for key in keys:
        fout.write('%s\t%s\n' % (key, sexdict[key]))
    fout.close()

#################################
# DS9 REGIONS FILES

def saveregions1(x, y, filename, coords='image', color="green", symbol="circle", size=0, width=0):
    """SAVES POSITIONS AS A ds9 REGIONS FILE"""
    fout = open(filename, 'w')
    fout.write('global color=' + color + ' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    #fout.write("image\n")
    fout.write(coords+"\n")
    n = len(x)
    for i in range(n):
        if not size and not width:
            fout.write("%s point %s %s\n" % (symbol, x[i], y[i]))
        else:
            sout = '%s %s %s' % (symbol, x[i], y[i])
            if size:
                sout += ' %d' % size
            if width:
                sout += ' # width = %d' % width
            sout += '\n'
            fout.write(sout)
##         if size:
##             fout.write("%s %6.1f %6.1f %d\n" % (symbol, x[i], y[i], size))
##         else:
##             fout.write("%s point %6.1f %6.1f\n" % (symbol, x[i], y[i]))
    
    fout.close()

def saveregions(x, y, filename, labels=[], precision=1, coords='image', color="green", symbol="circle", size=0, width=0):
    """SAVES POSITIONS AND LABELS AS A ds9 REGIONS FILE"""
    fout = open(filename, 'w')
    fout.write('global color=' + color + ' font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    #fout.write("image\n")
    fout.write(coords+"\n")
    n = len(x)
    for i in range(n):
        if not size and not width:
            sout = '%s point %s %s' % (symbol, x[i], y[i])
        else:
            sout = '%s %s %s' % (symbol, x[i], y[i])
            if size:
                sout += ' %d' % size
            if width:
                sout += ' # width = %d' % width
        if i < len(labels):
            label = "%%.%df" % precision % labels[i]
            sout += ' # text={%s}' % label
        print(sout)
        sout += '\n'
        fout.write(sout)
    
    fout.close()

def savelabels(x, y, labels, filename, coords='image', color="green", symbol="circle", precision=1, fontsize=12, bold=1):
    """SAVES POSITIONS AS A ds9 REGIONS FILE"""
    if type(labels) in [int, float]:
        labels = arange(len(x)) + 1
    fout = open(filename, 'w')
    fout.write('global color=%s font="helvetica %d %s" select=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n' % (color, fontsize, ['normal', 'bold'][bold]))
    fout.write(coords+"\n")
    n = len(x)
    for i in range(n):
        label = labels[i]
        #if type(label) == int:  # IntType
        #    label = "%d" % label
        #elif type(label) == float: # FloatType
        #    label = "%%.%df" % precision % label
        label = "%%.%df" % precision % label
        fout.write("text %d %d {%s}\n" % (x[i], y[i], label))
    fout.close()


#################################
# FITS FILES

def savefits(data, filename, dir="", silent=0):
    """SAVES data (A 2-D ARRAY) AS filename (A .fits FILE)"""
    # THIS PROGRAM HAS FEWER OPTIONS THAN writefits IN fitsio, SO IT GETS THE JOB DONE EASILY!
    filename = capfile(filename, '.fits')
    filename = dirfile(filename, dir)
    if not silent:
        print('savefits:', filename)
    #print type(data)
    if os.path.exists(filename):
        os.remove(filename)
    # UNLESS $NUMERIX IS SET TO numpy, pyfits(v1.1b) USES NumArray
    if not pyfitsusesnumpy:
        data = numarray.array(data.tolist())
    pyfits.writeto(filename, data)

def loadfits(filename, dir=""):
    """READS in the data of a .fits file (filename)"""
    filename = capfile(filename, '.fits')
    filename = dirfile(filename, dir)
    if os.path.exists(filename):
        # CAN'T RETURN data WHEN USING memmap
        # THE POINTER GETS MESSED UP OR SOMETHING
        #return pyfits.open(filename, memmap=1)[0].data
        data = pyfits.open(filename)[0].data
        # UNLESS $NUMERIX IS SET TO numpy, pyfits(v1.1b) USES NumArray
        if not pyfitsusesnumpy:
            data = array(data)  # .tolist() ??
        return data
    else:
        print()
        print(filename, "DOESN'T EXIST")
        FILE_DOESNT_EXIST[9] = 3

def fitssize(filename):
    """RETURNS (ny, nx)"""
    filename = capfile(filename, '.fits')
    return pyfits.open(filename, memmap=1)[0]._dimShape()

## def fits2int(filename):
##     """CONVERTS A FITS FILE TO INTEGER"""
##     filename = capfile(filename, '.fits')
##     if os.path.exists(filename):
##         fitsio.writefits(filename, fitsio.readfits(filename), 16)
##     else:
##         print filename, "DOESN'T EXIST"

## def fits2float(filename):
##     """CONVERTS A FITS FILE TO FLOAT"""
##     filename = capfile(filename, '.fits')
##     if os.path.exists(filename):
##         fitsio.writefits(filename, fitsio.readfits(filename), -32)
##     else:
##         print filename, "DOESN'T EXIST"

## def fits8to16(filename):
##     """CONVERTS A FITS int8 FILE TO int16"""
##     filename = capfile(filename, '.fits')
##     if os.path.exists(filename):
##         data = loadfits(filename)
##         #data = where(less(data, 0), 256+data, data)
##         data = data % 256
##         oldfile = filename[:-5] + '8.fits'
##         os.rename(filename, oldfile)
##         savefits(data, filename)
##     else:
##         print filename, "DOESN'T EXIST"

def txt2fits(textfile, fitsfile):
    """CONVERTS A TEXT FILE DATA ARRAY TO A FITS FILE"""
    savefits(loaddata(textfile), fitsfile)

## def fitsheadval(file, param):
##     return fitsio.parsehead(fitsio.gethead(file), param)

## def enlargefits(infits, outsize, outfits):
##     """ENLARGES A FITS FILE TO THE DESIRED SIZE, USING BI-LINEAR INTERPOLATION
##     outsize CAN EITHER BE A TUPLE/LIST OR A FILE TO GRAB THE SIZE FROM"""
##     # RUNS SLOW!!
##     # WOULD RUN A LITTLE QUICKER IF WE COULD FILL EACH BOX BEFORE MOVING ON TO THE NEXT ONE
##     # RIGHT NOW I'M DOING ROW BY ROW WHICH DOES EACH BOX SEVERAL TIMES...

##     if type(outsize) == StringType:  # GRAB SIZE FROM A FILE
##         print "DETERMINING SIZE OF OUTPUT fits FILE"
##         fits = fitsio.readfits(outsize)
##         [nxout, nyout] = fits['naxes'][0:2]
##     else:  # SIZE GIVEN
##         (nxout, nyout) = outsize

##     print "CREATING OUTPUT DATA ARRAY..."
##     dataout = FltArr(nyout, nxout) * 0
        
##     print "READING IN INPUT fits FILE"
##     fits = fitsio.readfits(infits)
##     [nxin, nyin] = fits['naxes'][0:2]
##     datain = fits['data']

##     print "Interpolating data to full-size fits file (size [%d, %d])... " % (nxout, nyout),
##     print "%4d, %4d" % (0, 0),

##     # TRANSLATION FROM in COORDS TO out COORDS
##     dx = 1. * (nxout - 1) / (nxin - 1)
##     dy = 1. * (nyout - 1) / (nyin - 1)

##     # PRINT CREATE BOXES (4 POINTS) IN in SPACE, TRANSLATED TO out COORDS
##     byout = array([0., dy]) 
##     # GO THROUGH out COORDS
##     iyin = 0
##     for iyout in range(nyout):
##         if iyout > byout[1]:
##             byout = byout + dy
##             iyin = iyin + 1
##         bxout = array([0., dx])
##         ixin = 0
##         for ixout in range(nxout):
##             print "\b" * 11 + "%4d, %4d" % (ixout, iyout),
##             if ixout > bxout[1]:
##                 bxout = bxout + dx
##                 ixin = ixin + 1
##             # MAYBE BETTER IF I DON'T HAVE TO CALL bilin:
##             #lavg = ( (y - datay[0]) * data[1,0] + (datay[1] - y) * data[0,0] ) / (datay[1] - datay[0])
##             #ravg = ( (y - datay[0]) * data[1,1] + (datay[1] - y) * data[0,1] ) / (datay[1] - datay[0])
##             #return ( (x - datax[0]) * ravg + (datax[1] - x) * lavg ) / (datax[1] - datax[0])
##             dataout[iyout, ixout] = bilin(ixout, iyout, datain[iyin:iyin+2, ixin:ixin+2], bxout, byout)

##     #fits['data'] = fitsdata
##     print
##     #print "Writing out .fits file ", outfits, "...\n"
##     savefits(dataout, outfits)
##     #fitsio.writefits(outfits,fits)

def loadpixelscale(image):
    if os.path.exists('temp.txt'):
        os.remove('temp.txt')
    print('imsize ' + capfile(image, '.fits') + ' > temp.txt')
    os.system('imsize ' + capfile(image, '.fits') + ' > temp.txt')
    s = loadfile('temp.txt')[0]
    if s.find('"/pix') == -1:
        print('PIXEL SCALE NOT GIVEN IN IMAGE HEADER OF', capfile(image, '.fits'))
        pixelscale = 0
    else:
        s = s.split('"/pix')[0]
        s = s.split('/')[1]
        pixelscale = float(s)
    os.remove('temp.txt')
    return pixelscale


# p.14 postscript, native gtk and native wx do not support alpha or antialiasing.
# You can create an arbitrary number of axes images inside a single axes, and these will be composed via alpha blending. However, if you want to blend several images, you must make sure that the hold state is True and that the alpha of the layered images is less than 1.0; if alpha=1.0 then the image on top will totally obscure the images below. Because the image blending is done using antigrain (regardless of your backend choice), you can blend images even on backends which don't support alpha (eg, postscript). This is because the alpha blending is done in the frontend and the blended image is transferred directly to the backend as an RGB pixel array. See Recipe 9.4.2 for an example of how to layer images.

Ellipse = matplotlib.patches.Ellipse

def ellpatch1(x, y, w, h, ang, fc, ec, alpha=1, zorder=1, fill=1):
    patch = Ellipse((x,y), w, h, ang*180/pi, fill=fill)
    patch.set_fc(fc)
    patch.set_ec(ec)
    patch.set_alpha(alpha)
    patch.set_zorder(zorder)
    return patch

def polypatch1(x, y, fc, ec, alpha=1, zorder=1):
    patch = Polygon(list(zip(x,y)))
    patch.set_fc(fc)
    patch.set_ec(ec)
    patch.set_alpha(alpha)
    patch.set_zorder(zorder)
    return patch

def rectpatch1(x, y, dx, dy, fc, ec, alpha=1, zorder=1):
    patch = Rectangle((x-dx, y-dy), 2*dx, 2*dy)
    patch.set_fc(fc)
    patch.set_ec(ec)
    patch.set_alpha(alpha)
    patch.set_zorder(zorder)
    return patch

# /Library/Frameworks/Python.framework/Versions/Current/lib/python2.5/site-packages/matplotlib/contour.py
def contour_set(CS, alpha=None, zorder=None, color=None, sh=1):
    """contourf alpha doesn't work!
    CS = contourf(x, y, z)
    contour_set_alpha(CS, 0.5)"""
    for collection in CS.collections:
        if alpha  != None: collection.set_alpha(alpha)
        if zorder != None: collection.set_zorder(zorder)
        if color != None: collection.set_color(color)
    if sh: show()

def contour_set_alpha(CS, alpha, sh=1):
    """contourf alpha doesn't work!
    CS = contourf(x, y, z)
    contour_set_alpha(CS, 0.5)"""
    for collection in CS.collections:
        collection.set_alpha(alpha)
    if sh: show()

def zlabel(s, vert=False, fontsize=None, x=0.88, **other):
    """Places label on colorbar"""
    if vert:
        if fontsize == None:
            fontsize = 14
        figtext(x, 0.5, s, rotation='vertical', va='center', fontsize=fontsize, **other)
    else:
        figtext(x, 0.5, s, **other)

def plotsort(x, **other):
    i = arange(len(x)) / (len(x) - 1.)
    plot(i, sort(x), **other)

def makelines(x1, y1, x2, y2):
    """(x1, y1), (x2, y2) LIST ALL CONNECTIONS
    CONNECT THE DOTS AND RETURN LISTS OF LISTS: x, y"""
    n = len(x1)
    
    i = 0
    j = i
    while (x1[i+1] == x2[i]) and (y1[i+1] == y2[i]):
        i += 1
        if i > n - 2:
            break
    
    x = x1[j:i+1].tolist()
    y = y1[j:i+1].tolist()
    x.append(x2[i])
    y.append(y2[i])
    
    xx = [x]
    yy = [y]
    
    while i < n-2:
        i += 1
        j = i
        while (x1[i+1] == x2[i]) and (y1[i+1] == y2[i]):
            i += 1
            if i > n - 2:
                break
            
        x = x1[j:i+1].tolist()
        y = y1[j:i+1].tolist()
        x.append(x2[i])
        y.append(y2[i])
        #
        xx.append(x)
        yy.append(y)
    
    return xx, yy

#################################
# ZOOM IN ON DATASET
# SEE ~/glens/h0limis/results4ae.py

def xyrcut(xr, yr, i, fac=0.8):
    dx = p2p(xr)
    dy = p2p(yr)
    if i == 0:    # RIGHT
        xr = xr[0], xr[0] + fac * dx
    elif i == 1:  # LEFT
        xr = xr[1] - fac * dx, xr[1]
    elif i == 2:  # TOP
        yr = yr[0], yr[0] + fac * dy
    elif i == 3:  # BOTTOM
        yr = yr[1] - fac * dy, yr[1]
    return xr, yr

def catxyrcut(fac, ct, xr, yr, i, justn=0):
    if fac <= 0:
        return inf
    elif fac >= 1:
        return 0
    
    xr2, yr2 = xyrcut(xr, yr, i, 1-fac)
    ct2 = ct.between(xr2[0], 'x', xr2[1])
    ct2 = ct2.between(yr2[0], 'y', yr2[1])
    #print ct.len()
    #print ct2.len()
    
    n = (ct.len() - ct2.len()) / float(ct.len())
    n = divsafe(fac, n)
    if justn:
        return n
    else:
        return n, xr2, yr2, ct2

def funcy(x, func, args=(), y=0):
    #print x, func(x, *args), 'result'
    return abs(func(x, *args) - y)


def zoom(x, y, nfac=30, fac=0.2, margin=0.02):
    cat = VarsClass()
    cat.add('x', x)
    cat.add('y', y)
    cat.updatedata()
    
    xr = minmax(cat.x)
    yr = minmax(cat.y)
    
    for i in range(4):
        cat2 = cat
        n = 1e30
        while n > nfac:
            cat = cat2
            xr = minmax(cat.x)
            yr = minmax(cat.y)
            n, xr2, yr2, cat2 = catxyrcut(fac, cat, xr, yr, i)
            print(i, n, cat2.len())
    
    xlim(prange(xr, margin=margin))
    ylim(prange(yr, margin=margin))

#################################

def hline(v=0, c='k', ls='-', **other):
    """HORIZONTAL LINE THAT ALWAYS SPANS THE AXES"""
    axhline(v, c=c, ls=ls)

yline = hline

def vline(v=0, c='k', ls='-', **other):
    """VERTICAL LINE THAT ALWAYS SPANS THE AXES"""
    axvline(v, c=c, ls=ls, **other)

xline = vline

def axlines(x=0, y=0, c='k', ls='-', **other):
    """VERTICAL LINE THAT ALWAYS SPANS THE AXES"""
    axvline(x, c=c, ls=ls, **other)
    axhline(y, c=c, ls=ls, **other)

# from MLab_coe:
def singlevalue(x):
    """IS x A SINGLE VALUE?  (AS OPPOSED TO AN ARRAY OR LIST)"""
    return type(x) in [NoneType, float, float32, float64, int, int0, int8, int16, int32, int64]  # THERE ARE MORE TYPECODES IN Numpy

# log x and/or y axes with nice tick labels
# formatter = FuncFormatter(log_10_product)
# FOR FURTHER USE INSTRUCTIONS, SEE e.g.,
#  ~/LensPerfect/A1689/analysis/NFWfitWSplot.py
# http://www.thescripts.com/forum/thread462268.html
def log_10_product(x, pos):
    """The two args are the value and tick position.
    Label ticks with the product of the exponentiation"""
    #return '%1i' % (x)
    ndec1 = ndec(x)
    if ndec1 == 0:
        format = '%d'
    else:
        format = '%%.%df' % ndec1
    #print format, x
    return format % (x)

def savepdf(figname, saveeps=1):
    if figname[:-4] == '.pdf':
        figname = figname[:-4]
    savefig(figname+'.eps')
    #os.system('epstopdf %s.eps' % figname)
    os.system('pstopdf %s.eps' % figname)
    if not saveeps:
        os.remove(figname+'.eps')

def savepngpdf(figname, saveeps=1):
    if len(figname) > 4:
        if figname[-4] == '.':
            figname = figname[:-4]
    savefig(figname+'.png')
    savepdf(figname, saveeps=saveeps)

def savepng(figname):
    if len(figname) > 4:
        if figname[-4] == '.':
            figname = figname[:-4]
    savefig(figname+'.png')

def ploterrorbars1(x, y, dy, ymax=None, color='k', xfac=1, **other):
    if ymax == None:
        ymin = y - dy
        ymax = y + dy
    else:
        ymin = dy
        ymax = ymax
    
    dx = 0.005 * xfac * (xlim()[1] - xlim()[0])
    itemp = isinteractive()
    xtemp = xlim()
    ytemp = ylim()
    ioff()
    for i in range(len(x)):
        plot([x[i], x[i]], [ymin[i], ymax[i]], color=color, **other)
        plot([x[i]-dx, x[i]+dx], [ymax[i], ymax[i]], color=color, **other)
        plot([x[i]-dx, x[i]+dx], [ymin[i], ymin[i]], color=color, **other)
    
    if itemp:
        ion()
        show()
    
    xlim(xtemp[0], xtemp[1])
    ylim(ytemp[0], ytemp[1])

def ploterrorbars(x, y, dy, ymax=None, color='k', xfac=1, ax=None, **other):
    if ax == None:
        ax = gca()

    if ymax == None:
        ymin = y - dy
        ymax = y + dy
    else:
        ymin = dy
        ymax = ymax
    
    dx = 0.005 * xfac * (xlim()[1] - xlim()[0])
    itemp = isinteractive()
    xtemp = xlim()
    ytemp = ylim()
    ioff()
    for i in range(len(x)):
        ax.plot([x[i], x[i]], [ymin[i], ymax[i]], color=color, **other)
        ax.plot([x[i]-dx, x[i]+dx], [ymax[i], ymax[i]], color=color, **other)
        ax.plot([x[i]-dx, x[i]+dx], [ymin[i], ymin[i]], color=color, **other)
    
    if itemp:
        ion()
        show()
    
    ax.set_xlim(xtemp[0], xtemp[1])
    ax.set_ylim(ytemp[0], ytemp[1])

#################################
# ABILITY AVAILABLE IN RECENT VERSION OF matplotlib


def setaspect(xr=None, yr=None, aspect=5/7.):
    xr = xlim()
    yr = ylim()
    
    dx = xr[1] - xr[0]
    dy = yr[1] - yr[0]
    
    if dy/dx < aspect:
        yr = mean(yr) + dx * aspect * (arange(2)-0.5)
    else:
        xr = mean(xr) + dy / aspect * (arange(2)-0.5)
    
    xlim(xr[0], xr[1])
    ylim(yr[0], yr[1])

def xlim(lo=None,hi=None):
    if lo == None and hi == None:
        return xlim1()
    else:
        if singlevalue(lo):
            lo1, hi1 = xlim1()
            if lo == None:
                lo = lo1
            if hi == None:
                hi = hi1
        else:
            lo, hi = lo
        xlim1(lo,hi)

def ylim(lo=None,hi=None):
    if lo == None and hi == None:
        return ylim1()
    else:
        if singlevalue(lo):
            lo1, hi1 = ylim1()
            if lo == None:
                lo = lo1
            if hi == None:
                hi = hi1
        else:
            lo, hi = lo
        ylim1(lo,hi)

#################################


def len0(x):
    try:
        n = len(x)
    except:
        n = 0
    return n

# PROBLEM WITH THIS IS I CAN'T DO plot(3, 'o')
#from pylab import plot as _plot
#def plot(x, **other):
#    if not len0(x):
#       x = [x]
#    if 'y' not in other.keys():
#       _plot(x, **other)
#    else:
#       y = other['y']
#       if not len0(y):
#           y = [y]
#       _plot(x, y, **other)

#plot(1)
# FROM THE matplotlib MANUAL: "Because most GUIs have a mainloop, they become unresponsive to input outside of their mainloop once they are launched."  
# THEY SUGGEST SOME WORKAROUNDS, BUT THIS WORKS, TOO:
def killplot():
    plot([1])
    title('KILL THIS WINDOW!')
    show()
    # IMAGE WILL DISPLAY, BUT YOU WON'T BE ABLE TO USE THE Python SHELL
    # KILL THIS WINDOW, THEN...
    ioff()
    
if 0:
    killplot()


# WHEN YOU show() THE NEXT PLOT,
# YOU'LL STILL BE ABLE TO USE THE Python SHELL
# PLUS, THE PLOT WILL BE AUTOMATICALLY UPDATED WITH EACH COMMAND YOU GIVE

clear = cla  # cla() CLEAR PLOT



def closefig(num=None):
    if num == None:
        closecurrentfig()
    else:
        figure(num)
        closecurrentfig()

figclose = closefig

def smallfig(num, fac=2, reopen=0):
    if reopen:
        closefig(num)
    
    figure(num, figsize=(8/fac, 6/fac))

# now x, y ranges can be input: lo & hi dictate axes
def showarr(a, showcolorbar=1, nan=None, valrange=[None, None], sq=0, 
            cmap='jet', x=None, y=None, cl=1, sh=1):
    if valrange[0] or valrange[1]:
        if valrange[0] == None:
            valrange[0] = min(a)
        if valrange[1] == None:
            valrange[1] = max(a)
        a = clip2(a, valrange[0], valrange[1])
    if nan != None:
        a = where(isnan(a), nan, a)
    if cl:
        clf()
    ioff()
    if x != None and y != None:
        xlo, xhi = minmax(x)
        ylo, yhi = minmax(y)
        dx = (xhi - xlo) / (len(x) - 1.)
        dy = (yhi - ylo) / (len(y) - 1.)
        extent = (xlo-dx/2., xhi+dx/2., ylo-dy/2., yhi+dy/2.)
    else:
        extent = None
    cmap = cm.get_cmap(cmap)
    aspect = ['auto', 1][sq]  # 'preserve'
    im = imshow(a, origin='lower', interpolation='nearest', cmap=cmap,
                aspect=aspect, extent=extent)
    if showcolorbar:
        colorbar()
    if sh:
        show()
        ion()
    return im

def showxyz(x, y, z, **other):
    showarr(z, x=x, y=y, **other)

def rectangle(lolimits, hilimits, fillit=0, **other):
    [xmin,ymin] = lolimits
    [xmax,ymax] = hilimits
    if not fillit:
        return plot([xmin, xmin, xmax, xmax, xmin], [ymin, ymax, ymax, ymin, ymin], **other)
    else:
        if 'color' in list(other.keys()):
            color = other['color']
            del other['color']
            color = color2hex(color)
            return fill([xmin, xmin, xmax, xmax, xmin], [ymin, ymax, ymax, ymin, ymin], color, **other)
        else:
            return fill([xmin, xmin, xmax, xmax, xmin], [ymin, ymax, ymax, ymin, ymin], **other)

# FIX THE AXES ON A PLOT OF AN ARRAY
def retick(lo, hi, N, ndec=1, ytx=None):
    N = N - 1.

    if ytx == None:
        ylocs = arange(0, N+.001, N/4.)
        ytx = ylocs / float(N) * (hi - lo) + lo
    else:
        ylocs = (ytx - lo) * N / float(hi - lo)
    
    ytxs = []
    for ytick in ytx:
        format = '%%.%df' % ndec
        ytxs.append(format % ytick)
        #ytxs.append('%.1f' % ytick)
    
    yticks(ylocs, ytxs)
    xticks(ylocs, ytxs)
    
    xlim(0,N)
    ylim(0,N)
    show()

# FIX ONE AXIS ON A PLOT OF AN ARRAY
def retick1(x, axs, ntx=4, ndec=1, ytx=None, relim=True):
    N = len(x) - 1
    lo, hi = minmax(x)
    
    if ytx == None:
        ylocs = arange(0, N+.001, N/float(ntx))
        ytx = ylocs / float(N) * (hi - lo) + lo
    else:
        ylocs = (ytx - lo) * N / float(hi - lo)
    
    ytxs = []
    for ytick in ytx:
        format = '%%.%df' % ndec
        ytxs.append(format % ytick)
        #ytxs.append('%.1f' % ytick)
    
    if axs == 'x':
        p = xticks(ylocs, ytxs)
        if relim: p = xlim(-0.5, N+0.5)
    elif axs == 'y':
        p = yticks(ylocs, ytxs)
        if relim: p = ylim(-0.5, N+0.5)
    else:
        print('Sorry, which axis? --retick1')
        
    show()

def retick2(xlo, xhi, Nx, ylo, yhi, Ny, Nxtx=4, Nytx=4, ndec=1):
    Nx = Nx - 1.
    Ny = Ny - 1.
    
    xlocs = arange(0, Nx+.001, Nx/float(Nxtx))
    xtx = xlocs / float(Nx) * (xhi - xlo) + xlo
    xtxs = []
    for xtick in xtx:
        format = '%%.%df' % ndec
        xtxs.append(format % xtick)
    
    ylocs = arange(0, Ny+.001, Ny/float(Nytx))
    ytx = ylocs / float(Ny) * (yhi - ylo) + ylo
    ytxs = []
    for ytick in ytx:
        format = '%%.%df' % ndec
        ytxs.append(format % ytick)
    
    p = xticks(xlocs, xtxs)
    p = yticks(ylocs, ytxs)
    
    #xlim(0,Nx)
    #ylim(0,Ny)
    xlim(-0.5,Nx+0.5)
    ylim(-0.5,Ny+0.5)
    show()

def retick3(x, y, Nxtx=4, Nytx=4, ndec=1):
    xlo, xhi = minmax(x)
    ylo, yhi = minmax(y)
    Nx = len(x)
    Ny = len(y)
    retick2(xlo, xhi, Nx, ylo, yhi, Ny, Nxtx=Nxtx, Nytx=Nytx, ndec=ndec)

def reticklog1(x, axs, relim=True):
    """Automatcially places tickmarks at log intervals (0.001, 0.01, etc.)
    x contains data for either the x or y axis
    axs = 'x' or 'y' (which axis?)"""
    eps = 1e-7
    lo = min(log10(x))
    hi = max(log10(x))
    #lo = floor(lo+eps)
    #tx = arange(lo,hi+eps,1)
    tx = multiples(lo,hi)
    txs = []
    for logx in tx:
        n = 10 ** logx
        tx1 = '%g' % n
        tx1 = string.replace(tx1, '-0', '-')
        tx1 = string.replace(tx1, '+0', '+')
        txs.append(tx1)
    
    ntx = len(tx)
    nx = len(x) - 1
    #locs = mgrid[0:nx:ntx*1j]
    locs = interp(tx, array([lo, hi]), array([0, nx]))
    if axs == 'x':
        p = xticks(locs, txs)
        if relim: p = xlim(-0.5, nx+0.5)
    elif axs == 'y':
        p = yticks(locs, txs)
        if relim: p = ylim(-0.5, nx+0.5)
    else:
        print('Sorry, which axis? --reticklog1')

def reticklog(x, y, relim=True):
    reticklog1(x, 'x', relim=relim)
    reticklog1(y, 'y', relim=relim)

def fillbetween(x1, y1, x2, y2, **other):
    # MAKE SURE IT'S NOT A LIST, THEN IN CASE IT'S A numpy ARRAY, CONVERT TO LIST, THEN CONVERT TO numarray ARRAY
    if type(x1) != list:
        x1 = x1.tolist()
    if type(y1) != list:
        y1 = y1.tolist()
    if type(x2) != list:
        x2 = x2.tolist()
    if type(y2) != list:
        y2 = y2.tolist()
    x = x1[:]
    x[len(x):] = x2[::-1]
    y = y1[:]
    y[len(y):] = y2[::-1]
    return fill(x, y, **other)

def sqrange(x, y, margin=None, set=1):
    xmax = max(abs(array(x)))
    ymax = max(abs(array(y)))
    xymax = max([xmax, ymax])
    xyr = -xymax, xymax
    if margin != None:
        xyr = prange(xyr)
    if set:
        xlim(xyr)
        ylim(xyr)
    return xyr

def prange(x, xinclude=None, margin=0.05):
    """RETURNS GOOD RANGE FOR DATA x TO BE PLOTTED IN.
    xinclude = VALUE YOU WANT TO BE INCLUDED IN RANGE.
    margin = FRACTIONAL MARGIN ON EITHER SIDE OF DATA."""
    xmin = min(x)
    xmax = max(x)
    if xinclude != None:
        xmin = min([xmin, xinclude])
        xmax = max([xmax, xinclude])
    
    dx = xmax - xmin
    if dx:
        xmin = xmin - dx * margin
        xmax = xmax + dx * margin
    else:
        xmin = xmin - margin
        xmax = xmax + margin
    return [xmin, xmax]

def prangelog(x, xinclude=None, margin=0.05):
    """RETURNS GOOD RANGE FOR DATA x TO BE PLOTTED IN.
    xinclude = VALUE YOU WANT TO BE INCLUDED IN RANGE.
    margin = FRACTIONAL MARGIN ON EITHER SIDE OF DATA."""
    xmin = min(x)
    xmax = max(x)
    if xinclude != None:
        xmin = min([xmin, xinclude])
        xmax = max([xmax, xinclude])
    
    fac = xmax / xmin
    xmin = xmin / (fac ** margin)
    xmax = xmax * (fac ** margin)
    
    return [xmin, xmax]

# whiskerplot.py
def vectorplot(vx, vy, x=[], y=[], xyfactor=[], heads=1, clear=1, color='k', **other):
    #print 'vx, vy', minmax(ravel(vx)), minmax(ravel(vy))
    nmax = 10
    if clear:
                clf()
    ioff()  # DON'T UPDATE PLOT WITH EACH ADDITIION (WAIT UNTIL END...)
    # IF TOO MANY VECTORS, SAMPLE:
    if len(vx.shape) == 2:
        ny, nx = vx.shape
        if (ny > nmax) or (nx > nmax):
            yall, xall = indices((ny, nx))
            nyred = min([ny, nmax])
            nxred = min([nx, nmax])
            dy = ny / float(nyred)
            dx = nx / float(nxred)
            vxred = zeros((nyred, nxred), float)
            vyred = zeros((nyred, nxred), float)
            for iy in range(nyred):
                y = int(iy * dy)
                for ix in range(nxred):
                    x = int(ix * dx)
                    vxred[iy,ix] = vx[y,x]
                    vyred[iy,ix] = vy[y,x]
                    #print (iy,ix), (y,x)
                    #print (vy[y,x], vx[y,x]), (vyred[iy,ix], vxred[iy,ix])
            vy = vyred
            vx = vxred
        y, x = indices(vx.shape)
        vx = ravel(vx)
        vy = ravel(vy)
        x = ravel(x)
        y = ravel(y)
    #print 'RED vx, vy', minmax(ravel(vx)), minmax(ravel(vy))
    try:
        xfactor, yfactor = xyfactor
    except:
        if xyfactor != []:
            xfactor = yfactor = xyfactor
        else:
            xfactor = 1.
            yfactor = 1.
    #ny, nx = vx.shape
    vr = hypot(vx, vy)
    #print max(ravel(vx))
    #print max(ravel(vy))
    #maxvr = max(ravel(vr))
    maxvr = max(vr)
    vx = vx / maxvr * 0.1 * ptp(x) * xfactor
    vy = vy / maxvr * 0.1 * ptp(y) * yfactor
    #print maxvr, xfactor, yfactor
    for i in range(len(vx)):
        xo = x[i]
        yo = y[i]
        dx = vx[i]
        dy = vy[i]
        #print dx, dy
        #pause('TO PLOT...')
        dr = hypot(dx, dy)
        # IF IT'S TOO SMALL, DRAW A DOT
        if dr < xfactor / 100.:
            plot([xo], [yo], 'o', markerfacecolor=color)
            continue
        # OTHERWISE, DRAW A STICK
        plot([xo-dx, xo+dx], [yo-dy, yo+dy], color=color)
        if heads:
            # NOW FOR THE HEAD OF THE ARROW
            hx = -dx + -dy
            hy =  dx + -dy
            #hr = hypot(hx, hy)
            #hx = 0.1 * hx / hr * xfactor
            #hy = 0.1 * hy / hr * yfactor
            hx = 0.2 * hx
            hy = 0.2 * hy
            plot([xo+dx, xo+dx+hx], [yo+dy, yo+dy+hy], color=color)
            plot([xo+dx, xo+dx-hy], [yo+dy, yo+dy+hx], color=color)
    #
    show() # NOW SHOW THE PLOT
    ion()  # AND TURN INTERACTIVE MODE BACK ON


def atobplot(xa, ya, xb, yb, color='k', linetype='arrow', showplot=1, **other):
    """DRAWS LINES FROM a TO b"""
    n = len0(xa)
    nb = len0(xb)
    if not n and not nb:
        n = 1
        xa = [xa]
        ya = [ya]
        xb = [xb]
        yb = [yb]
    elif not n:
        n = nb
        xa = [xa] * n
        ya = [ya] * n
    elif not nb:
        xb = [xb] * n
        yb = [yb] * n
    isint = isinteractive()
    ioff()
    for i in range(n):
        plot([xa[i], xb[i]], [ya[i], yb[i]], color=color, **other)
        if linetype=='arrow':
            # NOW FOR THE HEAD OF THE ARROW
            dx = xb[i] - xa[i]
            dy = yb[i] - ya[i]
            hx = -dx + -dy
            hy =  dx + -dy
            #hr = hypot(hx, hy)
            #hx = 0.1 * hx / hr * xfactor
            #hy = 0.1 * hy / hr * yfactor
            hx = 0.1 * hx
            hy = 0.1 * hy
            plot([xb[i], xb[i]+hx], [yb[i], yb[i]+hy], color=color, **other)
            plot([xb[i], xb[i]-hy], [yb[i], yb[i]+hx], color=color, **other)
        if linetype=='xo':
                plot(xa, ya, 'x', markerfacecolor=color, markeredgecolor=color, **other)
                plot(xb, yb, 'o', markerfacecolor=color, **other)
    if isint and showplot:
        ion()
        show()

#atheta = interp2Dtheta(x, y, xim, yim, athetaim, pow=3.)

#vy, vx = indices((3,3))
#vectorplot(vx, vy)

#vectorplot(ax, ay)


        #yred, xred = indices((nyred, nxred))
        #yred = yred * (ny / float(nyred))
        #xred = xred * (nx / float(nxred))
        #vx = interpxy(xred, yred, xall, yall, vx)
        #vy = interpxy(xred, yred, xall, yall, vy)


def colormaprgb(val, valrange=[0.,1.], cmap='jet', silent=0):
    if valrange != [0.,1.]:
        lo = float(valrange[0])
        hi = float(valrange[1])
        val = (val - lo) / (hi - lo)
    
    try:
        n = len(val)
    except:  # SINGLE VALUE
        val = array([val])
    
    cmapa = colormap_map[cmap]
    
    xa = arange(len(cmapa)) / float(len(cmapa)-1)
    ra, ga, ba = transpose(cmapa)
    
    r = interpn(val, xa, ra, silent)
    g = interpn(val, xa, ga, silent)
    b = interpn(val, xa, ba, silent)
    rgb = ravel(array([r, g, b]))
    return rgb

# colormap IS ALREADY BUILT IN, BUT NAMED scatter
# AND THIS WAY, YOU CAN ADD A colorbar()!
def colormap(x, y, z, showcolorbar=1, ticks=None, **other):
    x = ravel(x)
    y = ravel(y)
    z = ravel(z)
    scatter(x, y, c=z, **other)
    if showcolorbar:
        colorbar(ticks=ticks)

# SEE colormap.py
def colormap_obsolete(x, y, z, valrange=[None, None], cmap='jet', markersize=5):
    ioff()
    
    if valrange[0] == None:
        valrange[0] = min(z)
    if valrange[1] == None:
        valrange[1] = max(z)
    n = len(x)
    
    for i in range(n):
        if not (isNaN(z[i]) or (valrange[0] == valrange[1])):
            plot( [x[i]], [y[i]], 'o', markerfacecolor=colormaprgb(z[i], valrange=valrange), markersize=markersize )
        else:
            p.add(Point( x[i], y[i], 'ko', markersize=markersize ))
    
    colorbar()
    ion()
    show()


def test():
    x = arange(10)
    y = arange(10)
    z = arange(10)
    colormap(x, y, z)

#test()

# THEY MUST'VE WRITTEN THIS, BUT I CAN'T FIND IT!
def circle(xc, yc, r, n=100, **other):
    ang = arange(n+1) / float(n) * 2*pi
    x = xc + r * cos(ang)
    y = yc + r * sin(ang)
    plot(x, y, **other)

def circles(xc, yc, r, n=100, **other):
    if singlevalue(xc):
        xc = xc * ones(len(r))
    if singlevalue(yc):
        yc = yc * ones(len(r))
    for i in range(len(xc)):
        circle(xc[i], yc[i], r[i], n=n, **other)

# COULD ALSO USE matplotlib.patches.Ellipse
def ellipse(xc, yc, a, b, ang=0, n=100, **other):
    t = arange(n+1) / float(n) * 2*pi
    x = a * cos(t)
    y = b * sin(t)
    if ang:
        x, y = rotdeg(x, y, ang)
    x = x + xc
    y = y + yc
    plot(x, y, **other)

def ellipses(xc, yc, a, b, ang=None, n=100, **other):
    if ang==None:
        ang = zeros(len(a))
    for i in range(len(xc)):
        ellipse(xc[i], yc[i], a[i], b[i], ang[i], n=n, **other)

def texts(x, y, z, format='%d', showit=True, **other):
    isint = isinteractive()
    ioff()
    if singlevalue(z):
        z = z * ones(len(x))
    for i in range(len(x)):
        t = text(x[i], y[i], format % z[i], **other)
    if isint and showit:
        ion()
        show()

# FROM MLab_coe.py
# COLOR INPUT: fc = face color; ec = edge color
def bargraph(x, y, fill=1, color='black', zeroedges=1, **other):
    n = len(x)
    xx = repeat(x, 2)
    y = y.astype(float)
    z = array([0.])
    yy = repeat(y, 2)
    if zeroedges:
        yy = concatenate([z, repeat(y, 2), z])
    else:
        xx = xx[1:-1]
    zz = yy*0
    
    if fill:
        fc = color
        ec = ['', color, 'k'][fill]
        if 'fc' in list(other.keys()):
            fc = other['fc']
            del(other['fc'])
            if (fill == 1) and 'ec' not in list(other.keys()):
                ec = fc
        if 'ec' in list(other.keys()):
            ec = other['ec']
            del(other['ec'])
        return fillbetween(xx, yy, xx, zz, fc=fc, ec=ec, **other)
    else:
        return plot(xx, yy, color=color, **other)

#        p.add(FillBetween(xx, yy, xx, zz, color=color))
#        p = fillbetween(xx, yy, xx, zz, **other)
#        if fill == 2:
#            p0 = plot(xx, yy, color='k')

# ALL LINES EXTEND DOWN TO ZERO
def bargraph2(x, y, fill=1, **other):
    n = len(x)
    xx = repeat(x, 3)[1:-1]
    y = y.astype(float)
    z = array([0.])
    yy = zeros(len(xx))
    put(yy, arange(len(y))*3+1, y)
    put(yy, arange(len(y))*3+2, y)
    zz = yy*0
    
    if fill:
        fillbetween(xx, yy, xx, zz, color=color)
    else:
        plot(xx, yy, color=color, **other)

# ALL LINES EXTEND DOWN TO ZERO
# color input works
def bargraph3(x, y, fill=1, **other):
    n = len(x)
    xx = repeat(x, 3)[1:-1]
    y = y.astype(float)
    z = array([0.])
    yy = zeros(len(xx))
    put(yy, arange(len(y))*3+1, y)
    put(yy, arange(len(y))*3+2, y)
    zz = yy*0
    
    if fill:
        fillbetween(xx, yy, xx, zz)
    else:
        plot(xx, yy, **other)


def matrix_multiply(MM):
    """Multiplies a list of matrices: M[0] * M[1] * M[2]..."""
    P = MM[0]
    for M in MM[1:]:
        P = dot(P, M)
    return P

def sinn(x):
    """
    x < 0: sin
    x > 0: sinh
    """
    if x < 0:
        return sin(x)
    else:
        return sinh(x)

def multiples(lo, hi, x=1, eps=1e-7):
    """Returns an array of the multiples of x between [lo,hi] inclusive"""
    l = ceil((lo-eps)/x)*x
    return arange(l, hi+eps, x)

def multiples2(lohi, x=1, eps=1e-7):
    """Returns an array of the multiples of x between [lo,hi] inclusive"""
    lo, hi = lohi
    return multiples(lo, hi, x, eps)

def multipleslog(lo, hi):
    """Returns an array of the log multiples between [lo,hi] inclusive.
    That didn't make sense, but what I'm trying to say is:
    multipleslog(2, 30) = 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30"""
    loglo = log10(lo)
    loghi = log10(hi)
    ll = multiples(loglo, loghi)
    ll = concatenate([[loglo], ll, [loghi]])
    mm = []
    for i in range(len(ll)-1):
        lo = 10 ** ll[i]
        hi = 10 ** ll[i+1]
        ex = 10 ** floor(ll[i])
        m1 = multiples(lo, hi, ex)
        if len(mm):
            if close(m1[0], mm[-1]):
                m1 = m1[1:]
        mm = concatenate([mm, m1])
    return mm

def multiples2log(lohi):
    """Returns an array of the log multiples between [lo,hi] inclusive.
    That didn't make sense, but what I'm trying to say is:
    multipleslog(2, 30) = 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30"""
    lo, hi = lohi
    return multipleslog(lo, hi)


def onlyids(data, ids):
    """ALTERS ARRAY data TO INCLUDE ONLY NUMBERS IN ids
    ALL OTHER VALUES SET TO zero"""
    keys = arange(data.size)
    
    keysc = compress(data.flat, keys)
    valsc = compress(data.flat, data.flat)
    
    mask = zeros(data.size)
    for id in ids:
        ee = equal(valsc, id)
        mask = logical_or(mask, ee)
    
    keyscm = compress(mask, keysc)
    valscm = compress(mask, valsc)
    
    datanew = zeros(data.shape)
    datanew.put(keyscm, valscm)
    
    return datanew

def cliplohi(xlo, xhi, xmin, xmax):
    return max([xlo, xmin]), min([xhi, xmax])

def base(b, nums):
    """base(10, [1, 2, 3]) RETURNS 123"""
    if not isinstance(nums, list):
        nums = nums.tolist()
    nums.reverse()
    x = 0
    for i, num in enumerate(nums):
        x += array(num) * b**i
    return x

def strbegin(str, phr):  # coetools.py
    return str[:len(phr)] == phr

def minsec(x, format=(), precision=None):
    """
    CONVERTS decimal degrees/hours to degrees/hours : minutes : seconds
    minsec(13.52340987)
    minsec(13.52340987, ':')
    minsec(13.52340987, 'hms')
    minsec(13.52340987, 'dms')
    minsec(13.52340987, 'dms', 1)
    """
    f, i = math.modf(x)
    i = int(i)
    m = 60 * f
    s, m = math.modf(m)
    m = int(m)
    s = 60 * s
    if type(format) == str:
        if precision == None:
            s = '%f' % s
        else:
            fmt = '%%.%df' % precision
            s = fmt % s
        if strbegin(s, '60'):  # rounded up
            s = '0'
            m = m + 1
        m = '%d' % m
        if m == '60':  # rounded up
            m = '0'
            i += 1
        i = '%d' % i
        ims = (i,m,s)
        if len(format) == 1:
            out = string.join(ims, format)
        elif len(format) == 3:
            out = i+format[0] + m+format[1] + s+format[2]
    else:
        out = (i, m, s)
    return out

def sec2hms(x, precision=0, mpersist=True):
    """
    CONVERTS decimal seconds to hours : minutes : seconds
    """
    out = ''
    if x > 60:
        if x > 3600:
            h = int(x / 3600)
            out = '%d:' % h
            x = x - 3600 * h
        m = int(x / 60)
        out += '%d:' % m
        x = x - 60 * m
    elif mpersist:
        out = '0:'
    if precision == None:
        fmt = '%g'
    elif precision == 0:
        fmt = '%d'
    else:
        fmt = '%%.%df' % precision
    s = fmt % x
    if (x < 10) and mpersist:
        s = '0' + s
    out += s
    return out

def prange(x, xinclude=None, margin=0.05):
    """RETURNS GOOD RANGE FOR DATA x TO BE PLOTTED IN.
    xinclude = VALUE YOU WANT TO BE INCLUDED IN RANGE.
    margin = FRACTIONAL MARGIN ON EITHER SIDE OF DATA."""
    xmin = min(x)
    xmax = max(x)
    if xinclude != None:
        xmin = min([xmin, xinclude])
        xmax = max([xmax, xinclude])
    
    dx = xmax - xmin
    if dx:
        xmin = xmin - dx * margin
        xmax = xmax + dx * margin
    else:
        xmin = xmin - margin
        xmax = xmax + margin
    return [xmin, xmax]

def minmax(x, range=None):
    if range:
        lo, hi = range
        good = between(lo, x, hi)
        x = compress(good, x)
    return min(x), max(x)

def rescale(x, lohi):
    lo, hi = lohi
    xlo, xhi = minmax(x)
    dx = xhi - xlo
    dy = hi - lo
    y = x / dx * dy + lo
    return y

def inrange(x, r):
    lo, hi = minmax(r)
    return between(lo, x, hi)

def pairs(x):
    p = []
    for i in range(len(x)):
        for j in range(i+1, len(x)):
            p.append((x[i], x[j]))
    return p

def Psig(P, nsigma=1):
    """(ir, il) bound central nsigma of P
    -- edges contain equal amounts of P"""
    Pn = P / total(P)
    g = gausst(nsigma)
    Pl = cumsum(Pn)
    Pr = cumsum(Pn[::-1])
    n = len(P)
    i = arange(n)
    il = interp(g, Pl, i)
    ir = interp(g, Pr, i)
    ir = n - ir
    return il, ir

def xsig(x, P, nsigma=1):
    return p2p(take(x, Psig(P, nsigma))) / 2.

def gaussin(nsigma=1):
    """FRACTION WITHIN nsigma"""
    return erf(nsigma / sqrt(2))

def gaussp(nsigma=1):
    """FRACTION INCLUDED UP TO nsigma"""
    return 0.5 + gaussin(nsigma) / 2.

def gaussbtw(nsig1, nsig2):
    """FRACTION BETWEEN nsig1, nsig2"""
    return abs(gaussp(nsig2) - gaussp(nsig1))

#gaussbtw(0, 3)


sigma = gaussin

def gausst(nsigma=1):
    """FRACTION IN TAIL TO ONE SIDE OF nsigma"""
    return 1 - gaussp(nsigma)

def pick(x):
    n = len(x)
    i = random_integers(n)
    return x[i-1]

def randrange(N=1):
    return (2 * random(N) - 1)

def randrange2(lo, hi, N=1):
    return ((hi - lo) * random(N) + lo)

class PDraw:
    def __init__(self, x, P):
        self.x = x
        self.P = P
        self.Pcum = cumsum(P)
        self.N = self.Pcum[-1]
    def draw(self, n=1):
        r = self.N * random(n)
        i = searchsorted(self.Pcum, r)
        return take(self.x, i)

def hypotsq(dx, dy):
    return dx**2 + dy**2

def hypotn(x):
    return sqrt(sum(x**2))

def hypotnn(*x):
    return hypotn(array(x))

#hypotnn(3, 4, 5)

def hypotxy(x1, y1, x2, y2):
    return hypot(x1-x2, y1-y2)

def subtend(x1, y1, x2, y2):
    """ANGLE SUBTENDED BY TWO VECTORS (wrt THE ORIGIN)"""
    # v1 (dot) v2 = |v1| |v2| cos(theta)
    # d = r1 r2 cos(theta)
    d = dot([x1, y1], [x2, y2])
    r1 = hypot(x1, y1)
    r2 = hypot(x2, y2)
    costheta = d / (r1 * r2)
    theta = arccos(costheta)
    return theta

def subtends(x, y):
    n = len(x)
    dd = []
    for i in range(n-1):
        for j in range(i+1,n):
            dd.append(subtend(x[i], y[i], x[j], y[j]))
    return array(dd)

def distances(x, y):
    n = len(x)
    dd = []
    for i in range(n-1):
        for j in range(i+1,n):
            dd.append(hypot(x[i]-x[j], y[i]-y[j]))
    return array(dd)

def nrange(x, n=100):
    """n EQUALLY-SPACED SAMPLES ON THE RANGE OF x"""
    return arange(n) / (n-1.) * (max(x) - min(x)) + min(x)

def range01(n=100):
    """n EQUALLY-SPACED SAMPLES ON THE RANGE OF [0,1]"""
    return arange(n) / (n-1.)

def middle(x):
    return (max(x) + min(x)) / 2.

def within(A, xc, yc, ro, yesorno=0):  # --DC
    """RETURNS WHETHER EACH PIXEL OF AN ARRAY IS WITHIN A CIRCLE
    DEFINED ON THE ARRAY'S COORDINATES.
    FRACTIONAL MEMBERSHIP IS ALSO ESTIMATED
    BY THE FRACTION OF THE BOX CROSSED BY THE CIRCLE AT THAT ANGLE.
    IT'S LIKE ANTI-ALIASING.
    THESE FRACTIONS ARE SLIGHTLY OVERESTIMATED
    BUT ARE AN IMPROVEMENT OVER NOT USING THEM AT ALL!
    TO TURN OFF FRACTIONS AND JUST RETURN True/False, SET yesorno=1"""
    ny, nx = A.shape
    a = ones((ny,nx))
    y = arange(ny)
    x = arange(nx)
    x, y = meshgrid(x, y)
    x = x-xc + 0.
    y = y-yc + 0.
    r = hypot(x,y)
    xy = abs(divsafe(x, y, nan=0))
    yx = abs(divsafe(y, x, nan=0))
    m = min([xy, yx])
    dr = hypot(1, m)  # = 1 ON AXES, sqrt(2) ON DIAGONALS
    
    if (ro - xc > 0.5) or (ro - yc > 0.5) \
            or (ro + xc > nx - 0.5) or (ro + yc > ny - 0.5):
        print('WARNING: CIRCLE EXTENDS BEYOND BOX IN MLab_coe.within')
    
    if yesorno:
        v = less_equal(r, ro)  # TRUE OR FALSE, WITHOUT FRACTIONS
    else:
        v = less_equal(r, ro-0.5*dr) * 1
        v = v + between(ro-0.5*dr, r, ro+0.5*dr) * (ro+0.5*dr - r) / dr
    
    #if showplot:  matplotlib NOT LOADED IN MLab_coe
    if 0:
        matshow(v)
        circle(xc+0.5, yc+0.5, ro, color='k', linewidth=2)
    
    return v

#def sumwithin(A, xc, yc, ro, showplot=0):
#    return total(A * within(A, xc, yc, ro, showplot=showplot))

def sumwithin(A, xc, yc, ro):
    """RETURNS SUM OF ARRAY WITHIN CIRCLE DEFINED ON ARRAY'S COORDINATES"""
    return total(A * within(A, xc, yc, ro))

def floatin(x, l, ndec=3):
    """IS x IN THE LIST l?
    WHO KNOWS WITH FLOATING POINTS!"""
    x = int(x * 10**ndec + 0.1)
    l = (array(l) * 10**ndec + 0.1).astype(int).tolist()
    return x in l

def floatindex(x, l, ndec=3):
    """IS x IN THE LIST l?
    WHO KNOWS WITH FLOATING POINTS!"""
    x = int(x * 10**ndec + 0.1)
    l = (array(l) * 10**ndec + 0.1).astype(int).tolist()
    return l.index(x)

def integral(f, x1, x2):
    return quad(f, x1, x2)[0]

def magnify(a, n):
    """MAGNIFIES A MATRIX BY n
    YIELDING, FOR EXAMPLE:
    >>> magnify(IndArr(3,3), 2)
    001122
    001122
    334455
    334455
    667788
    667788
    """
    ny, nx = a.shape
    a = repeat(a, n**2)
    a = reshape(a, (ny,nx,n,n))
    a = transpose(a, (0, 2, 1, 3))
    a = reshape(a, (n*ny, n*nx))
    return a

def demagnify(a, n, func='mean'):
    """DEMAGNIFIES A MATRIX BY n
    YIELDING, FOR EXAMPLE:
    >>> demagnify(magnify(IndArr(3,3), 2), 2)
    012
    345
    678
    """
    ny, nx = array(a.shape) / n
    a = a[:ny*8,:nx*8]  # Trim if not even multiples
    a = reshape(a, (ny, n, nx, n))
    a = transpose(a, (0, 2, 1, 3))
    a = reshape(a, (ny, nx, n*n))
    a = transpose(a, (2, 0, 1))
    exec('a = %s(a)' % func)
    return a

# Elementary Matrices

# zeros is from matrixmodule in C
# ones is from Numeric.py



def listo(x):
    if singlevalue(x):
        x = [x]
    return x

# ~/glens/h0limits/scatterrea.py
def insidepoly1(xp, yp, x, y):
    """DETERMINES WHETHER THE POINT (x, y)
    IS INSIDE THE CONVEX POLYGON DELIMITED BY (xp, yp)"""
    xp, yp = CCWsort(xp, yp)
    xp = xp.tolist()
    yp = yp.tolist()
    if xp[-1] != xp[0]:
        xp.append(xp[0])
        yp.append(yp[0])
    
    xo = mean(xp)
    yo = mean(yp)
    inpoly = 1
    xa = [xo, x]
    ya = [yo, y]
    for j in range(len(xp)-1):
        xb = xp[j:j+2]
        yb = yp[j:j+2]
        if linescross2(xa, ya, xb, yb):
            inpoly = 0
            break
    
    return inpoly

# ~/glens/h0limits/scatterrea.py
def insidepoly(xp, yp, xx, yy):
    """DETERMINES WHETHER THE POINTS (xx, yy)
    ARE INSIDE THE CONVEX POLYGON DELIMITED BY (xp, yp)"""
    xp, yp = CCWsort(xp, yp)
    xx = ravel(listo(xx))
    yy = ravel(listo(yy))
    inhull = []
    for i in range(len(xx)):
        if i and not (i % 10000):
            print('%d / %d' % (i, len(xx)))
        inhull1 = insidepoly1(xp, yp, xx[i], yy[i])
        inhull.append(inhull1)
    
    return array(inhull).astype(int)


# TESTED IN ~/glens/lenspoints/optdefl/sourceconstraints/testconvexhull.py
#  testinsidepoly() -- NEVER QUITE GOT THE TEST TO WORK HERE
def insidepolyshwag(xp, yp, xx, yy):
    """DETERMINES WHETHER THE POINTS (xx, yy)
    ARE INSIDE THE CONVEX POLYGON DELIMITED BY (xp, yp)"""
    xp, yp = CCWsort(xp, yp)  # NEEDED
    xp = xp.tolist()
    yp = yp.tolist()
    if xp[-1] != xp[0]:
        xp.append(xp[-1])  # SHOULD BE [0]
        yp.append(yp[-1])  # SHOULD BE [0]
    
    xo = mean(xp)
    yo = mean(yp)
    xx = ravel(listo(xx))
    yy = ravel(listo(yy))
    inhull = ones(len(xx)).astype(int)
    for i in range(len(xx)):
        if i and not (i % 10000):
            print('%d / %d' % (i, len(xx)))
        xa = [xo, xx[i]]
        ya = [yo, yy[i]]
        for j in range(len(xp)-2):
            xb = xp[j:j+2]
            yb = yp[j:j+2]
            if linescross2(xa, ya, xb, yb):
                inhull[i] = 0
                break
    
    return inhull

def testinsidepoly():
    #from numpy.random import random
    N = 40
    x = random(50) * N
    y = random(50) * N
    xh, yh = convexhull(x, y)
    zz = arange(N)
    xx, yy = meshgrid(zz, zz)
    xx = ravel(xx)
    yy = ravel(yy)
    inhull = insidepoly(xh, yh, xx, yy)
    figure(11)
    clf()
    plot(xh, yh)
    ioff()
    for i in range(len(XX)):
        color = ['r', 'g'][ininin[i]]
        p = plot([xx[i]], [yy[i]], 'o', mfc=color)
    
    show()

def p2p(x):  # DEFINED AS ptp IN MLab (BELOW)
    return max(x) - min(x)

def rotate(x, y, ang):
    """ROTATES (x, y) BY ang RADIANS CCW"""
    x2 = x * cos(ang) - y * sin(ang)
    y2 = y * cos(ang) + x * sin(ang)
    return x2, y2

def rotdeg(x, y, ang):
    """ROTATES (x, y) BY ang DEGREES CCW"""
    return rotate(x, y, ang/180.*pi)

def linefit(x1, y1, x2, y2):
    """y = mx + b FIT TO TWO POINTS"""
    if x2 == x1:
        m = Inf
        b = NaN
    else:
        m = (y2 - y1) / (x2 - x1)
        b = y1 - m * x1
    return m, b

def linescross(xa, ya, xb, yb):
    """
    DO THE LINES CONNECTING A TO B CROSS?
    A: TWO POINTS: (xa[0], ya[0]), (xa[1], ya[1])
    B: TWO POINTS: (xb[0], yb[0]), (xb[1], yb[1])
    DRAW LINE FROM A0 TO B0 
    IF A1 & B1 ARE ON OPPOSITE SIDES OF THIS LINE, 
    AND THE SAME IS TRUE VICE VERSA,
    THEN THE LINES CROSS
    """
    if xa[0] == xb[0]:
        xb = list(xb)
        xb[0] = xb[0] + 1e-10
    
    if xa[1] == xb[1]:
        xb = list(xb)
        xb[1] = xb[1] + 1e-10
    
    m0, b0 = linefit(xa[0], ya[0], xb[0], yb[0])
    ya1 = m0 * xa[1] + b0
    yb1 = m0 * xb[1] + b0
    cross1 = (ya1 > ya[1]) != (yb1 > yb[1])
    
    m1, b1 = linefit(xa[1], ya[1], xb[1], yb[1])
    ya0 = m1 * xa[0] + b1
    yb0 = m1 * xb[0] + b1
    cross0 = (ya0 > ya[0]) != (yb0 > yb[0])
    
    return cross0 and cross1

def linescross2(xa, ya, xb, yb):
    """
    DO THE LINES A & B CROSS?
    DIFFERENT NOTATION:
    LINE A: (xa[0], ya[0]) -> (xa[1], ya[1])
    LINE B: (xb[0], yb[0]) -> (xb[1], yb[1])
    DRAW LINE A
    IF THE B POINTS ARE ON OPPOSITE SIDES OF THIS LINE, 
    AND THE SAME IS TRUE VICE VERSA,
    THEN THE LINES CROSS
    """
    if xa[0] == xa[1]:
        xa = list(xa)
        xa[1] = xa[1] + 1e-10
    
    if xb[0] == xb[1]:
        xb = list(xb)
        xb[1] = xb[1] + 1e-10
    
    ma, ba = linefit(xa[0], ya[0], xa[1], ya[1])
    yb0 = ma * xb[0] + ba
    yb1 = ma * xb[1] + ba
    crossb = (yb0 > yb[0]) != (yb1 > yb[1])
    
    mb, bb = linefit(xb[0], yb[0], xb[1], yb[1])
    ya0 = mb * xa[0] + bb
    ya1 = mb * xa[1] + bb
    crossa = (ya0 > ya[0]) != (ya1 > ya[1])
    
    return crossa and crossb

def linescross2test():
    # from numpy.random import random
    xa = random(2)
    ya = random(2)
    xb = random(2)
    yb = random(2)
    
    figure(1)
    clf()
    plot(xa, ya)
    plot(xb, yb)
    title('%s' % linescross2(xa, ya, xb, yb))
    show()

def linescrosstest():
    # from random import random
    xa = random(), random()
    ya = random(), random()
    xb = random(), random()
    yb = random(), random()
    
    figure(1)
    clf()
    atobplot(xa, ya, xb, yb, linetype='')
    title('%s' % linescross(xa, ya, xb, yb))
    show()

def outside(x, y, xo, yo):
    """GIVEN 3 POINTS a, b, c OF A POLYGON 
    WITH CENTER xo, yo
    DETERMINE WHETHER b IS OUTSIDE ac,
    THAT IS, WHETHER abc IS CONVEX"""
    # DOES o--b CROSS a--c ?
    #      A--B       A--B
    xa, xb, xc = x
    ya, yb, yc = y
    xA = (xo, xa)
    yA = (yo, ya)
    xB = (xb, xc)
    yB = (yb, yc)
    return linescross(xA, yA, xB, yB)

# TESTED IN ~/glens/lenspoints/optdefl/sourceconstraints/testconvexhull.py
def convexhull(x, y, rep=1, nprev=0):
    """RETURNS THE CONVEX HULL OF x, y
    THAT IS, THE EXTERIOR POINTS"""
    x = x.astype(float)
    y = y.astype(float)
    x, y = CCWsort(x, y)
    xo = mean(x)
    yo = mean(y)
    x = x.tolist()
    y = y.tolist()
    dmax = max([p2p(x), p2p(y)])
    ngood = 0
    while ngood < len(x)+1:
        dx = x[1] - xo
        dy = y[1] - yo
        dr = hypot(dx, dy)
        dx = dx * dmax / dr
        dy = dy * dmax / dr
        x1 = xo - dx
        y1 = yo - dy
        if not outside(x[:3], y[:3], x1, y1):
            del x[1]
            del y[1]
        else: # ROTATE THE COORD LISTS
            x.append(x.pop(0))
            y.append(y.pop(0))
            ngood += 1
    
    x = array(x)
    y = array(y)
    
    # REPEAT UNTIL CONVERGENCE
    if (nprev == 0) or (len(x) < nprev):
        x, y = convexhull(x, y, nprev=len(x))
    
    if rep:
        x = concatenate((x, [x[0]]))
        y = concatenate((y, [y[0]]))
    
    return x, y

def gauss(r, sig=1., normsum=1):
    """GAUSSIAN NORMALIZED SUCH THAT AREA=1"""
    r = clip(r/float(sig), 0, 10)
    G = exp(-0.5 * r**2)
    G = where(less(r, 10), G, 0)
    if normsum:
        G = G * 0.5 / (pi * sig**2)
    return G

def gauss1(r, sig=1.):
    """GAUSSIAN NORMALIZED SUCH THAT PEAK AMPLITUDE = 1"""
    return gauss(r, sig, 0)

def atanxy(x, y, degrees=0):
    """ANGLE CCW FROM x-axis"""
    theta = arctan(divsafe(y, x, inf=1e30, nan=0))
    theta = where(less(x, 0), theta + pi, theta)
    theta = where(logical_and(greater(x, 0), less(y, 0)), theta + 2*pi, theta)
    if degrees:
        theta = theta * 180. / pi
    return theta


def chebyshev(x,n):
    if n == 0:
        return x ** 0
    elif n == 1:
        return x
    elif n == 2:
        return 2 * x ** 2 - 1
    elif n == 3:
        return 4 * x ** 3 - 3 * x
    elif n == 4:
        return 8 * x ** 4 - 8 * x ** 2
    elif n == 5:
        return 16 * x ** 5 - 20 * x ** 3 + 5 * x
    elif n == 6:
        return 32 * x ** 6 - 48 * x ** 4 + 18 * x ** 2 - 1

def chebyshev2d(x,y,a):
    A = x * 0
    ncy, ncx = a.shape
    for iy in range(ncy):
        for ix in range(ncx):
            if a[iy][ix]:
                A = A + a[iy][ix] * chebyshev(x,ix) * chebyshev(y,iy)
    return A

def crossprod(a, b):
    """CROSS PRODUCT (PROBABLY DEFINED IN SOME BUILT-IN MODULE!)"""
    return a[0] * b[1] - a[1] * b[0]

def dotprod(a, b):
    """DOT PRODUCT (PROBABLY DEFINED IN SOME BUILT-IN MODULE!)"""
    return a[0] * b[0] + a[0] * b[0]

def triarea(x, y, dir=0):
    """RETURNS THE AREA OF A TRIANGLE GIVEN THE COORDINATES OF ITS VERTICES
    A = 0.5 * | u X v |
    where u & v are vectors pointing from one vertex to the other two
    and X is the cross-product
    The dir flag lets you retain the sign (can tell if triangle is flipped)"""
    ux = x[1] - x[0]
    vx = x[2] - x[0]
    uy = y[1] - y[0]
    vy = y[2] - y[0]
    A = 0.5 * (ux * vy - uy * vx)
    if not dir:
        A = abs(A)
    return A

def CCWsort(x, y):
    """FOR A CONVEX SET OF POINTS, 
    SORT THEM SUCH THAT THEY GO AROUND IN ORDER CCW FROM THE x-AXIS"""
    xc = mean(x)
    yc = mean(y)
    ang = atanxy(x-xc, y-yc)
    SI = array(argsort(ang))
    x2 = x.take(SI, 0)
    y2 = y.take(SI, 0)
    return x2, y2

def polyarea(x, y):
    """RETURNS THE AREA OF A CONVEX POLYGON 
    GIVEN ITS COORDINATES (IN ANY ORDER)"""
    A = 0.
    x, y = CCWsort(x, y)
    for i in range(1, len(x)-1):
        xtri = x.take((0, i, i+1), 0)
        ytri = y.take((0, i, i+1), 0)
        A += triarea(xtri, ytri)
    return A

def odd(n):
    """RETURNS WHETHER AN INTEGER IS ODD"""
    return n & 1

def even(n):
    """RETURNS WHETHER AN INTEGER IS EVEN"""
    return 1 - odd(n)

def fpart(x):
    """FRACTIONAL PART OF A REAL NUMBER"""
    if type(x) in [array, list]:
        if len(x) == 1:
            x = x[0]
    return math.modf(x)[0]

def sigrange(x, nsigma=1):
    lo = percentile(gausst(nsigma), x)
    hi = percentile(gaussp(nsigma), x)
    return lo, hi

def sqrtsafe(x):
    """sqrt(x) OR 0 IF x < 0"""
    x = clip2(x, 0, None)
    return sqrt(x)

def sgn(a):
    return where(a, where(greater(a, 0), 1, -1), 0)

def sym8(a):
    """OKAY, SO THIS ISN'T QUITE RADIAL SYMMETRY..."""
    x = a + flipud(a) + fliplr(a) + transpose(a) + rot90(transpose(a),2) + rot90(a,1) + rot90(a,2) + rot90(a,3)
    return x / 8.

#def divsafe(a, b, inf=1e30, nan=0.):
def divsafe(a, b, inf=Inf, nan=NaN):
    """a / b with a / 0 = inf and 0 / 0 = nan"""
    a = array(a).astype(float)
    b = array(b).astype(float)
    asgn = greater_equal(a, 0) * 2 - 1.
    bsgn = greater_equal(b, 0) * 2 - 1.
    xsgn = asgn * bsgn
    sgn = where(b, xsgn, asgn)
    sgn = where(a, xsgn, bsgn)
    babs = clip(abs(b), 1e-200, 1e9999)
    bb = bsgn * babs
    #return where(b, a / bb, where(a, Inf, NaN))
    return where(b, a / bb, where(a, sgn*inf, nan))

def expsafe(x):
    x = array(x)
    y = []
    for xx in x:
        if xx > 708:
            y.append(1e333) # inf
        elif xx < -740:
            y.append(0)
        else:
            y.append(exp(xx))
    if len(y) == 1:
        return y[0]
    else:
        return array(y)

def floorint(x):
    return(int(floor(x)))

def ceilint(x):
    return(int(ceil(x)))

def roundint(x):
    if singlevalue(x):
        return(int(round(x)))
    else:
        return asarray(x).round().astype(int)

intround = roundint

def singlevalue(x):
    """IS x A SINGLE VALUE?  (AS OPPOSED TO AN ARRAY OR LIST)"""
    return type(x) in [type(None), float, float32, float64, int, int0, int8, int16, int32, int64]  # THERE ARE MORE TYPECODES IN Numpy

def roundn(x, ndec=0):
    if singlevalue(x):
        fac = 10.**ndec
        return roundint(x * fac) / fac
    else:
        rr = []
        for xx in x:
            rr.append(roundn(xx, ndec))
        return array(rr)

def percentile(p, x):
    x = sort(x)
    i = p * (len(x) - 1.)
    return interp(i, arange(len(x)), x)

def logical(x):
    return where(x, 1, 0)

def element_or(*l):
    """l is a list/tuple of arrays
    USAGE: x = element_or(a, b, c)"""
    x = where(l[0], l[0], l[1])
    for i in range(2,len(l)):
        x = where(x, x, l[2])
    return x

def log2(x, loexp=''):
    if loexp != '':
        x = clip2(x, 2**loexp, None)
    return log10(x) / log10(2)

def log10clip(x, loexp, hiexp=None):
    if hiexp==None:
        return log10(clip2(x, 10.**loexp, None))
    else:
        return log10(clip2(x, 10.**loexp, 10.**hiexp))

def lnclip(x, loexp):
    return log(clip2(x, e**loexp, None))

def linreg(X, Y):
    # written by William Park
    # http://www.python.org/topics/scicomp/recipes_in_python.html
    """ Returns coefficients to the regression line 'y=ax+b' from x[] and y[]. 
    Basically, it solves Sxx a + Sx b = Sxy Sx a + N b = Sy 
    where Sxy = \sum_i x_i y_i, Sx = \sum_i x_i, and Sy = \sum_i y_i. 
    The solution is a = (Sxy N - Sy Sx)/det b = (Sxx Sy - Sx Sxy)/det 
    where det = Sxx N - Sx^2. 
    In addition, 
    Var|a| = s^2 |Sxx Sx|^-1 
    = s^2 | N -Sx| / det |b| |Sx N | |-Sx Sxx| s^2
    = {\sum_i (y_i - \hat{y_i})^2 \over N-2} 
    = {\sum_i (y_i - ax_i - b)^2 \over N-2} 
    = residual / (N-2) R^2 
    = 1 - {\sum_i (y_i - \hat{y_i})^2 \over \sum_i (y_i - \mean{y})^2} 
    = 1 - residual/meanerror 
    It also prints to &lt;stdout&gt; 
    few other data, N, a, b, R^2, s^2, 
    which are useful in assessing the confidence of estimation. """
    #from math import sqrt
    if len(X) != len(Y):
        raise ValueError('unequal length')
    N = len(X)
    if N == 2: # --DC
        a = (Y[1] - Y[0]) / (X[1] - X[0])
        b = Y[0] - a * X[0]
    else:
        Sx = Sy = Sxx = Syy = Sxy = 0.0
        for x, y in map(None, X, Y):
            Sx = Sx + x
            Sy = Sy + y
            Sxx = Sxx + x*x
            Syy = Syy + y*y
            Sxy = Sxy + x*y
        det = Sxx * N - Sx * Sx
        a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
        meanerror = residual = 0.0
        for x, y in map(None, X, Y):
            meanerror = meanerror + (y - Sy/N)**2
            residual = residual + (y - a * x - b)**2
        RR = 1 - residual/meanerror
        ss = residual / (N-2)
        Var_a, Var_b = ss * N / det, ss * Sxx / det
    print("y=ax+b")
    print("N= %d" % N)
    if N == 2:
        print("a= ", a)
        print("b= ", b)
    else:
        print("a= %g \\pm t_{%d;\\alpha/2} %g" % (a, N-2, sqrt(Var_a)))
        print("b= %g \\pm t_{%d;\\alpha/2} %g" % (b, N-2, sqrt(Var_b)))
        print("R^2= %g" % RR)
        print("s^2= %g" % ss)
    return a, b


def linregrobust(x, y):
    n = len(x)
    a, b = linreg(x, y)
    dy = y - (a * x + b)
    #s = std2(dy)
    s = std(dy)
    good = less(abs(dy), 3*s)
    x, y = compress(good, (x, y))
    ng = len(x)
    if ng < n:
        print('REMOVED %d OUTLIER(S), RECALCULATING linreg' % (n - ng))
        a, b = linreg(x, y)
    return a, b


def close(x, y, rtol=1.e-5, atol=1.e-8):
    """JUST LIKE THE Numeric FUNCTION allclose, BUT FOR SINGLE VALUES.  (WILL IT BE QUICKER?)"""
    return abs(y - x) < (atol + rtol * abs(y))

def wherein(x, vals):
    """RETURNS 1 WHERE x IS IN vals"""
    try:
        good = zeros(len(x), int)
    except:
        good = 0
    for val in vals:
        good = logical_or(good, close(x, val))
    return good

def wherenotin(x, vals):
    """RETURNS 1 WHERE x ISN'T IN vals"""
    return logical_not(wherein(x, vals))

def count(a):
    """RETURNS A DICTIONARY WITH THE NUMBER OF TIMES EACH ID OCCURS"""
    bins = norep(a)
    h = histogram(a, bins)
    d = {}
    for i in range(len(h)):
        d[bins[i]] = h[i]
    return d

def rep(a):
    """RETURNS A DICTIONARY WITH THE NUMBER OF TIMES EACH ID IS REPEATED
    1 INDICATES THE VALUE APPEARED TWICE (WAS REPEATED ONCE)"""
    a = sort(a)
    d = a[1:] - a[:-1]
    c = compress(logical_not(d), a)
    if c.any():
        bins = norep(c)
        h = histogram(c, bins)
        d = {}
        for i in range(len(h)):
            d[bins[i]] = h[i]
        return d
    else:
        return {}

def norep(a):
    """RETURNS a w/o REPETITIONS, i.e. THE MEMBERS OF a"""
    a = sort(a)
    d = a[1:] - a[:-1]
    c = compress(d, a)
    x = concatenate((c, [a[-1]]))
    return x
##     l = []
##     for x in ravel(a):
##         if x not in l:
##             l.append(x)
##     return array(l)

def norepxy(x, y, tol=1e-8):
    """REMOVES REPEATS IN (x,y) LISTS -- WITHIN tol EQUALS MATCH"""
    if type(x) == type(array([])):
        x = x.tolist()
        y = y.tolist()
    else:  # DON'T MODIFY ORIGINAL INPUT LISTS
        x = x[:]
        y = y[:]
    i = 0
    while i < len(x)-1:
        j = i + 1
        while j < len(x):
            dist = hypot(x[i] - x[j], y[i] - y[j])
            if dist < tol:
                del x[j]
                del y[j]
            else:
                j += 1
        i += 1
    return x, y


def isseq(a):
    """TELLS YOU IF a IS SEQUENTIAL, LIKE [3, 4, 5, 6]"""
    return (alltrue(a == arange(len(a)) + a[0]))

def between(lo, x, hi):  # --DC
    # RETURNS 1 WHERE lo < x < hi
    # (can also use that syntax "lo < x < hi")
    if lo in [None, '']:
        try:
            good = ones(len(x)).astype(int)
        except:
            good = 1
    else:
        good = greater(x, lo)
    if hi not in [None, '']:
        good = good * less(x, hi)
    return good

def divisible(x, n): # --DC
    return (x / float(n) - x / n) < (0.2 / n)

def ndec(x, max=3):  # --DC
    """RETURNS # OF DECIMAL PLACES IN A NUMBER"""
    for n in range(max, 0, -1):
        if round(x, n) != round(x, n-1):
            return n
    return 0  # IF ALL ELSE FAILS...  THERE'S NO DECIMALS

def qkfmt(x, max=8):
    n = ndec(x, max=max)
    if n:
        fmt = '%%.%df' % n
    else:
        fmt = '%d'
    return fmt % x

def interp(x, xdata, ydata, silent=0, extrap=0):  # NEW VERSION!
    """DETERMINES y AS LINEAR INTERPOLATION OF 2 NEAREST ydata"""
    SI = argsort(xdata)
    xdata = xdata.take(SI, 0)
    ydata = ydata.take(SI, 0)
    ii = searchsorted(xdata, x)
    if singlevalue(ii):
        ii = array([ii])
    # 0 = before all
    # len(xdata) = after all
    n = len(xdata)
    if extrap:
        i2 = clip(ii,   1, n-1)
        i1 = i2 - 1
    else:
        i2 = clip(ii,   0, n-1)
        i1 = clip(ii-1, 0, n-1)
    
    x2 = take(xdata, i2)
    x1 = take(xdata, i1)
    y2 = take(ydata, i2)
    y1 = take(ydata, i1)
    # m = (y2 - y1) / (x2 - x1)
    m = divsafe(y2 - y1, x2 - x1, nan=0)
    b = y1 - m * x1
    y = m * x + b
    if len(y) == 0:
        y = y[0]
    return y

interpn = interp

def interp1(x, xdata, ydata, silent=0):  # --DC
    """DETERMINES y AS LINEAR INTERPOLATION OF 2 NEAREST ydata"""
    SI = argsort(xdata)
    # NEW numpy's take IS ACTING FUNNY
    # NO DEFAULT AXIS, MUST BE SET EXPLICITLY TO 0
    xdata = xdata.take(SI, 0).astype(float).tolist()
    ydata = ydata.take(SI, 0).astype(float).tolist()
    if x > xdata[-1]:
        if not silent:
            print(x, 'OUT OF RANGE in interp in MLab_coe.py')
        return ydata[-1]
    elif x < xdata[0]:
        if not silent:
            print(x, 'OUT OF RANGE in interp in MLab_coe.py')
        return ydata[0]
    else:
        # i = bisect(xdata, x)  # SAME UNLESS EQUAL
        i = searchsorted(xdata, x)
        if xdata[i] == x:
            return ydata[i]
        else:
            [xlo, xhi] = xdata[i-1:i+1]
            [ylo, yhi] = ydata[i-1:i+1]
            return ((x - xlo) * yhi + (xhi - x) * ylo) / (xhi - xlo)

def interpn1(x, xdata, ydata, silent=0):  # --DC
    """DETERMINES y AS LINEAR INTERPOLATION OF 2 NEAREST ydata
    interpn TAKES AN ARRAY AS INPUT"""
    yout = []
    for x1 in x:
        yout.append(interp(x1, xdata, ydata, silent=silent))
    return array(yout)

def interp2(x, xdata, ydata):  # --DC
    """LINEAR INTERPOLATION/EXTRAPOLATION GIVEN TWO DATA POINTS"""
    m = (ydata[1] - ydata[0]) / (xdata[1] - xdata[0])
    b = ydata[1] - m * xdata[1]
    y = m * x + b
    return y

def bilin(x, y, data, datax, datay):  # --DC
    """ x, y ARE COORDS OF INTEREST
    data IS 2x2 ARRAY CONTAINING NEARBY DATA
    datax, datay CONTAINS x & y COORDS OF NEARBY DATA"""
    lavg = ( (y - datay[0]) * data[1,0] + (datay[1] - y) * data[0,0] ) / (datay[1] - datay[0])
    ravg = ( (y - datay[0]) * data[1,1] + (datay[1] - y) * data[0,1] ) / (datay[1] - datay[0])
    return ( (x - datax[0]) * ravg + (datax[1] - x) * lavg ) / (datax[1] - datax[0])

def bilin2(x, y, data):  # --DC
    """ x, y ARE COORDS OF INTEREST, IN FRAME OF data - THE ENTIRE ARRAY"""
    # SHOULD BE CHECKS FOR IF x, y ARE AT EDGE OF data
    ny, nx = data.shape
    ix = int(x)
    iy = int(y)
    if ix == nx-1:
        x -= 1e-7
        ix -= 1
    if iy == ny-1:
        y -= 1e-7
        iy -= 1
    if not ((0 <= ix < nx-1) and (0 <= iy < ny-1)):
        val = 0
    else:
        stamp = data[iy:iy+2, ix:ix+2]
        datax = [ix, ix+1]
        datay = [iy, iy+1]
        # print x, y, stamp, datax, datay
        val = bilin(x, y, stamp, datax, datay)
    return val

def rand(*args):
        """rand(d1,...,dn) returns a matrix of the given dimensions
        which is initialized to random numbers from a uniform distribution
        in the range [0,1).
        """
        return RandomArray.random(args)

def eye(N, M=None, k=0, dtype=None):
        """eye(N, M=N, k=0, dtype=None) returns a N-by-M matrix where the 
        k-th diagonal is all ones, and everything else is zeros.
        """
        if M == None: M = N
        if type(M) == type('d'): 
                typecode = M
                M = N
        m = equal(subtract.outer(arange(N), arange(M)),-k)
        return asarray(m,dtype=typecode)

def tri(N, M=None, k=0, dtype=None):
        """tri(N, M=N, k=0, dtype=None) returns a N-by-M matrix where all
        the diagonals starting from lower left corner up to the k-th are all ones.
        """
        if M == None: M = N
        if type(M) == type('d'): 
                typecode = M
                M = N
        m = greater_equal(subtract.outer(arange(N), arange(M)),-k)
        return m.astype(typecode)
        
# Matrix manipulation

def diag(v, k=0):
        """diag(v,k=0) returns the k-th diagonal if v is a matrix or
        returns a matrix with v as the k-th diagonal if v is a vector.
        """
        v = asarray(v)
        s = v.shape
        if len(s)==1:
                n = s[0]+abs(k)
                if k > 0:
                        v = concatenate((zeros(k, v.dtype.char),v))
                elif k < 0:
                        v = concatenate((v,zeros(-k, v.dtype.char)))
                return eye(n, k=k)*v
        elif len(s)==2:
                v = add.reduce(eye(s[0], s[1], k=k)*v)
                if k > 0: return v[k:]
                elif k < 0: return v[:k]
                else: return v
        else:
                raise ValueError("Input must be 1- or 2-D.")
        

def fliplr(m):
        """fliplr(m) returns a 2-D matrix m with the rows preserved and
        columns flipped in the left/right direction.  Only works with 2-D
        arrays.
        """
        m = asarray(m)
        if len(m.shape) != 2:
                raise ValueError("Input must be 2-D.")
        return m[:, ::-1]

def flipud(m):
        """flipud(m) returns a 2-D matrix with the columns preserved and
        rows flipped in the up/down direction.  Only works with 2-D arrays.
        """
        m = asarray(m)
        if len(m.shape) != 2:
                raise ValueError("Input must be 2-D.")
        return m[::-1]
        
# reshape(x, m, n) is not used, instead use reshape(x, (m, n))

def rot90(m, k=1):
        """rot90(m,k=1) returns the matrix found by rotating m by k*90 degrees
        in the counterclockwise direction.
        """
        m = asarray(m)
        if len(m.shape) != 2:
                raise ValueError("Input must be 2-D.")
        k = k % 4
        if k == 0: return m
        elif k == 1: return transpose(fliplr(m))
        elif k == 2: return fliplr(flipud(m))
        elif k == 3: return fliplr(transpose(m))

def rot180(m):
    return rot90(m, 2)

def rot270(m):
    return rot90(m, 3)

def tril(m, k=0):
        """tril(m,k=0) returns the elements on and below the k-th diagonal of
        m.  k=0 is the main diagonal, k > 0 is above and k < 0 is below the main
        diagonal.
        """
        return tri(m.shape[0], m.shape[1], k=k, dtype=m.dtype.char)*m

def triu(m, k=0):
        """triu(m,k=0) returns the elements on and above the k-th diagonal of
        m.  k=0 is the main diagonal, k > 0 is above and k < 0 is below the main
        diagonal.
        """     
        return (1-tri(m.shape[0], m.shape[1], k-1, m.dtype.char))*m 

# Data analysis

# Basic operations
def max(m):
        """max(m) returns the maximum along the first dimension of m.
        """
        return maximum.reduce(m)

def min(m):
        """min(m) returns the minimum along the first dimension of m.
        """
        return minimum.reduce(m)

# Actually from BASIS, but it fits in so naturally here...

def ptp(m):
        """ptp(m) returns the maximum - minimum along the first dimension of m.
        """
        return max(m)-min(m)

def mean1(m):
        """mean(m) returns the mean along the first dimension of m.  Note:  if m is
        an integer array, integer division will occur.
        """
        return add.reduce(m)/len(m)

def mean(m, axis=0):
        """mean(m) returns the mean along the first dimension of m.  Note:  if m is
        an integer array, integer division will occur.
        """
        return add.reduce(m, axis=axis) / m.shape[axis]

def meangeom(m):
    return product(m) ** (1. / len(m))

# sort is done in C but is done row-wise rather than column-wise
def msort(m):
        """msort(m) returns a sort along the first dimension of m as in MATLAB.
        """
        return transpose(sort(transpose(m)))

def median(m):
        """median(m) returns the median of m along the first dimension of m.
        """
        if m.shape[0] & 1:
            return msort(m)[m.shape[0]/2]  # ODD # OF ELEMENTS
        else:
            return (msort(m)[m.shape[0]/2] + msort(m)[m.shape[0]/2-1]) / 2.0  # EVEN # OF ELEMENTS
            

def rms(m):
    """Root-Mean-Squared, as advertised.
    std (below) first subtracts by the mean
    and later divides by N-1 instead of N"""
    return sqrt(mean(m**2))

def std(m):
        """std(m) returns the standard deviation along the first
        dimension of m.  The result is unbiased meaning division by len(m)-1.
        """
        mu = mean(m)
        return sqrt(add.reduce(pow(m-mu,2)))/sqrt(len(m)-1.0)

stddev = std

def meanstd(m):
        """meanstd(m) returns the mean and uncertainty = std / sqrt(N-1)
        """
        mu = mean(m)
        dmu = sqrt(add.reduce(pow(m-mu,2)))/(len(m)-1.0)
        return mu, dmu

def avgstd2(m): # --DC
        """avgstd2(m) returns the average & standard deviation along the first dimension of m.
        avgstd2 ELIMINATES OUTLIERS
        The result is unbiased meaning division by len(m)-1.
        """
        done = ''
        while not done:
            n = len(m)
            mu = mean(m)
            sig = sqrt(add.reduce(pow(m-mu,2)))/sqrt(n-1.0)
            good = greater(m, mu-3*sig) * less(m, mu+3*sig)
            m = compress(good, m)
            done = sum(good) == n
            
        return [mu, sqrt(add.reduce(pow(m-mu,2)))/sqrt(len(m)-1.0)]

def std2(m): # --DC
        """std2(m) returns the standard deviation along the first dimension of m.
        std2 ELIMINATES OUTLIERS
        The result is unbiased meaning division by len(m)-1.
        """
        [a, s] = avgstd2(m)
        return s

stddev = std

def weightedavg(x, w):
    return sum(x * w) / sum(w)

weightedmean = weightedavg

## def thetaavgstd1(theta):
##     """SHWAG VERSION: WON'T WORK IF THETA SPANS A RANGE > pi
##     CALCULATES THE AVERAGE & STANDARD DEVIATION IN A LIST (OR 1-D ARRAY) OF THETA (ANGLE) MEASUREMENTS
##     RETURNS THE LIST [avg, std]    
##     NEED A NEW CODE TO HANDLE THAT: ?INCREASING WEIGHTED AVERAGES (2 POINTS AT A TIME)?"""
##     if len(theta) == 1:
##      return([theta[0], 999])
##     else:
##      # PUT ALL theta IN [0, 2 * pi]
##      for i in range(len(theta)):
##          if theta[i] < 0:
##              theta[i] = theta[i] + 2 * pi
##      if max(theta) - min(theta) > pi:
##          # "PUT ALL THETA IN [-pi, pi]"
##          for i in range(len(theta)):
##              if theta[i] > pi:
##                  theta[i] = theta[i] - 2 * pi
##      #print theta
##      if max(theta) - min(theta) > pi:
##          print "THETA RANGE TOO BIG FOR thetaavg"
##          return([999, 999])
##         else:
##          thavg = mean(theta)
##          thstd = sqrt( sum( (theta - thavg) ** 2 ) / (len(theta) - 1.) )
##          return([thavg, thstd])

def thetaavgstd(theta):
    """CALCULATES THE AVERAGE & STANDARD DEVIATION IN A LIST (OR 1-D ARRAY) OF THETA (ANGLE) MEASUREMENTS
    RETURNS THE LIST [avg, std]
    CAN HANDLE ANY RANGE OF theta
    USES INCREASING WEIGHTED AVERAGES (2 POINTS AT A TIME)"""
    n = len(theta)
    if n == 1:
        return([theta[0], 999])
    else:
        thavg = theta[0]
        for i in range(1,n):
            th = theta[i]
            if thavg - th > pi:
                thavg = thavg - 2 * pi
            elif th - thavg > pi:
                th = th - 2 * pi
            thavg = ( i * thavg + th ) / (i+1)
        for i in range(n):
            if theta[i] > thavg + pi:
                theta[i] = theta[i] - 2 * pi
        thstd = std(theta)
        return([thavg, thstd])
                


def clip2(m, m_min=None, m_max=None):
    if m_min == None:
        m_min = min(m)
    if m_max == None:
        m_max = max(m)
    return clip(m, m_min, m_max)


## def sum(m):
##      """sum(m) returns the sum of the elements along the first
##      dimension of m.
##      """
##      return add.reduce(m)
sum = add.reduce  # ALLOWS FOR AXIS TO BE INPUT --DC

def total(m):
    """RETURNS THE TOTAL OF THE ENTIRE ARRAY --DC"""
##     t = m
##     while not(type(t) in [type(1), type(1.)]):
##      t = sum(t)
##     return t
    return sum(ravel(m))

def size(m):
    """RETURNS THE TOTAL SIZE OF THE ARRAY --DC"""
    s = m.shape
    x = 1
    for n in s:
        x = x * n
    return x

def cumsum(m, axis=0):
        """cumsum(m) returns the cumulative sum of the elements along the
        first dimension of m.
        """
        return add.accumulate(m, axis=axis)

def prod(m):
        """prod(m) returns the product of the elements along the first
        dimension of m.
        """
        return multiply.reduce(m)

def cumprod(m):
        """cumprod(m) returns the cumulative product of the elments along the
        first dimension of m.
        """
        return multiply.accumulate(m)

def trapz(y, x=None):
        """trapz(y,x=None) integrates y = f(x) using the trapezoidal rule.
        """
        if x == None: d = 1
        else: d = diff(x)
        return sum(d * (y[1:]+y[0:-1])/2.0)

def cumtrapz(y, x=None, axis=0):
        """trapz(y,x=None) integrates y = f(x) using the trapezoidal rule. --DC"""
        if x == None: d = 1
        else: d = diff(x)
        if axis == 0:
            return cumsum(d * (y[1:]+y[0:-1])/2.0)
        elif axis == 1:
            return cumsum(d * (y[:,1:]+y[:,0:-1])/2.0, axis=1)
        else:
            print('YOUR VALUE OF axis = %d IS NO GOOD IN MLab_coe.cumtrapz' % axis)

def xbins(x):
    """[-0.5, 0.5, 1] --> [-1, 0, 0.75, 1.25]"""
    d = shorten(x)
    da = x[1] - x[0]
    db = x[-1] - x[-2]
    d = concatenate(([x[0] - da/2.], d, [x[-1] + db/2.]))
    return d

def diff(x, n=1):
        """diff(x,n=1) calculates the first-order, discrete difference
        approximation to the derivative.
        """
        if n > 1:
            return diff(x[1:]-x[:-1], n-1)
        else:
            return x[1:]-x[:-1]

def shorten(x, n=1):
        """shorten(x,n=1) 
        SHORTENS x, TAKING AVG OF NEIGHBORS, RECURSIVELY IF n > 1
        """
        a = (x[1:] + x[:-1]) / 2.
        if n > 1:
            return avg(a, n-1)
        else:
            return a

def lengthen(x, n):
    """lengthen([0, 1, 5], 4) ==> 0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5"""
    x = array(x)
    d = diff(x)
    i = arange(n) / float(n)
    o = outer(i, d)
    o = o + x[:-1]
    o = ravel(transpose(o))
    o = concatenate((o, [x[-1]]))
    return o

def powerlaw(x, y):
    """RETURNS EXPONENT n TO POWER LAW FIT y ~ x^n
    AT POINTS ON AVERAGED x"""
    logx = log10(x)
    logy = log10(y)
    
    dlogx = diff(logx)
    dlogy = diff(logy)
    
    dd = dlogy / dlogx
    x2 = (x[1:] + x[:-1]) / 2
    
    return x2, dd


def grad(m):
    """Calculates the gradient of the matrix m using the finite difference method
    The result will be 2 arrays, one for each of the axes x & y, respectively,
    with each having dimension (N-2, N-2), where m was (N, N).
    The coordinates will be in between of those of m.  --DC"""
    ay = (m[2:]   - m[:-2]) / 2.       # (N-2, N)
    ax = (m[:,2:] - m[:,:-2]) / 2.     # (N,   N-2)
    ay = ay[:,1:-1]                    # (N-2, N-2)
    ax = ax[1:-1,:]
    return array([ax, ay])

def laplacian(m):
    """Calculates the laplacian of the matrix m
    using the finite differencing method.
    The result will have dimension (ny-2, nx-2) where m had (ny, nx).
    see Fig. 2 of Bradac & Schneider 2005
    (Strong & Weak Lensing United I)
    although theirs is a factor of 1/2 too low.
    """
    ny, nx = m.shape
    center = m[1:-1,1:-1]
    
    sides = zeros(center.shape, float)
    for dx,dy in [(-1,0), (0,1), (1,0), (0,-1)]:
        sides = sides + m[1+dy:ny-1+dy, 1+dx:nx-1+dx]
    
    corners = zeros(center.shape, float)
    for dx,dy in [(-1,-1), (-1,1), (1,1), (1,-1)]:
        corners = corners + m[1+dy:ny-1+dy, 1+dx:nx-1+dx]
    
    return (2*corners - sides - 4*center) / 3.

def corrcoef(x, y=None):
        """The correlation coefficients
        """
        c = cov(x, y)
        d = diag(c)
        return c/sqrt(multiply.outer(d,d))

def cov(m,y=None):
        m = asarray(m)
        mu = mean(m)
        if y != None: m = concatenate((m,y))
        sum_cov = 0.0
        for v in m:
                sum_cov = sum_cov+multiply.outer(v,v)
        return (sum_cov-len(m)*multiply.outer(mu,mu))/(len(m)-1.0)

# Added functions supplied by Travis Oliphant
#import numpy.linalg.old as LinearAlgebra
def squeeze(a):
    "squeeze(a) removes any ones from the shape of a"
    b = asarray(a.shape)
    reshape (a, tuple (compress (not_equal (b, 1), b)))
    return

def kaiser(M,beta):
    """kaiser(M, beta) returns a Kaiser window of length M with shape parameter
    beta. It depends on the cephes module for the modified bessel function i0.
    """
    import cephes
    n = arange(0,M)
    alpha = (M-1)/2.0
    return cephes.i0(beta * sqrt(1-((n-alpha)/alpha)**2.0))/cephes.i0(beta)

def blackman(M):
    """blackman(M) returns the M-point Blackman window.
    """
    n = arange(0,M)
    return 0.42-0.5*cos(2.0*pi*n/M) + 0.08*cos(4.0*pi*n/M)


def bartlett(M):
    """bartlett(M) returns the M-point Bartlett window.
    """
    n = arange(0,M)
    return where(less_equal(n,M/2.0),2.0*n/M,2.0-2.0*n/M)

def hanning(M):
    """hanning(M) returns the M-point Hanning window.
    """
    n = arange(0,M)
    return 0.5-0.5*cos(2.0*pi*n/M)

def hamming(M):
    """hamming(M) returns the M-point Hamming window.
    """
    n = arange(0,M)
    return 0.54-0.46*cos(2.0*pi*n/M)

def sinc(x):
    """sinc(x) returns sin(pi*x)/(pi*x) at all points of array x.
    """
    return where(equal(x,0.0),1.0,sin(pi*x)/(pi*x))

#def eig(v):
#    """[x,v] = eig(m) returns the the eigenvalues of m in x and the corresponding
#    eigenvectors in the rows of v.
#    """
#    return LinearAlgebra.eigenvectors(v)

#def svd(v):
#    """[u,x,v] = svd(m) return the singular value decomposition of m.
#    """
#    return LinearAlgebra.singular_value_decomposition(v)


def histogram(a, bins):
    n = searchsorted(sort(a), bins)
    n = concatenate([n, [len(a)]])
    return n[1:]-n[:-1]

def cumhisto(a,da=1.,amin=[],amax=[]): # --DC
    """
    Histogram of 'a' defined on the bin grid 'bins'
       Usage: h=histogram(p,xp)
    """
    if amin == []:
        amin = min(a)
    if amax == []:
        amax = max(a)
    nnn = (amax - amin) / da
    if less(nnn - int(nnn), 1e-4):
        amax = amax + da
    bins = arange(amin,amax+da,da)
    n=searchsorted(sort(a),bins)
    n=array(list(map(float,n)))
    return n[1:]

def cumHisto(a,da=1.,amin=[],amax=[]): # --DC
    if amin == []:
        amin = min(a)
    if amax == []:
        amax = max(a)
    h = cumhisto(a, da, amin, amax)
    return Histogram(h, amin, da)

def plotcumhisto(a,da=1.,amin=[],amax=[]): # --DC
    p = FramedPlot()
    p.add(cumHisto(a, da, amin, amax))
    p.show()
    return p

# from useful_coe.py
def histo(a,da=1.,amin=[],amax=[]): # --DC
    """
    Histogram of 'a' defined on the bin grid 'bins'
       Usage: h=histogram(p,xp)
    """
    if amin == []:
        amin = min(a)
    if amax == []:
        amax = max(a)
    nnn = (amax - amin) / da
    if less(nnn - int(nnn), 1e-4):
        amax = amax + da
    bins = arange(amin,amax+da,da)
    n=searchsorted(sort(a),bins)
#    n=concatenate([n,[len(a)]])
    n=array(list(map(float,n)))
##     print a
##     print bins
##     print n
    return n[1:]-n[:-1]
#    return hist(a, bins)

def Histo(a,da=1.,amin=[],amax=[], **other): # --DC
    if amin == []:
        amin = min(a)
    if amax == []:
        amax = max(a)
    try:
        amin = amin[0]
    except:
        pass
##     print 'hi'
##     print da
##     print amin
##     print amax
    h = histo(a, da, amin, amax)
##     print h
    return Histogram(h, amin, da, **other)

def plothisto(a,da=1.,amin=[],amax=[]): # --DC
    p = FramedPlot()
    p.add(Histo(a, da, amin, amax))
    p.show()

def bargraphbiggles(x, y, fill=1, color='black', **other):
    n = len(x)
    xx = repeat(x, 2)
    y = y.astype(float)
    z = array([0.])
    yy = concatenate([z, repeat(y, 2), z])
    zz = yy*0
    
    p = FramedPlot()
    if fill:
        p.add(FillBetween(xx, yy, xx, zz, color=color))
    else:
        p.add(Curve(xx, yy, color=color, **other))
    p.show()

def BarGraph(x, y, fill=1, color='black', bottom=0, **other):
    n = len(x)
    xx = repeat(x, 2)
    y = y.astype(float)
    z = array([0.])
    yy = concatenate([z, repeat(y, 2), z])
    zz = yy*0 + bottom
    if fill:
        return FillBetween(xx, yy, xx, zz, color=color)
    else:
        return Curve(xx, yy, color=color, **other)

def histob(a,da=1.,amin=[],amax=[]): # --DC
    # NOTE searchsorted can't be counted on to act consistently
    #   when bin values are equal to data values
    # for example, neither 0.04 or 0.05 gets put in the 0.04-0.05 bin
    #   0.04 gets put in the bin below, but 0.05 gets put in the bin above
    # So it's good to stagger your bins values when necessary (0.035, 0.045, 0.055)
    """
    Histogram of 'a' defined on the bin grid 'bins'
       Usage: h=histogram(p,xp)
    """
    if amin == []:
        amin = min(a)
    if amax == []:
        amax = max(a)
    # MAKE SURE 18 GOES IN THE 18-18.9999 bin (for da=1 anyway)
    amin = amin - 1e-4
    amax = amax + 1e-4
    #if less(abs(amax - a[-1]), da*1e-4):
    nnn = (amax - amin) / da
    if less(nnn - int(nnn), 1e-4):
        amax = amax + da
    #bins = arange(amin,amax+da,da)
    bins = arange(amin,amax+da,da)
    n=searchsorted(sort(a),bins)
    n=array(list(map(float,n)))
    n = n[1:]-n[:-1]
    return (bins, n)

def Histob(a, da=1., amin=[], amax=[], fill=1, color='black', bottom=0):
    bins, n = histob(a, da, amin, amax)
    return BarGraph(bins, n, fill=fill, color=color, bottom=bottom)

def histov(a, bins, v, presorted=0):
    """Total of values (v) in bins
    (other historgrams just count number of elements in bins)"""
    if not presorted:
        SI = argsort(a)
        a = take(a, SI)
        v = take(v, SI)
    vcum = cumsum(v)
    i = searchsorted(a, bins)
    i = i[1:] - 1
    vcumi = vcum.take(i)
    vcumi = concatenate([[0], vcumi])
    vb = vcumi[1:] - vcumi[:-1]
    return vb

#def isNaN(x):
#    return (x == 1) and (x == 0)

def isNaN(x):
    return not (x < 0) and not (x > 0) and (x != 0)

def isnan(x):
    l = less(x, 0)
    g = greater(x, 0)
    e = equal(x, 0)
    n = logical_and(logical_not(l), logical_not(g))
    n = logical_and(n, logical_not(e))
    return n
    
#from coeplot2a import *
#testinsidepoly()




pwd = os.getcwd
die = sys.exit

def color1to255(color):
    return tuple((array(color) * 255. + 0.49).astype(int).tolist())  # CONVERT TO 0-255 SCALE

def color255to1(color):
    return tuple((array(color) / 255.).tolist())  # CONVERT TO 0-255 SCALE

def color2hex(color):
    if 0:  # 0 < max(color) <= 1:  # 0-1 SCALE
        # BUT EVERY ONCE IN A WHILE, YOU'LL GET A (0,0,1) OUT OF 255...
        color = color1to255(color)
    colorhex = '#'
    for val in color:
        h = hex(val)[2:]
        if len(h) == 1:
            h = '0'+h
        colorhex += h
    return colorhex

###

def keyvals(k, keys, vals):
    """GIVEN {keys: vals}, RETURNS VALUES FOR k
    THERE MUST BE A BUILT-IN WAY OF DOING THIS!"""
    d = dict(list(zip(keys, vals)))
    d[0] = 0
    f = lambda x: d[x]
    v = list(map(f, ravel(k)))
    if type(k) == type(array([])):
        v = array(v)
        v.shape = k.shape
    return v

def printmult(x, n):
    if not (x % n):
        print(x)

def cd(dir):
    if len(dir) > 2:
        if dir[0:2] == '~/':
            dir = os.path.join(home, dir[2:])
    os.chdir(dir)

def cdmk(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
    os.chdir(dir)

def splitparagraphs(txt):
    paragraphs = ['']
    for line in txt:
        line = line.strip()
        if not line:
            line = '\n'
        if line[-1] != '\n':
            line += '\n'
        if line == '\n':
            paragraphs.append('')
        else:
            #paragraphs[-1].append(line)
            paragraphs[-1] += line
    if paragraphs[-1] == '':
        paragraphs = paragraphs[:-1]
    return paragraphs

def echo(word):
    cmd = 'echo ' + word
    subproc = popen2.Popen4(cmd)
    out = subproc.fromchild.readlines()  # SExtractor output
    out = out[0][:-1]  # LIST OF 1 STRING WITH \n AT END
    return out

home = os.environ.get('HOME', '')

def singlevalue(x):
    """IS x A SINGLE VALUE?  (AS OPPOSED TO AN ARRAY OR LIST)"""
    # return type(x) in [float, int]  THERE ARE MORE TYPECODES IN Numpy
    return type(x) in [float, float32, float64, int, int0, int8, int16, int32, int64]  # THERE ARE MORE TYPECODES IN Numpy
##     try:
##         a = x[0]
##         singleval = False
##     except:
##         singleval = True
##     return singleval

def comma(x, ndec=0):
    if ndec:
        format = '%%.%df' % ndec
        s = format % x
        si, sf = s.split('.')
        sf = '.' + sf
    else:
        s = '%d' % x
        si = s
        sf = ''
    ss = ''
    while len(si) > 3:
        ss = ',' + si[-3:] + ss
        si = si[:-3]
    ss = si + ss + sf
    return ss

# print comma(9812345.67)

def th(n):
    """RETURNS 0th, 1st, 2nd, 3rd, 4th, 5th, etc."""
    if n == 1:
        return '1st'
    elif n == 2:
        return '2nd'
    elif n == 3:
        return '3rd'
    else:
        return '%dth' % n

nth = th

def num2str(x, max=3):
    try:
        n = ndec(x, max)
        if n:
            return "%%.%df" % n % x
        else:
            return "%d" % x
    except:
        return x

def str2num(stri, rf=0):
    """CONVERTS A STRING TO A NUMBER (INT OR FLOAT) IF POSSIBLE
    ALSO RETURNS FORMAT IF rf=1"""
    try:
        num = stri.atoi()
        format = 'd'
    except:
        try:
            num = stri.atof()
            format = 'f'
        except:
            if not stri.strip():
                num = None
                format = ''
            else:
                num = stri
                format = 's'
    if rf:
        return (num, format)
    else:
        return num

def minmax(x, range=None):
    if range:
        lo, hi = range
        good = between(lo, x, hi)
        x = compress(good, x)
    return min(x), max(x)

#############################################################################
# ARRAYS
#
# PYTHON USES BACKWARDS NOTATION: a[row,column] OR a[iy,ix] OR a[iy][ix]
# NEED TO MAKE size GLOBAL (I THINK) OTHERWISE, YOU CAN'T CHANGE IT!
# COULD HAVE ALSO USED get_data IN ~txitxo/Python/useful.py

def FltArr(n0,n1):
    """MAKES A 2-D FLOAT ARRAY"""
    #a = ones([n0,n1], dtype=float32)
    a = ones([n0,n1], float32)
    return(a[:])


def IndArr(n0,n1):
    """MAKES A 2-D INTEGER ARRAY WITH INCREASING INDEX"""
    a = arange(n0*n1)
    return resize(a, [n0,n1])

#################################
# STRINGS, INPUT

def striskey(stri):
    """IS stri AN OPTION LIKE -C or -ker
    (IT'S NOT IF IT'S -2 or -.9)"""
    iskey = 0
    if stri:
        if stri[0] == '-':
            iskey = 1
            if len(stri) > 1:
                iskey = stri[1] not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.']
    return iskey
    

def pause(text=''):
    inp = input(text)

def wait(seconds):
    t0 = time()
    t1 = time()
    while (t1 - t0) < seconds:
        t1 = time()

def inputnum(question = ''):
    done = 0
    while not done:
        rinp = input(question)
        try: 
            x = rinp.atof()
            done = 1
        except: 
            pass
    try: x = rinp.atoi()
    except: pass
    return x

def stringsplitatoi(stri, separator=''):
    if separator:
        words = stri.split(separator)
    else:
        words = stri.split()
    vals = []
    for word in words:
        vals.append(word.atoi())
    return vals

def stringsplitatof(stri, separator=''):
    if separator:
        words = stri.split(separator)
    else:
        words = stri.split()
    vals = []
    for word in words:
        vals.append(word.atof())
    return vals

def stringsplitstrip(stri, separator=''):
    # SPLITS BUT ALSO STRIPS EACH ITEM OF WHITESPACE
    if separator:
        words = stri.split( separator)
    else:
        words = stri.split()
    vals = []
    for word in words:
        vals.append(word.strip())
    return vals

def strbegin(stri, phr):
    return stri[:len(phr)] == phr

def strend(stri, phr):
    return stri[-len(phr):] == phr

def strfindall(stri, phr):
    """FIND ALL INSTANCES OF phr IN stri
    RETURN LIST OF POSITIONS WHERE phr IS FOUND
    (OR RETURN [] IF NOT FOUND)"""
    pos = []
    start = -1
    while 1:
        start = stri.find(phr, start+1)
        if start > -1:
            pos.append(start)
        else:
            break
    return pos

def strbtw1(s, left, right=None):
    """RETURNS THE PART OF STRING s BETWEEN left & right
    EXAMPLE strbtw('det_lab.reg', '_', '.') RETURNS 'lab'
    EXAMPLE strbtw('det_{a}.reg', '{}') RETURNS 'a'"""
    out = None
    if right == None:
        if len(left) == 1:
            right = left
        elif len(left) == 2:
            left, right = left
    i1 = s.find(left)
    if (i1 > -1):
        i1 += len(left) - 1
        i2 = s.find(right, i1+1)
        if (i2 > i1):
            out = s[i1+1:i2]
    #out = string.split(s, left)[1]
    #out = string.split(out, right)[0]
    return out

def strbtw(s, left, right=None, r=False):
    """RETURNS THE PART OF STRING s BETWEEN left & right
    EXAMPLE strbtw('det_lab.reg', '_', '.') RETURNS 'lab'
    EXAMPLE strbtw('det_{a}.reg', '{}') RETURNS 'a'
    EXAMPLE strbtw('det_{{a}, b}.reg', '{}', r=1) RETURNS '{a}, b'"""
    out = None
    if right == None:
        if len(left) == 1:
            right = left
        elif len(left) == 2:
            left, right = left
    i1 = s.find(left)
    if (i1 > -1):
        i1 += len(left) - 1
        if r:  # search from the right
            i2 = s.rfind(right, i1+1)
        else:
            i2 = s.find(right, i1+1)
        if (i2 > i1):
            out = s[i1+1:i2]
    #out = string.split(s, left)[1]
    #out = string.split(out, right)[0]
    return out

def getanswer(question=''):
    ans = -1
    while ans == -1:
        inp = input(question)
        if inp:
            if inp[0].upper() == 'Y':
                ans = 1
            if inp[0].upper() == 'N':
                ans = 0
    return ans

ask = getanswer

#################################
# LISTS

def putids(selfvalues, selfids, ids, values):
    """ selfvalues = INITIAL ARRAY -OR- A DEFAULT VALUE FOR UNput ELEMENTS """
    try:
        n = len(selfvalues)
    except:
        n = len(selfids)
        selfvalues = zeros(n, int) + selfvalues
    indexlist = zeros(max(selfids)+1, int) - 1
    put(indexlist, array(selfids).astype(int), arange(len(selfids)))
    indices = take(indexlist, array(ids).astype(int))
    put(selfvalues, indices, values)
    return selfvalues

def takelist(a, ind):
    l = []
    for i in ind:
        l.append(a[i])
    return l

def common(id1, id2):
    # ASSUME NO IDS ARE NEGATIVE
    id1 = array(id1).astype(int)
    id2 = array(id2).astype(int)
    n = max((max(id1), max(id2)))
    in1 = zeros(n+1, int)
    in2 = zeros(n+1, int)
    put(in1, id1, 1)
    put(in2, id2, 1)
    inboth = in1 * in2
    ids = arange(n+1)
    ids = compress(inboth, ids)
    return ids

# FROM sparse.py ("sparse3")
def census(a, returndict=1):
    a = sort(ravel(a))
    if returndict:
        i = arange(min(a), max(a)+2)
    else:
        i = arange(max(a)+2)
    s = searchsorted(a, i)
    s = s[1:] - s[:-1]
    i = i[:-1]
    if returndict:
        print(i)
        print(s)
        #i, s = compress(s, (i, s))
        i = compress(s, i)
        s = compress(s, s)
        print('is')
        print(i)
        print(s)
        d = {}
        for ii in range(len(i)):
            d[i[ii]] = s[ii]
        return d
    else:
        return s

# ALSO CONSIDER: set(all) - set(ids)
def invertselection(ids, all):
    if type(all) == int:  # size input
        all = arange(all) + 1
        put(all, array(ids)-1, 0)
        all = compress(all, all)
        return all
    else:
        out = []
        for val in all:
            #if val not in ids:
            if not floatin(val, ids):
                out.append(val)
        return out

def mergeids(id1, id2):
    # ASSUME NO IDS ARE NEGATIVE
    id1 = array(id1).astype(int)
    id2 = array(id2).astype(int)
    idc = common(id1, id2)
    id3 = invertselection(idc, id2)
    return concatenate((id1, id3))


def findmatch1(x, xsearch, tol=1e-4):
    """RETURNS THE INDEX OF x WHERE xsearch IS FOUND"""
    i = argmin(abs(x - xsearch))
    if abs(x[i] - xsearch) > tol:
        print(xsearch, 'NOT FOUND IN findmatch1')
        return -1
    else:
        return i

def findmatch(x, y, xsearch, ysearch, dtol=4, silent=0, returndist=0, xsorted=0):
    """FINDS AN OBJECT GIVEN A LIST OF POSITIONS AND SEARCH COORDINATE
    RETURNS INDEX OF THE OBJECT OR n IF NOT FOUND"""
    
    n = len(x)
    if silent < 0:
        print('n=', n)
    if not xsorted:
        SI = argsort(x)
        x = take(x, SI)
        y = take(y, SI)
    else:
        SI = arange(n)
    
    dist = 99  # IN CASE NO MATCH IS FOUND
    
    # SKIP AHEAD IN CATALOG TO x[i] = xsearch - dtol
    #print "SEARCHING..."
    if xsearch > dtol + max(x):
        done = 'too far'
    else:
        done = ''
        i = 0
        while xsearch - x[i] > dtol:
            if silent < 0:
                print(i, xsearch, x[i])
            i = i + 1
    
    while not done:
        if silent < 0:
            print(i, x[i], xsearch)
        if x[i] - xsearch > dtol:
            done = 'too far'
        else:
            dist = sqrt( (x[i] - xsearch) ** 2 + (y[i] - ysearch) ** 2)
            if dist < dtol:
                done = 'found'
            elif i == n - 1:
                done = 'last gal'
            else:
                i = i + 1
        if silent < 0:
            print(done)
    
    if done == 'found':
        if not silent:
            print('MATCH FOUND %1.f PIXELS AWAY AT (%.1f, %.1f)' % (dist, x[i], y[i]))
        ii = SI[i]
    else:
        if not silent:
            print('MATCH NOT FOUND')
        ii = n
    if returndist:
        return ii, dist
    else:
        return ii

def findmatches2(x1, y1, x2, y2):
    """MEASURES ALL DISTANCES, FINDS MINIMA
    SEARCHES FOR 2 IN 1
    RETURNS INDICES AND DISTANCES"""
    dx = subtract.outer(x1, x2)
    dy = subtract.outer(y1, y2)
    d = sqrt(dx**2 + dy**2)
    i = argmin(d,0)
    
    n1 = len(x1)
    n2 = len(x2)
    j = arange(n2)
    di = n2*i + j
    dmin = take(d,di)
    return i, dmin


def xref(data, ids, idcol=0, notfoundval=None):
    """CROSS-REFERENCES 2 DATA COLUMNS
    data MAY EITHER BE A 2-COLUMN ARRAY, OR A FILENAME CONTAINING THAT DATA
    ids ARE THE KEYS -- THE VALUES CORRESPONDING TO THESE (IN data's OTHER COLUMN) ARE RETURNED
    idcol TELLS WHICH COLUMN THE ids ARE IN (0 OR 1)"""
    if type(data) == stri:
        data = transpose(loaddata(data))
    iddata = data[idcol].astype(int)
    xrefdata = data[not idcol].astype(int)

    dict = {}
    for i in range(len(iddata)):
        dict[iddata[i]] = xrefdata[i]

    xrefs = []
    for id in ids:
        xrefs.append(dict.get(id, notfoundval))

    return array(xrefs)


def takeid(data, id):
    """TAKES data COLUMNS CORRESPONDING TO id.
    data's ID's ARE IN ITS FIRST ROW"""
    dataids = data[0].astype(int)
    id = int(id)
    outdata = []
    i = 0
    while id != dataids[i]:
        i += 1
    return data[:,i]

def takeids(data, ids, idrow=0, keepzeros=0):
    """TAKES data COLUMNS CORRESPONDING TO ids.
    data's ID's ARE IN idrow, ITS FIRST ROW BY DEFAULT"""
    dataids = data[idrow].astype(int)
    ids = ids.astype(int)
    outdata = []
    n = data.shape[1]
    for id in ids:
        gotit = 0
        for i in range(n):
            if id == dataids[i]:
                gotit = 1
                break
        if gotit:
            outdata.append(data[:,i])
        elif keepzeros:
            outdata.append(0. * data[:,0])
    return transpose(array(outdata))
    

#################################
# FLUX, BPZ

bpzpath = os.environ.get('BPZPATH', '')

def bpzsedname(tb, seds, interp=2):
    if type(seds) == stri:
        seds = loadfile(bpzpath + '/SED/' + seds)
    rb = roundint(tb)
    name = seds[rb-1]
    if abs(rb - tb) > 0.1:
        rb2 = roundint((tb - rb) * 3 + rb)
        name = name[:-4] + '-' + seds[rb2-1]
    return name

def bpztypename(tb, tbs, interp=2):
    rb = roundint(tb)
    name = tbs[rb-1]
    if abs(rb - tb) > 0.1:
        rb2 = roundint((tb - rb) * 3 + rb)
        name += '-' + tbs[rb2-1]
    return name

def addmags(m1, m2, dm1=0, dm2=0):
    # F = 10 ** (-0.4 * m)
    # dF = -0.4 * ln(10) * 10 ** (-0.4 * m) * dm = -0.921034 * F * dm
    # somehow this is wrong, should be:
    # dF / F = 10 ** (0.4 * dm) - 1  (as in bpz_tools.e_frac2mag)
    if (m1 >= 99 and m2 >= 99) or (dm1 >= 99 and dm2 >= 99):
        m = 99
        dm = 99
    elif m1 >= 99 or dm1 >= 99:
        m = m2
        dm = dm2
    elif m2 >= 99 or dm2 >= 99:
        m = m1
        dm = dm1
    else:  # NORMAL SITUATION
        F1 = 10 ** (-0.4 * m1)
        F2 = 10 ** (-0.4 * m2)
        F = F1 + F2
        m = -2.5 * log10(F)
        #dF1 = 0.921034 * F1 * dm1
        #dF2 = 0.921034 * F2 * dm2
        #dF = sqrt(dF1 ** 2 + dF2 ** 2)
        #dm = dF / F / 0.921034
        dm = sqrt( (F1 * dm1) ** 2 + (F2 * dm2) ** 2 ) / F
    output = (m, dm)
    
    return output

def addfluxes(F1, F2, dF1=0, dF2=0):
    F = F1 + F2
    dF = sqrt(dF1 ** 2 + dF2 ** 2)
    output = (F, dF)
    
    return output



#################################
# FROM Txitxo's bpz_tools.py

def sex2bpzmags(f,ef,zp=0.,sn_min=1.):
    """
    This function converts a pair of flux, error flux measurements from SExtractor
    into a pair of magnitude, magnitude error which conform to BPZ input standards:
    - Nondetections are characterized as mag=99, errormag=+m_1sigma
      - corrected error in previous version: was errormag=-m_1sigma
    - Objects with absurd flux/flux error combinations or very large errors are
      characterized as mag=-99 errormag=0.
    """

    nondetected=less_equal(f,0.)*greater(ef,0) #Flux <=0, meaningful phot. error
    nonobserved=less_equal(ef,0.) #Negative errors
    #Clip the flux values to avoid overflows
    f=clip(f,1e-100,1e10)
    ef=clip(ef,1e-100,1e10)
    nonobserved+=equal(ef,1e10)
    nondetected+=less_equal(f/ef,sn_min) #Less than sn_min sigma detections: consider non-detections
    
    detected=logical_not(nondetected+nonobserved)

    m=zeros(len(f), float)
    em=zeros(len(ef), float)

    m = where(detected,-2.5*log10(f)+zp,m)
    m = where(nondetected,99.,m)
    m = where(nonobserved,-99.,m)

    em = where(detected,2.5*log10(1.+ef/f),em)
    #em = where(nondetected,2.5*log10(ef)-zp,em)
    em = where(nondetected,zp-2.5*log10(ef),em)
    #print "NOW WITH CORRECT SIGN FOR em"
    em = where(nonobserved,0.,em)
    return m,em


# NOTE PLACEMENT OF THIS LINE IS IMPORTANT
# coeio ALSO IMPORTS FROM coetools (THIS MODULE)
# SO TO AVOID AN INFINITE LOOP, coeio ONLY LOADS FROM coetools
#  THOSE FUNCTIONS DEFINED BEFORE coeio IS LOADED
#from smooth import *

