#!/usr/bin/env python
import numpy
import os
import plottools
import pylab
from astropy.io import fits
from glob import glob

from matplotlib import rcParams
for tick in ('xtick', 'ytick'):
    rcParams['{0}.major.size'.format(tick)] = 8
    rcParams['{0}.minor.size'.format(tick)] = 4
    rcParams['{0}.major.width'.format(tick)] = 2
    rcParams['{0}.minor.width'.format(tick)] = 2
    rcParams['{0}.labelsize'.format(tick)] = 20
rcParams['axes.linewidth'] = 2
rcParams['axes.labelsize'] = 22
rcParams['font.size'] = 28
rcParams['legend.fontsize'] = 18
rcParams['lines.linewidth'] = 2
rcParams['mathtext.fontset'] = 'stix'
rcParams['pdf.use14corefonts'] = True
rcParams['text.usetex'] = True
rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']

red = (1,0,0)
green = (0.2,0.6,0)
blue = (0,0,1)
yellow = (1,1,0.2)
magenta = (1,0.4,0.6)
cyan = (0.2,0.7,1)
orange = (1,0.5,0)
purple = (0.8,0,0.4)

"""
Make a corner plot of the three models on top of each other

"""

def main(burn=100000):
    suff = '-redshift_0.0_0.5-z_0.00_0.15-rmax_3.0.fits'
    files = ['fullnfw{0}'.format(suff),
             'fullnfw_halfduffy{0}'.format(suff),
             'fullnfw_twiceduffy{0}'.format(suff)]
    files = [os.path.join('output', f) for f in files]
    data = [fits.open(f)[1].data for f in files]
    params = [[d['Msat1'][burn:], d['chost'][burn:], d['Mhost1'][burn:]]
              for d in data]
    labels = (r'$\log\,M_{\rm sub}$', r'$c_{\rm host}$',
              r'$\log\,M_{\rm host}$')
    names = ('Duffy+', r'$1/2\times$Duffy+', r'$2\times$Duffy+')
    plottools.corner(params, labels=labels, names=names,
                     bcolor=(orange, cyan, purple),
                     medians1d=False, percentiles1d=False)
    output = 'plots/corner_all.pdf'
    pylab.savefig(output, format='pdf')
    return

main()
