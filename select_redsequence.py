#!/usr/bin/env python
import numpy
import os
import pylab
import readfile
from astro import cosmology
from astro.clusters import redsequence
from glob import glob
from itertools import count, izip
from matplotlib import cm, colors as mplcolors, ticker

# local
from pytools import savefig

# until a newer matplotlib is available
import colormaps

from matplotlib import rcParams
for tick in ('xtick', 'ytick'):
    rcParams['{0}.major.size'.format(tick)] = 8
    rcParams['{0}.minor.size'.format(tick)] = 4
    rcParams['{0}.major.width'.format(tick)] = 2
    rcParams['{0}.minor.width'.format(tick)] = 2
    rcParams['{0}.labelsize'.format(tick)] = 20
rcParams['axes.linewidth'] = 2
rcParams['axes.labelsize'] = 22
rcParams['font.size'] = 22
rcParams['legend.fontsize'] = 18
rcParams['lines.linewidth'] = 2
rcParams['mathtext.fontset'] = 'stix'
rcParams['pdf.use14corefonts'] = True
rcParams['text.usetex'] = True
rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']

red = (1,0,0)
green = (0.2,0.6,0)
blue = (0,0,1)
yellow = (1,1,0.4)
magenta = (1,0.4,0.6)
cyan = (0,1,1)
orange = (1,0.7,0)
purple = (0.8,0,0.4)


def main():
    # includes both "all" and "good"
    ls = sorted(glob('catalogs/20160716/*/*_udgs.cat'))
    for cat in ls:
        run(cat)
        return
    return


def run(catalog):
    """
    Need to correct for extinction?
    """
    cluster = os.path.split(catalog)[1].split('_')[0]
    print cluster
    hdr = readfile.header(catalog)
    udgs = readfile.dict(catalog)
    print '    {0} UDGs'.format(udgs['x'].size)
    # galaxy catalogs
    gals = os.path.join('..', '..', 'catalogs', 'photometry', 'galaxies',
                        '{0}-galaxies_1Mpc.cat'.format(cluster))
    gals = readfile.dict(gals)
    rs = os.path.join('..', '..', 'catalogs', 'photometry', 'redsequence',
                      '0.10Lstar', '{0}-rs.cat'.format(cluster))
    Ag, Ar = readfile.table(rs, include='##', cols=2)[1:3]
    print '    Ag = {0:.3f} | Ar = {1:.3f}'.format(Ag, Ar)
    udgs['mag'] -= Ar
    udgs['grcolor'] += Ag - Ar
    rs = readfile.dict(rs)

    # plot UDGs on top of the red sequence
    fig, ax = pylab.subplots(figsize=(5,4))
    ax.plot(rs['rmag'], rs['color'], 'r.')
    ax.plot(udgs['mag'], udgs['grcolor'], 'b.')
    ax.set_xlabel('r-band magnitude')
    ax.set_ylabel('g-r color')
    fig.tight_layout(pad=0.4)
    output = 'plots/redsequence/{0}.pdf'.format(cluster)
    savefig(output, fig)
    return


main()

