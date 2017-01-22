#!/usr/bin/env python
import numpy
import os
import pylab
import readfile
from astLib import astCoords, astWCS
from astro import cosmology
from glob import glob
from lnr import to_log
from matplotlib import cm, colors as mplcolors, ticker
from numpy import log10

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
    dyn = readfile.dict('../../dynamics/output/dynamics.out')
    catalogs = sorted(glob('catalogs/20160716/good/*_udgs.cat'))
    # keys: ['rating', 'logmstar', 'cluster', 'sersicn', 'mag',
           # 'y', 'x', 'Rsat', 'grcolor', 'id', 'meansurfb', 'reff']
    ncl = len(catalogs)
    fig, ax = pylab.subplots(figsize=(5.5,5))
    clusters = []
    texlines = []
    z, m200, m200_err, r200, r200_err = numpy.zeros((5,ncl))
    nudg = numpy.zeros(ncl, dtype=int)
    # average reff per cluster
    size, size_lo, size_hi = numpy.zeros((3,ncl))
    for i, cat in enumerate(catalogs):
        cluster = os.path.split(cat)[1].split('_')[0]
        clusters.append(cluster)
        j = (dyn['cluster'] == cluster)
        z[i] = dyn['z'][j]
        m200[i] = dyn['m200'][j]
        m200_err[i] = dyn['m200_err'][j]
        r200[i] = dyn['r200'][j]
        r200_err[i] = dyn['r200_err'][j]
        texline = r'{0:<10s} & {1:.3f} & {2:>4.1f} & {3:.1f}'.format(
                        clusters[i], z[i], m200[i], m200_err[i])
        texline = r'{0} & ${1:.1f}\pm{2:.1f}$'.format(
                        texline, r200[i], r200_err[i])
        full = cat.replace('good/', 'all/').replace('_good', '')
        good = readfile.dict(cat)
        print '{0:<8s}  {1:.3f}  {2:.2f}'.format(cluster, z[i],
                                                 good['reff'].min())
        nudg[i] = good['reff'].size
        if nudg[i] > 100:
            texline = r'{0} & {1:8d} \\'.format(texline, nudg[i])
        else:
            texline = r'{0} & \,\,\,{1:d} \\'.format(texline, nudg[i])
        texlines.append(texline)
        size[i] = numpy.median(good['reff'])
        size_lo[i], size_hi[i] = numpy.absolute(numpy.percentile(
                    good['reff'], [16,84]) - size[i])
        #size_err[i] = (((good['reff']-size[i])**2).sum() / (2*nudg[i]))**0.5
        # color also by redshift?
        ax.plot(good['reff'], good['logmstar'], 'k.')
        full = readfile.dict(full)
        ax.plot(full['reff'], full['logmstar'], 'k+', ms=5, zorder=-1)
    # any clusters I want to mask (too low mass?)
    j = (m200 > 0)
    ncl = z[j].size
    print 'Using {0} UDGs in {1} clusters'.format(nudg[j].sum(), ncl)
    m200 *= 1e14
    m200_err *= 1e14

    averages(m200, m200_err, z, nudg)

    # to get some sort of error on the median
    size_lo /= (2*nudg)**0.5
    size_hi /= (2*nudg)**0.5
    # write tex table
    j = sort_clusters(clusters)
    outtex = 'clusters.tex'
    tex = open(outtex, 'w')
    for i in j:
        print >>tex, texlines[i]
    tex.close()
    # finish reff vs logmstar plot
    ax.set_xlabel(r'$r_\mathrm{eff}\,(\mathrm{kpc})$')
    ax.set_ylabel(r'$\log m_\star\,(\mathrm{M_\odot})$')
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    savefig('plots/reff_logmstar.pdf', fig)
    # plot m200 vs Nudg
    fig, ax = pylab.subplots(figsize=(6,4))
    ax.errorbar(m200, nudg, xerr=m200_err, fmt='ko', capsize=0)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(5e13, 2e15)
    ax.set_ylim(2, 200)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))
    ax.set_xlabel(r'$M_{200}\,(\mathrm{M_\odot})$')
    ax.set_ylabel(r'$N_\mathrm{UDG}(R<r_{200})$')
    savefig('plots/m200_nudg.pdf', fig)
    # plot m200 vs <reff>
    fig, ax = pylab.subplots(figsize=(5.5,5))
    ax.errorbar(m200, size, xerr=m200_err, yerr=(size_lo,size_hi),
                fmt='ko', capsize=0)
    ax.set_xscale('log')
    ax.set_xlim(5e13, 2e15)
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.set_xlabel(r'$M_{200}\,(\mathrm{M_\odot})$')
    ax.set_ylabel(r'$\langle r_\mathrm{eff} \rangle (\mathrm{kpc})$')
    savefig('plots/m200_meanreff.pdf', fig)
    return


def averages(m200, m200_err, z, nudg):
    w = 1 / m200_err**2
    mavg = (m200*w).sum() / w.sum()
    mavg_err = w.sum()**-0.5
    logmavg, logmavg_err = to_log(mavg, mavg_err)
    print 'inverse-variance:'
    print '  log Mavg = {0:.3f} +/- {1:.3f}'.format(logmavg, logmavg_err)
    w = nudg
    mavg = (m200*w).sum() / w.sum()
    mavg_err = 1 / w.sum()**0.5
    logmavg, logmavg_err = to_log(mavg, mavg_err)
    print 'Nudg-weighted:'
    print '  log Mavg = {0:.3f} +/- {1:.3f}'.format(logmavg, logmavg_err)
    w = nudg / m200_err**2
    mavg = (m200*w).sum() / w.sum()
    mavg_err = 1 / w.sum()**0.5
    logmavg, logmavg_err = to_log(mavg, mavg_err)
    print 'Both:'
    print '  log Mavg = {0:.3f} +/- {1:.3f}'.format(logmavg, logmavg_err)
    return


def savefig(output, fig=None, close=True, verbose=True, pad=0.4, **kwargs):
    if fig is None:
        fig = pylab
    fig.tight_layout(pad=pad)
    fig.savefig(output)
    if close:
        pylab.close()
    if verbose:
        print 'Saved to {0}'.format(output)
    return


def sort_clusters(clusters):
    ordered = numpy.argsort(clusters)
    abell_nr = [i for i in ordered if clusters[i][0] == 'A']
    abell = numpy.array([clusters[i][1:] for i in abell_nr], dtype=int)
    abell_sorted = numpy.argsort(abell)
    ordered[abell_nr] = abell_sorted
    return ordered

main()

