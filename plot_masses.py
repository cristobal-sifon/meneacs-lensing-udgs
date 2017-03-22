#!/usr/bin/env python
import numpy
import os
import pylab
import readfile
import sys
from itertools import count, izip
from matplotlib import ticker

from astro.clusters import conversions
from kids_ggl_pipeline.halomodel import nfw
from kids_ggl_pipeline.halomodel.utils import density_average

from astro import cosmology
cosmology.h = 0.7
cosmology.Omega_M = 0.315
cosmology.Omega_L = 1 - cosmology.Omega_M

import plottools
plottools.update_rcParams()

# local
from model_components.concentration import moline16
import pytools

# until a newer matplotlib is available
import colormaps

red = (1,0,0)
green = (0.2,0.6,0)
blue = (0,0,1)
yellow = (1,1,0.4)
magenta = (1,0.4,0.6)
cyan = (0,1,1)
orange = (1,0.7,0)
purple = (0.8,0,0.4)


def main():
    file = sys.argv[1]
    hdr = file.replace('.out', '.hdr')
    hdr = open(hdr)
    for line in hdr:
        if line[:9] == 'renclosed':
            r = numpy.array(line.split()[-1].split(','), dtype=float)
            break
    hdr.close()
    mcmc = readfile.dict(file, include='Menclosed',
                         cols=('median','p95'))
    mass = 10**mcmc['median']
    p95 = 10**mcmc['p95']
    plot_path = os.path.join('mcmcplots', file.split('/')[-1][:-4])
    plot(mass, p95, r, plot_path)
    return
    

def plot(mass, p95, radius, plot_path):
    """
    Use different symbols for different techniques

    There seems to be a bug in errorbar() when using limits
    """
    fig, ax = pylab.subplots(figsize=(5.5,4))
    # DF17
    df17x = [2.4, 8.1]
    df17y = [2.6e9, 4.5e9]
    df17xerr = [0.1, 0.1]
    df17yerr = [[[1.3e9],[3.4e9]], 2.8e9]
    ax.errorbar(df17x[0], df17y[0], xerr=df17xerr[0], yerr=df17yerr[0],
                fmt='s', mec='k', ecolor='k', mfc='none', mew=1, lw=1,
                capsize=0, label='DF17 (B+16)')
    ax.errorbar(df17x[1], df17y[1], yerr=df17yerr[1], xerr=df17xerr[1],
                fmt='s', mec='k', ecolor='k', mfc='none', mew=1, lw=1,
                capsize=0, label='_none_')
    ax.plot(df17x, df17y, 'k--', lw=1)
    # DF44, velocity dispersion (van Dokkum et al. 2016)
    ax.errorbar(4.6, 7e9, xerr=0.2, yerr=((2e9,),(3e9,)), fmt='o',
                mec='k', ecolor='k', mfc='k', mew=1, lw=1, capsize=0,
                #label='_none_')
                label=r'DF44 (vD+16)')
    yerr = p95 * numpy.ones(mass.size)
    ax.errorbar(radius, mass, yerr=yerr-mass, fmt=',', color=red, ecolor=red,
                lw=2, mew=2, ms=12, capsize=5, lolims=True, zorder=10,
                label='_none_')
    # bug in errorbar?
    verts = [[0,0], [0,-0.8], [-0.15,-0.6], [0,-0.8], [0.15,-0.6], [0,-0.8]]
    ax.scatter(radius, p95, color=red, marker=None, s=1000, verts=verts, lw=2,
               zorder=10, label='_none_')
    # NFW profiles
    t = numpy.logspace(-3, -0.1, 1000)
    m200 = numpy.array([1e11, 5e11, 1e12])
    if not numpy.iterable(m200):
        m200 = numpy.array([m200])
    z = 0.023
    rho_m = density_average(z, h=cosmology.h, Om=cosmology.Omega_M,
                            Ol=cosmology.Omega_L)
    r200 = conversions.rsph(m200, z, ref='200a')
    verts_r200 = [[1,-1], [0,0], [-1,-1], [0,0], [0,-5], [0,0]]

    ## first option: various concentrations, fixed mass
    #m200 = m200[2]
    #c = numpy.array([5, 10, 20])
    #for ci, dashes in izip(c, ((5,5), (10,0), (10,8))):
        #rs = r200 / ci
        #for rsi, r, m in izip(rs, r200, m200):
            #note = r'$\log m_{{200}}/\mathrm{{M}}_\odot={0:.1f}$'.format(
                #numpy.log10(m))
            #ax.annotate(note, xy=(0.95,0.05), xycoords='axes fraction',
                        #fontsize=16, ha='right', va='bottom')
            #sigma = rsi * nfw.delta(ci) * rho_m
            #menc = nfw.mass_enclosed(t/rsi, rsi, sigma)
            #ax.plot(1e3*t, menc, 'k', dashes=dashes, lw=1)
            ##mtot = nfw.mass_enclosed(ci, rsi, sigma)
            ##ax.plot(1e3*ci*rsi, mtot, 'kx', ms=10, mew=2)
            #ax.scatter(1e3*r, 0.9*m, color='k', marker=None, s=200,
                       #verts=verts_r200, lw=1, label='_none_')
        #ax.plot([], [], 'k', dashes=dashes, label='$c={0}$'.format(ci))

    # second option: various masses with concentration scatter
    # concentration from Moline et al. (2016)
    # for m = (1e11, 5e11, 1e12) it's c = (12.8, 10.7, 9.9)
    cscatter = 10**0.1
    for i, m, r in izip(count(), m200, r200):
        cmol = moline16(m, 0.38)
        mmol = mass_profile(t, m, r, cmol, rho_m)
        ax.plot(1e3*t, mmol, 'k-', lw=1,# label='c={0:.1f}'.format(cmol))
                label='_none_')
        mlo = mass_profile(t, m, r, cmol/cscatter, rho_m)
        mup = mass_profile(t, m, r, cmol*cscatter, rho_m)
        color = str(0.2*(i+2))
        ax.fill_between(1e3*t, mlo, mup, color=color, lw=0, zorder=-10,
                        label='_none_')
    # finish up
    ax.legend(loc='upper left', numpoints=1, frameon=False)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1, 500)
    ax.set_ylim(1e9, 2e12)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%d$'))
    ax.set_xlabel(r'$r\,(\mathrm{kpc})$')
    ax.set_ylabel(r'$M(<r)\,(\mathrm{M}_\odot)$')
    fig.tight_layout(pad=0.4)
    pytools.savefig(os.path.join(plot_path, 'mass_literature.pdf'), fig)
    return


def mass_profile(t, m200, r200, c, rho_m):
    rs = r200 / c
    sigma = rs * nfw.delta(c) * rho_m
    return nfw.mass_enclosed(t/rs, rs, sigma)

if __name__ == '__main__':
    main()


