#!/usr/bin/env python
import numpy
import os
import pylab
import readfile
from glob import glob
from itertools import count, izip
from numpy import array

from astro.clusters import conversions
from kids_ggl_pipeline.halomodel import nfw
from kids_ggl_pipeline.halomodel.utils import density_average

# local
from models import delta, moline16
import pytools

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

"""
Assume a subhalo-to-stellar mass relation, calculate the mstar-weighted
signal and extract best-fit parameters from it - how accurate is the
mass inferred from the weighted signal?

"""

def main():
    path_data = 'data/logmstar/logmstar_7.0_10.0-all-z_0.00_0.15'
    path_data += '-sigma_0.25-zs_0.60-mrsource_20.0_24.5-rmax_2.0-maxsize_100'
    path_data += '-minsep_0-scale_logmstar_8.3_1_exp10-delmag-manual_fiducial'

    # get total masses assuming some SHSMR
    logmstar = readfile.table(os.path.join(path_data, 'hist-logmstar.dat'))
    msat = shsmr(logmstar[0])
    # cluster-centric distances
    Rsat = readfile.table(os.path.join(path_data, 'hist-logdistBCG.dat'))
    # host masses - lognormal distribution
    Mhost = [numpy.linspace(14.5, 15.1, 10)]
    Mhost.append(10**numpy.exp(-(Mhost[0]-14.7)**2/(2*0.2**2)))
    Mhost = array(Mhost)

    # make a grid in msat and Rsat to calculate the final signal
    # this assumes that there is no stellar mass segregation - almost true?
    #mgrid, Rgrid = numpy.meshgrid(msat, Rsat[0])
    weights = logmstar[1] * Rsat[1][:,numpy.newaxis]
    weights /= weights.sum()
    #extent = (logmstar[0].min(),logmstar[0].max(),Rsat[0].min(),Rsat[0].max())
    #pylab.imshow(weights, cmap='cubehelix_r', origin='lower', extent=extent,
                 #interpolation='none', aspect='auto')
    #pylab.colorbar()
    #pylab.show()

    R = numpy.logspace(-1.7, 0, 10)
    #calculate_signal(mgrid, Rgrid, Mhost)
    fig, ax = pylab.subplots(figsize=(5,4))
    #for c in (5, 10, 20, 40):
        #esd = calculate_signal(R, msat, Rsat[0], Mhost, weights, csat=c)
        #ax.plot(R, esd, label='c={0}'.format(c))
    esd_tot, esd_ind = calculate_signal(R, msat, Rsat[0], Mhost, weights)
    for i in esd_ind:
        ax.plot(R, i, color='0.5', alpha=0.5)
    ax.plot(R, esd_tot, color=blue, lw=3)
    ax.set_xscale('log')
    ax.legend(loc='upper right', frameon=False)
    fig.tight_layout()
    pylab.show()

    # now simulate data points with noise



    return


def calculate_signal(R, msat, Rsat, Mhost, weights, z=0.1,
                     csat=20, chost=5, h=1, Om=0.315):
    """
    Assuming the same concentration for all subhaloes for simplicity.

    ESDs do not include high-kappa correction

    """
    Ol = 1 - Om
    rho_m = density_average(z, h, Om, Ol)
    # lens-source separation bins
    # satellite signal - no need to run for each Rsat
    r200sat = conversions.rsph(msat, z, ref='200a')
    rssat = r200sat / csat
    Sigma_s = rssat * delta(csat) * rho_m
    # shape = (logmstar.size,R.size)
    esd_sat = array([nfw.esd(R/rs, sigma)
                     for rs, sigma in izip(rssat, Sigma_s)])
    print 'esd_sat =', esd_sat.shape
    print 'weights =', weights.shape

    # cluster signal, for a single cluster mass for now
    #r200host = conversions.rsph(Mhost, z, ref='200a')
    #rshost = rssat / chost
    #sigma_s = rshost * delta(chost) * rho_m
    #rsat_range = readfile.table('rsat_range.txt')
    #esd_host = array([nfw.esd_offset(R/rs, 10**Rsat[0]/rs, Rsat[1],

    # since I'm not using Rsat for now
    weights = numpy.sum(weights, axis=0)
    esd_sat_tot = numpy.sum(esd_sat*weights[:,numpy.newaxis], axis=0) / \
        weights.sum()
    print 'esd_sat =', esd_sat.shape
    return esd_sat_tot, esd_sat


def shsmr(logmstar, a=11.9, b=0.8, pivot=2e10):
    """
    extend to a double power law later
    
    """
    return 10**a * (10**logmstar/pivot)**b


main()


