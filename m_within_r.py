#!/usr/bin/env python
import numpy
import pyfits
import pylab
from astro import cosmology
from astro.clusters import conversions
from astro.ggl import nfw
from matplotlib import colors as mplcolors, cm, ticker

import lnr

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
#blue = '#2D7EDF'
#red = '#D34A1E'
#green = '#36911C'
#purple = '#650D52'
#yellow = '#F6B91F'
cNorm = mplcolors.Normalize(vmin=0, vmax=1)
cmap = pylab.get_cmap('viridis')
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)

filename = 'output/'
#filename += 'fullnfw-distBCG_0.0_1.0-all-z_0.00_0.15'
#filename += '-mrsource_20.0_24.5-delmag.fits'
filename += 'fullnfw-redshift_0.0_0.5-z_0.00_0.15-rmax_3.0.fits'
chain = pyfits.getdata(filename, ignore_missing_end=True)
print chain.names

def concentration(m, z):
    return 10.14 / (m/2e12)**0.081 / (1+z)**1.01

z = 0.05
print chain['Msat1'].size
m = chain['Msat1'][100000:]
c = concentration(m, z)
r200 = conversions.rsph(m, z)
rs = r200 / c
sigma_s = nfw.delta(c) * rs * cosmology.density(z, ref='average')
#print m, c, r200, rs
print numpy.median(rs)

R = numpy.array([10, 20]) / 1e3
x = numpy.array([Ri / rs for Ri in R])
mass = numpy.zeros((R.size,m.size))
for i in xrange(R.size):
    mass[i] = 4 * 3.14159265 * rs**2 * sigma_s * \
              (numpy.log(1+x[i]) - x[i] / (1+x[i]))
total_to_stellar =  mass / 10**8.2
#print mass

print ' ** Total masses **'
w = numpy.ones(m.size) / m.size
colors = scalarMap.to_rgba(numpy.linspace(0.1, 0.9, R.size))
fig, ax = pylab.subplots(figsize=(5.5,5))
for i in xrange(R.size):
    m =  numpy.median(mass[i])
    logm = numpy.log10(m)
    err = numpy.absolute(m - numpy.percentile(mass[i], [16,84]))
    print R[i], numpy.log10(m),
    #print numpy.log10(numpy.absolute(m - numpy.percentile(mass[i], [16,84])))
    #x1err = numpy.log10(numpy.array(x1)+numpy.array(x1err)) - logx1
    print -numpy.log10(m - err[0]) + logm,
    print numpy.log10(m + err[1]) - logm
    ax.hist(numpy.log10(mass[i]), bins=100, histtype='step', color=colors[i],
            label=r'$<%.0f\,{\rm kpc}$' %(1e3*R[i]), weights=w, lw=3)
#ax.hist(numpy.log10(m), bins=100, histtype='step', color=red, lw=3,
        #label='%.0f kpc' %(1e3*numpy.median(r200)), weights=w)
#ax.set_xscale('log')
arrowprops = dict(arrowstyle='->', linewidth=3, color=blue)
ax.annotate('', xy=(8.2,0), xytext=(8.2,0.018), arrowprops=arrowprops)
ax.set_xlim(8, 11)
ax.set_ylim(0, 0.1)
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax.set_xlabel(r'$\log [M(<R)/M_\odot]$')
ax.set_ylabel(r'$p(\log M)$')
fig.tight_layout(pad=0.4)
ax.legend(loc='upper left')
output = 'plots/m_within_r.png'
fig.savefig(output)
pylab.close()

print ' ** Mass ratios **'
fig, ax = pylab.subplots(figsize=(5.5,5))
for i in xrange(R.size):
    m =  numpy.median(total_to_stellar[i])
    print R[i], m,
    print numpy.absolute(m - numpy.percentile(total_to_stellar[i], [16,84]))
    ax.hist(total_to_stellar[i], bins=25,
            histtype='step', color=colors[i],
            label=r'$<%.0f\,{\rm kpc}$' %(1e3*R[i]), weights=w, lw=3)
    line = ax.axvline(m, color=colors[i], lw=2)
    line.set_dashes([10,10])
ax.set_ylim(0, 0.12)
#ax.set_xlabel(r'$\log [M(<R)/M_{\rm stellar}]$')
#ax.set_ylabel(r'$p(\log M/M_{\rm stellar})$')
ax.set_xlabel(r'$M(<R)/M_{\rm stellar}$')
ax.set_ylabel(r'$p(M/M_{\rm stellar})$')
fig.tight_layout(pad=0.4)
ax.legend(loc='upper right')
output = 'plots/mratio_within_r.png'
fig.savefig(output)
pylab.close()