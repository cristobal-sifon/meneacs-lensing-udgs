#!/usr/bin/env python
import numpy
import pylab

import plottools
from astro import cosmology
from astro.clusters import conversions
from kids_ggl_pipeline.halomodel import nfw
from kids_ggl_pipeline.halomodel.utils import density_average
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


def plot(mr, r, z, ax=None, cbins=numpy.linspace(10, 20, 11),
         m200bins=numpy.logspace(11, 12.3, 9)):
    """
    Make a grid of (c,m200) that covers the given [r,m(<r)].
    Could return the chi2 for each (c,m200) pair to draw confidence
    intervals.

    Note that r must be in Mpc
    """
    if ax is None:
        ax = pylab
    logm200bins = numpy.log10(m200bins)
    # calculate the enclosed mass for each (c,m200) pair
    mr_pred = numpy.zeros((cbins.size,m200bins.size))
    r200 = conversions.rsph(m200bins, z, ref='200a')
    rho_m = density_average(z, h=0.7, Om=0.315, Ol=0.685)
    for i, c in enumerate(cbins):
        rs = r200 / c
        sigma_s = rs * nfw.delta(c) * rho_m
        mr_pred[i] = nfw.mass_enclosed(r/rs, rs, sigma_s)
    extent = (cbins[0], cbins[-1], logm200bins[0], logm200bins[-1])
    fig, ax = pylab.subplots(figsize=(6,4))
    pylab.imshow(numpy.log10(mr_pred), origin='lower', extent=extent,
                 aspect='auto',
                 cmap=colormaps.viridis, interpolation='nearest')
    pylab.colorbar()
    fig.tight_layout(pad=0.4)
    output = os.path.join('plots/esd_overplot.pdf')
    plottools.savefi(output, fig=fig)
    return


if __name__ == '__main__':
    plot(7e9, 4.6e-3, z=0.023)


