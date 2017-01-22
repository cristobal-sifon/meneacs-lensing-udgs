#!/usr/bin/env python
import numpy
import os
import pylab
import readfile
import sys
from glob import glob
from matplotlib import ticker

# local
import pytools
import utils
import lensing_signal as signal

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
rcParams['lines.linewidth'] = 3
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
    root = 'redshift_0.04_0.07_0.1-all-z_0.00_0.15-sigma_0.25-zs_0.60'
    root += '-mrsource_20.0_24.5-rmax_2.0-maxsize_100-minsep_0'
    root += '-scale_logmstar_8.3_1_exp10-delmag-manual_fiducial'
    zbins = root.split('-')[0].split('_')[1:]
    zfull = '_'.join((zbins[0], zbins[-1]))
    zrange = '_'.join(zbins)
    #colors = utils.get_colors(vmin=0.2, vmax=0.6, n=2, cmap='viridis')
    colors = (green, orange)
    plot_path = os.path.join('plots', 'histograms_by_redshift',
                             'redshift_{0}'.format(zrange))
    if not os.path.isdir(plot_path):
        os.makedirs(plot_path)
    print 'Saving histograms to {0}'.format(plot_path)
    histograms = glob(os.path.join('data', 'redshift', root, 'hist-*'))
    for file1 in histograms:
        value = file1.split('/')[-1].split('-')[1][:-4]
        file2 = file1.replace(zrange, zfull)
        if not os.path.isfile(file2):
            continue
        hist1 = readfile.table(file1)
        if len(hist1) == 2:
            continue
        hist2 = readfile.table(file2)
        fig, ax = pylab.subplots(figsize=(5,4))
        ax.step(hist2[0], hist2[1], where='mid', color='k', lw=1, zorder=10)
        for i, color in enumerate(colors, 1):
            ax.step(hist1[0], hist1[i], where='mid', color=color,
                    #dashes=(4*(i+1),4*(i+1)),
                    label=r'${0} \leq z < {1}$'.format(zbins[i-1], zbins[i]))
        #ax.step(hist2[0], hist2[1], where='mid', color='k', dashes=(10,10))
        ax.legend(loc='upper right', frameon=False, fontsize=14)
        label = signal.get_label(value)
        ax.set_xlabel(label)
        ax.set_ylabel(r'$N({0})$'.format(label.replace('$', '')))
        ylim = get_ylim(value)
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(6))
        ax.xaxis.set_minor_locator(ticker.MaxNLocator(24))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
        ax.yaxis.set_minor_locator(ticker.MaxNLocator(20))
        fig.tight_layout(pad=0.4)
        pytools.savefig(os.path.join(plot_path, '{0}.pdf'.format(value)), fig,
                        verbose=False)
    return


def get_ylim(param):
    ylim = {'reff': (0, 200)}
    if param in ylim:
        return ylim[param]
    return
    

if __name__ == '__main__':
    main()


