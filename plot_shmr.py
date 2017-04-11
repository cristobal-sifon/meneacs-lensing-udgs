#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

import argparse
import numpy
import os
import pylab
from numpy import log10
from scipy.interpolate import interp1d
from uncertainties import unumpy

# my code
import plottools
import readfile
plottools.update_rcParams()

# local
from literature import leauthaud12


def main(xlim=(5e7,5e11)):
    args = parse_args()

    # average stellar mass in our sample
    mstar_udg = 2e8

    use_bound = (True, True, False, False)
    refs = ('leauthaud12-1,sifon17,mw',
            'leauthaud12-1,sifon17,udgs,mw',
            'leauthaud12-1,sifon17,udgs',
            'leauthaud12-1,sifon17,udgs,mw')
    include_m200 = (False, True, False, False)
    for ref, bound, show_m200 in zip(refs, use_bound, include_m200):
        # load data
        if bound:
            # this one includes too many uncertainties
            include = 'Msat_rbg'
            #include = 'Menclosed3'
        else:
            include = 'Msat'
        # median, e16, e84, p95 = udg
        udg = readfile.table(
            args.filename, include=include, cols=(2,3,4,5), whole=True)
        udg = numpy.append(mstar_udg, udg)
        print(include, udg[0], udg[4])
        if bound and show_m200:
            udg200 = readfile.table(
                args.filename, include='Msat', cols=(2,3,4,5), whole=True)
            udg200 = numpy.append(mstar_udg, udg200)
        else:
            udg200 = None
        do_plot(args, udg, xlim, ref, udg200=udg200, use_bound=bound)
        print()

    return


def do_plot(args, udg, xlim, ref, udg200=None, use_bound=True):
    if use_bound:
        label = 'm_\mathrm{{bound}}'
        name = 'mbound'
    else:
        label = 'm_{{200}}'
        name = 'm200'
    fig, ax = pylab.subplots(figsize=(5.5,5))
    x = numpy.logspace(numpy.log10(xlim[0]), numpy.log10(xlim[1]), 100)
    curves, curve_labels = literature(
        ax, ref, use_bound=use_bound, show_m200=(udg200 is not None))
    print('log udg = {0}'.format(udg[4]))
    pt = ax.errorbar(
        udg[0], 10**udg[4], yerr=0.7*10**udg[4], fmt=',', uplims=[True],
        color='C0', capsize=5, mew=3, elinewidth=3)
    curves[curves.index(None)] = pt
    curve_labels[curve_labels.index(None)] = 'MENeaCS UDGs (this work)'
    print('curve_labels =', curve_labels)
    if udg200 is not None:
        print('udg200 =', udg200)
        avg = (10**udg200[4]+10**udg[4]) / 2
        yerr = numpy.array([[avg-10**udg[4]],[10**udg200[4]-avg]])
        print('avg =', numpy.log10(avg), udg200[1], udg[1])
        print('yerr =', yerr, numpy.log10(yerr))
        ebar = ax.errorbar(
            udg200[0], avg, yerr=yerr,
            fmt=',', color='C0', capsize=5, mew=3, elinewidth=3)
        # eb1[-1][0] is the LineCollection object of the errorbar lines
        ebar[-1][0].set_linestyle('--')
        curves[curves.index(None)] = ebar
        curve_labels[curve_labels.index(None)] = \
            'MENeaCS UDGs m$_{200}$ (this work)'
    ax.legend(curves, curve_labels, loc='upper left', frameon=False,
              fontsize=14)

    ax.set_xlim(*xlim)
    ax.set_ylim(1e10, 1e14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$m_\star\,(\mathrm{M}_\odot)$')
    ax.set_ylabel(r'${0}\,(\mathrm{{M}}_\odot)$'.format(label))

    output = os.path.join(
        args.filename.replace('output/', 'mcmcplots/')[:-4],
        'shmr_literature_{0}_{1}.pdf'.format(name, len(ref.split(','))))
    if udg200 is not None:
        output = output.replace('.pdf', '_withm200.pdf')
    print('output =', output)
    plottools.savefig(output, fig=fig)
    print()

    return


def literature(ax, ref='leauthaud12-1,sifon17', use_bound=True, labelsize=17,
               show_m200=False, h=0.7):
    """
    Adapted from mcmcplots.plot_literature()

    """
    def log2lin(logm, dlogm_lo, dlogm_hi):
        m = 10**logm
        mlo = m - 10**(logm - dlogm_lo)
        mhi = 10**(logm + dlogm_hi) - m
        return m, mlo, mhi

    refs = ref.split(',')
    curves = []
    curve_labels = []
    if show_m200:
        n = 2
    else:
        n = 1
    indices = -numpy.ones(len(refs)+n, dtype=int)
    # Leauthaud et al. (2012) - COSMOS GGL + clustering + SMF
    print('indices =', indices)
    if 'leauthaud12' in ref:
        zbin = [i[-1] for i in refs if 'leauthaud12' in i][0]
        zleau = {'1': 0.36, '2': 0.66, '3': 0.88}
        label = 'Centrals (Leauthaud+12, z={0})'.format(zleau[zbin])
        x = numpy.logspace(9, 11.7, 100) * (0.72/h)**2
        y = 10**leauthaud12.shmr(x, zbin) * (0.72/h)
        yerr = unumpy.std_devs(y)
        y = unumpy.nominal_values(y)
        curves.append(ax.plot(x, y, 'k-', lw=2)[0])
        xlo = numpy.logspace(8, 9, 10) * (0.72/h)**2
        ylo = 10**leauthaud12.shmr(xlo, zbin) * (0.72/h)
        yloerr = unumpy.std_devs(ylo)
        ylo = unumpy.nominal_values(ylo)
        ax.plot(xlo, ylo, 'k--', lw=2)
        curve_labels.append(label)
        indices[refs.index('leauthaud12-{0}'.format(zbin))] = 0
    # Zu & Mandelbaum (2016, MNRAS, 457, 4360)
    print('indices =', indices)
    if 'zu16' in ref:
        filename = [os.path.join(path,
                                 'zu_mandelbaum_haloquench_%s.csv' %i)
                    for i in ('red', 'blue')]
        cent_red = readfile.table(filename[0])
        curves.append(
            ax.plot(cent_red[0]/h**2, cent_red[1]/h, '-',
                    color=red, zorder=1)[0])
        curve_labels.append('Z&M16 Red')
        indices[refs.index('zu16')] = indices.max()+1
        #cent_blue = readfile.table(filename[1])
        #ax.plot(cent_blue[0]/h**2, cent_blue[1]/h, '-',
                #color=blue, zorder=1, label='Z\&M16 Blue')
        #curves.append(interp1d(log10(cent_blue[0]/h**2),
                               #log10(cent_blue[1]/h)))
        #curve_labels.append('Z&M16 Blue')
    # My MENeaCS paper
    print('indices =', indices)
    if 'sifon17' in ref:
        label = r"MENeaCS satellites (Sif\'on+17)"
        # update to use the latest results once I have them
        # this table only has bound masses
        filename = 'literature/sifon17.txt'
        s17 = readfile.table(filename)
        s17[0] = 10**s17[0]
        s17[1:] = log2lin(*s17[1:])
        curves.append(
            ax.errorbar(s17[0], s17[1], yerr=(s17[2],s17[3]), fmt='ko',
                        ms=8, capsize=2, elinewidth=2, mew=2))
        curve_labels.append(label)
        indices[refs.index('sifon17')] = indices.max()+1
    #index_this_work = numpy.argmax(indices)
    #indices[index_this_work+1] = indices.max()+1
    # other UDGs
    print('indices =', indices)
    if 'udgs' in ref:
        ## DF44 (vDokkum+16)
        # extrapolate to 30 kpc only
        if use_bound:
            pass
        # use m200 extrapolation from van Dokkum
        else:
            curves.append(ax.plot(3e8, 1e12, 'wo', mec='0.3', ms=7, mew=2)[0])
            if use_bound:
                label = 'DF44 m$_{200}$ (van Dokkum+16)'
            else:
                label = 'DF44 (van Dokkum+16)'
        curve_labels.append(label)
        #ax.annotate('DF44', xy=(3.5e8,1.05e12), ha='left', va='bottom',
                    #color='0.3', fontsize=labelsize, fontweight='heavy')
        indices[refs.index('udgs')+n] = indices.max()+1
        # DF17 (Beasley+16)
        #ax.errorbar(2.8e7, 8e10, xerr=0.4e7, yerr=4e10, fmt=
    # Milky Way
    print('indices =', indices)
    if 'mw' in ref:
        curves.append(
            ax.plot(5e10, 1e12, '*', mec='0.3', mfc='w', ms=12, mew=2)[0])
        curve_labels.append('Milky Way')
        #ax.annotate('MW', xy=(5e10,6e11), ha='center', va='top', color='0.3',
                    #fontsize=labelsize, fontweight='heavy')
        indices[refs.index('mw')+n] = indices.max()+1
    print('indices =', indices)

    print('refs =', refs)
    print('curve_labels =', curve_labels)
    print('indices =', indices)
    print('len(curves) =', len(curves))
    curves = [curves[i] if i > -1 else None for i in indices]
    curve_labels = [curve_labels[i] if i > -1 else None for i in indices]
    print('curve_labels =', curve_labels)

    return curves, curve_labels


def axlabel(key):
    labels = {'Msat': '$m_{200}$',
              'Msat_rbg': r'$m_\mathrm{bound}$'}
    for i in xrange(1, 4):
        labels['Menclosed{0}'.format(10*i)] = \
            r'$m(<{0}\,\mathrm{{kpc}})$'.format(10*i)
    if key in labels:
        return labels[key]
    return key


def parse_args():
    parser = argparse.ArgumentParser()
    add = parser.add_argument
    add('filename')
    args = parser.parse_args()
    return args


main()


