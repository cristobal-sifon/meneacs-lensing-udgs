#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import argparse
import numpy
import os
import pylab
from astropy.cosmology import FlatLambdaCDM
from numpy import array, log10
from scipy.interpolate import interp1d
from scipy import optimize
from uncertainties import unumpy

# KiDS-GGl pipeline
from kids_ggl_pipeline.halomodel import nfw

# my code
import plottools
import readfile
from astro import cosmology
from astro.clusters import conversions, profiles
plottools.update_rcParams(
    {'text.usetex': True})

# local
from literature import ihod, leauthaud12

#cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
cosmology.h = 0.7
cosmology.Omega_M = 0.315
cosmology.Omega_L = 1 - cosmology.Omega_M


def main(xlim=(1e7,5e11)):
    args = parse_args()

    # average stellar mass in our sample
    mstar_udg = 2e8

    #refs = ('leauthaud12-1,sifon17,mw',
            #'leauthaud12-1,sifon17,udgs,mw',
            #'leauthaud12-1,sifon17,udgs',
            #'leauthaud12-1,sifon17,udgs,mw')
    #include_m200 = (False, True, False, False)
    # keys are the quantities to be plotted; each element is
    # (references
    layouts = {'Msat-1': 'leauthaud12-1,sifon17,mw',
               'Msat-2': 'leauthaud12-1,sifon17,udgs,mw',
               'Msat-3': 'zu15,sifon17,udgs,mw',
               'Msat_rbg': 'sifon17',
               'Menclosed3-1': 'sifon17,udgs',
               'Menclosed3-2': 'sifon17,udgs,mw'}
    for mass, ref in layouts.items():
        # load data
        include = mass.split('-')[0]
        # median, e16, e84, p95 = udg
        udg = readfile.table(
            args.filename, include=include, cols=(2,3,4,5), whole=True)
        udg = numpy.append(mstar_udg, udg)
        print(include, udg[0], udg[4])
        do_plot(args, udg, mass, xlim, ref)
        print()

    return


def do_plot(args, udg, name, xlim, ref, udg200=None, use_bound=True):
    ylabel = axlabel(name.split('-')[0])
    ylim = get_ylim(name, ref)
    fig, ax = pylab.subplots(figsize=(5.5,5))
    x = numpy.logspace(numpy.log10(xlim[0]), numpy.log10(xlim[1]), 100)
    curves, curve_labels = literature(
        ax, name, udg, ref, xlim=xlim, show_m200=(udg200 is not None))
    print('log udg = {0}'.format(udg[4]))
    if ylim[1] <= 1e13:
        errscale = 0.5
    else:
        errscale = 0.5
    pt = ax.errorbar(
        udg[0], 10**udg[4], yerr=errscale*10**udg[4], fmt=',', uplims=[True],
        color='C3', capsize=5, mew=3, elinewidth=3)
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
    ax.set_ylim(*ylim)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$m_\star\,(\mathrm{M}_\odot)$')
    ax.set_ylabel(r'${0}\,(\mathrm{{M}}_\odot)$'.format(ylabel))

    try:
        index = name.split('-')[1]
    except IndexError:
        index = 1
    output = os.path.join(
        args.filename.replace('output/', 'mcmcplots/')[:-4],
        'shmr_literature_{0}_{1}.pdf'.format(
            name.split('-')[0].lower(), index))
    if udg200 is not None:
        output = output.replace('.pdf', '_withm200.pdf')
    plottools.savefig(output, fig=fig)
    print()

    return


def literature(ax, name, udg, ref='leauthaud12-1,sifon17', labelsize=17,
               xlim=(5e7,1e12), show_m200=False, show_tidal=True, h=0.7):
    """
    Adapted from mcmcplots.plot_literature()

    """
    def log2lin(logm, dlogm_lo, dlogm_hi):
        m = 10**logm
        mlo = m - 10**(logm - dlogm_lo)
        mhi = 10**(logm + dlogm_hi) - m
        return m, mlo, mhi

    # quick fix to 3 UDGs
    if ',udgs' in ref:
        if show_tidal:
            nudg = 6
        elif 'Menclosed' in name:
            nudg = 5
        else:
            nudg = 0
        ref = ref.replace(',udgs', ',udgs'*nudg)
    else:
        nudg = 0
    refs = ref.split(',')
    print('refs =', refs)
    curves = []
    curve_labels = []
    if show_m200:
        n = 2 #* nudg
    else:
        n = 1 #* nudg
    indices = -numpy.ones(len(refs)+n, dtype=int)
    # Leauthaud et al. (2012) - COSMOS GGL + clustering + SMF
    if 'leauthaud12' in ref:
        zbin = [i[-1] for i in refs if 'leauthaud12' in i][0]
        zleau = {'1': 0.36, '2': 0.66, '3': 0.88}
        label = 'Centrals (Leauthaud+12, z={0})'.format(zleau[zbin])
        x = numpy.logspace(9, 11.7, 100) * (0.72/h)**2
        y = 10**leauthaud12.shmr(x, zbin) * (0.72/h)
        yerr = unumpy.std_devs(y)
        y = unumpy.nominal_values(y)
        curves.append(ax.plot(x, y, 'k-', lw=2)[0])
        xlo = numpy.logspace(7.3, 9, 10) * (0.72/h)**2
        ylo = 10**leauthaud12.shmr(xlo, zbin) * (0.72/h)
        yloerr = unumpy.std_devs(ylo)
        ylo = unumpy.nominal_values(ylo)
        ax.plot(xlo, ylo, 'k--', lw=2)
        curve_labels.append(label)
        indices[refs.index('leauthaud12-{0}'.format(zbin))] = 0
    if 'zu15' in ref:
        # this is the range probed by Z&M15
        logx = numpy.linspace(9.8, 12.3, 100)
        logy = ihod.logMh(logx, h=h)
        curves.append(ax.plot(10**logx, 10**logy, color='0.1', lw=2)[0])
        curve_labels.append('Centrals (Zu\&Mandelbaum15)')
        indices[refs.index('zu15')] = 0
        # continue a little to lower mass
        logx = numpy.linspace(7.3, 9.8, 20)
        logy = ihod.logMh(logx, h=h)
        ax.plot(10**logx, 10**logy, color='0.1', lw=2, ls='--')
    # Zu & Mandelbaum (2016, MNRAS, 457, 4360)
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
    # Edo's KiDS paper
    #if 'vanuitert16' in ref:
    # My MENeaCS paper
    if 'sifon17' in ref:
        label = r"MENeaCS satellites (Sif\'on+17)"
        filename = '../satellites/output_20170524/' \
                   'fullnfw_moline17_uniform_cMhost' \
                   '-logmstar_9.0_9.8_10.2_10.5_10.9_11.2-all-delmag' \
                   '-exclude_9_10_11_12.out.backup'
        #filename = 'literature/sifon17.txt'
        print(name)
        s17 = array(readfile.table(filename, include=name.split('-')[0]))
        s17_names = s17[0]
        if '_rbg' in name:
            s17 = array([i[2:5] for i in s17.T if '_rbg' in i[0]],
                        dtype=float).T
        else:
            s17 = array([i[2:5] for i in s17.T if '_rbg' not in i[0]],
                        dtype=float).T
        print('s17 =', s17)
        s17_mstar = readfile.table(
            filename.replace('.out', '.hdr'), include='Mstar', cols=2)
        s17_mstar = array(s17_mstar[0].split(','), dtype=float)
        s17 = log2lin(*s17)
        curves.append(
            ax.errorbar(s17_mstar, s17[0], yerr=(s17[1],s17[2]), fmt='ko',
                        ms=8, capsize=2, elinewidth=2, mew=2))
        curve_labels.append(label)
        indices[refs.index('sifon17')] = indices.max()+1

    ## other UDGs

    if 'udgs' in refs:
        # remember that I need to modify nudg above whenever I add or
        # remove objects here
        udgcolor = '0.3'
        # Make sure that the last element corresponds to the adopted
        # concentration
        c = array([5, 10])
        # Skipping these for now until I get confirmation from Nicola
        """
        # first, upper limit on average in Coma from globular cluster
        # counts. Converted m200c to m200a including 0.15 dex scatter
        # from Hudson+14 CFHTLenS, which was used to calibrate Ngc
        label = 'Coma UDGs (Amorisco+17)'
        mstar = 4.4e7
        zcoma = 0.023
        m200a = profiles.nfw(3.09e10, zcoma, ref_in='200c', ref_out='200a')
        r200a = conversions.rsph(m200a, zcoma, ref='200a')
        for ci in c:
            mass = udgmass(m200a, r200a, ci, zcoma, name)[0]
        coma = ax.errorbar(
            mstar, mass, yerr=0.5*mass, fmt=',', uplims=[True],
            color=udgcolor, mew=2, elinewidth=2, capsize=3)
        if 'Menclosed' in name:
            curves.append(coma)
            curve_labels.append(r'Coma UDGs (Amorisco+17)')
            indices[refs.index('udgs')+n] = indices.max()+1
        """

        ##
        zcoma = 0.023
        label = [r'DF17 (Beasley\&Trujillo16, Peng\&Lim16)',
                 r'DF44 (van Dokkum+16)',
                 r'VCC1287 (Beasley+16)',
                 r'UGC2162 (Trujillo+17)',]
        mstar = [8.4e7, 3e8, 2.8e7, 2e7]
        # converted from m200c for Amorisco+17 (including 0.15 dex
        # scatter from 
        # Ngc)
        mobs = [9e10, 7e9, 4.5e9, 4.6e9]
        merr = [4e10, 2.4e9, 2.8e9, 8e8]
        robs = [0, 4.6e-3, 8.1e-3, 5e-3]
        z = [zcoma, zcoma, 0.0036, 0.0039]
        # Using (m200,r200) for Amorisco+17 and DF17
        for i in xrange(len(mstar)):
            if robs[i] > 0:
                continue
            robs[i] = conversions.rsph(mobs[i], z[i], ref='200a')
        ls = 'so^v'
        print('indices =', indices)
        for i in xrange(len(label)):
            for ci in c:
                mass = udgmass(mobs[i], robs[i], ci, z[i], name)[0]
                print('{2} mass (c={0}): {1:.2f}'.format(
                    ci, numpy.log10(mass), label[i].split()[0]))
            if mstar[i] < xlim[0]:
                continue
            yerr = (merr[i]/mobs[i]) * mass
            err = ax.errorbar(
                mstar[i], mass, yerr=yerr, fmt=ls[i], mec=udgcolor,
                ecolor=udgcolor, mfc='w', ms=7, mew=2)
            if 'Menclosed' in name:
                curves.append(err)
                curve_labels.append(label[i])
                indices[refs.index('udgs')+n+i+1] = indices.max()+1
            print('refs =', refs)
            print('indices =', indices)

        ## NIHAO simulations (Di Cintio et al. 2017)
        color = '0.6'
        nihao = numpy.loadtxt('literature/Dicintio17_mass_profile_UDGs.txt')
        rnihao = nihao[50:]
        mnihao = nihao[:50]
        # at 30 kpc
        if 'Menclosed' in name:
            jnihao = numpy.argmin(abs(rnihao-30), axis=0)
        else:
            jnihao = -numpy.ones(nihao.shape[1])
        rnihao, mnihao = numpy.transpose(
            [[ri[ji], mi[ji]] for ri, mi, ji
             in zip(rnihao.T, mnihao.T, jnihao)])
        pts, = ax.plot(rnihao, mnihao, 'o', color=color, mec=color, ms=4)
        if 'Menclosed' in name:
            curves.append(pts)
            curve_labels.append('NIHAO (Di Cintio+17)')
            indices[refs.index('udgs')+n+i+1] = indices.max()+1

        ## Remco's lower limit
        mass_lo = udgmass(2e9, 6e-3, c[-1], 0.05, name)[0]
        print('Lower tidal limit (c={0}) = {1:.2e}'.format(c[-1], mass_lo))
        pt = ax.errorbar(
            udg[0], mass_lo, yerr=mass_lo, fmt=',', lolims=[True],
            color='C0', capsize=3, mew=2, elinewidth=2)
        if 'Menclosed' in name:
            curves.append(pt)
            curve_labels.append('Tidal argument (vdBurg+16)')
            indices[refs.index('udgs')+n+i+1] = indices.max()+1
        print('indices =', indices)
    # Milky Way
    print('indices =', indices)
    if 'mw' in ref:
        if 'Menclosed' in name:
            radius = 0.01 * float(name.split('-')[0][-1])
            mass_mw = mass_enclosed(1e12, 5, 0.01, radius, ref='200a')
        else:
            mass_mw = 1e12
        print('MW mass = {0:.2f}e12 Msun'.format(mass_mw/1e12))
        curves.append(
            ax.plot(5e10, mass_mw, '*', mec='0.3', mfc='w', ms=16, mew=2)[0])
        curve_labels.append('Milky Way')
        #ax.annotate('MW', xy=(5e10,6e11), ha='center', va='top', color='0.3',
                    #fontsize=labelsize, fontweight='heavy')
        indices[refs.index('mw')+n] = indices.max()+1

    print('refs =', refs)
    print('curve_labels =', curve_labels)
    print('indices =', indices)
    print('len(curves) =', len(curves))
    curves = [curves[i] if i > -1 else None for i in indices]
    curve_labels = [curve_labels[i] if i > -1 else None for i in indices]
    print('curve_labels =', curve_labels)

    return curves, curve_labels


def udgmass(mradius, radius, c, z, name):
    mass = find_mdelta(mradius, radius, c, z, ref='200a')
    # extrapolate to 10/20/30 kpc only
    if 'Menclosed' in name:
        radius_ref = 1e-3 * 10 * int(name.split('-')[0].split('_')[0][-1])
        mass = mass_enclosed(mass, c, z, radius_ref, ref='200a')
    return mass


def mass_enclosed(mdelta, c, z, radius, ref='200a'):
    """
    Wrapper of kids_ggl_pipeline.halomodel.nfw.mass_enclosed()

    Given an m200, a concentration and a redshift, calculate the mass
    at a fixed radius, given in Mpc

    """
    rho = cosmology.density(z, ref=ref[-1])
    rdelta = conversions.rsph(mdelta, z, ref=ref)
    rs = rdelta / c
    sigma = nfw.delta(c, overdensity=float(ref[:-1])) * rs * rho
    return nfw.mass_enclosed(radius/rs, rs, sigma)


def find_mdelta(mradius, radius, c, z, ref='200a', x0=1e11):
    """
    Loss function to perform the inverse of mass_enclosed(): given
    an overdensity mass (e.g., m200), calculate the mass at a fixed
    physical radius given a concentration.

    In order 

    """
    x = lambda mdelta, mo, radius, c, z, ref: \
        (mo - mass_enclosed(mdelta, c, z, radius, ref=ref))**2
    return optimize.fmin(x, x0, args=(mradius,radius,c,z,ref))


def axlabel(key):
    labels = {'Msat': 'm_{{200}}',
              'Msat_rbg': r'm_\mathrm{{bound}}'}
    for i in xrange(1, 4):
        labels['Menclosed{0}'.format(i)] = \
            r'm(<{0}0\,\mathrm{{kpc}})'.format(i)
    if key in labels:
        return labels[key]
    return key


def get_ylim(name, ref):
    if 'leauthaud12' in ref or 'zu15' in ref:
        ylim = (1e10, 1e14)
    elif name.split('_')[0] in ('Menclosed1', 'Menclosed2'):
        ylim = (5e9, 1e12)
    elif 'Menclosed3' in name:
        ylim = (5e9, 2e13)
    else:
        ylim = (1e10, 1e13)
    return ylim


def parse_args():
    parser = argparse.ArgumentParser()
    add = parser.add_argument
    add('filename')
    args = parser.parse_args()
    return args


main()


