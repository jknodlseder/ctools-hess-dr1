#! /usr/bin/env python
# ==========================================================================
# Show residuals of OFF observations
#
# Copyright (C) 2019 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import sys
import math
import gammalib
import ctools
import cscripts
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import numpy as np


# ======================= #
# Gaussian function class #
# ======================= #
class gaussian(gammalib.GPythonOptimizerFunction):

    # Constructor
    def __init__(self, x_vals, y_vals):

        # Call base class constructor
        gammalib.GPythonOptimizerFunction.__init__(self)

        # Set eval method
        self._set_eval(self.eval)

        # Set data
        self._x_vals = x_vals
        self._y_vals = y_vals

    # Methods
    def eval(self):
        """
        Evaluate function
        """
        # Recover parameters
        pars  = self._pars()
        norm  = pars[0].value()
        mean  = pars[1].value()
        sigma = pars[2].value()

        # Print warnings
        if norm == 0:
            print('Normalisation is zero!')
        if sigma == 0:
            print('Sigma is zero!')

        # Evaluate function values
        y = [norm * math.exp(-0.5*(x-mean)**2/(sigma**2)) for x in self._x_vals]

        # Compute weights (1/sqrt(y))
        weight = []
        for val in self._y_vals:
            if val > 0.0:
                weight.append(1.0/val)
            else:
                weight.append(0.0)

        # Compute Chi Square
        value = 0.0
        for i in range(len(self._x_vals)):
            arg    = self._y_vals[i] - y[i]
            value += arg * arg * weight[i]

        # Evaluate gradient and curvature
        sigma2 = sigma  * sigma
        sigma3 = sigma2 * sigma
        for i in range(len(self._x_vals)):

            # Evaluate function gradients
            dx     = self._x_vals[i] - mean
            dnorm  = y[i]         / norm   * pars[0].scale()
            dmean  = y[i] * dx    / sigma2 * pars[1].scale()
            dsigma = y[i] * dx**2 / sigma3 * pars[2].scale()

            # Setup gradient vector
            arg                 = (self._y_vals[i] - y[i]) * weight[i]
            self.gradient()[0] -= arg * dnorm
            self.gradient()[1] -= arg * dmean
            self.gradient()[2] -= arg * dsigma

            # Setup curvature matrix
            self.curvature()[0,0] +=  dnorm  * dnorm   * weight[i]
            self.curvature()[0,1] +=  dnorm  * dmean   * weight[i]
            self.curvature()[0,2] +=  dnorm  * dsigma  * weight[i]
            self.curvature()[1,0] +=  dmean  * dnorm   * weight[i]
            self.curvature()[1,1] +=  dmean  * dmean   * weight[i]
            self.curvature()[1,2] +=  dmean  * dsigma  * weight[i]
            self.curvature()[2,0] +=  dsigma * dnorm   * weight[i]
            self.curvature()[2,1] +=  dsigma * dmean   * weight[i]
            self.curvature()[2,2] +=  dsigma * dsigma  * weight[i]

        # Set value
        self._set_value(value)

        # Return
        return


# ============================== #
# Plot significance distribution #
# ============================== #
def plot_sigmap(sigmap, ax, nbins=80, sigma_min=-4.0, sigma_max=4.0):
    """
    Plot the significance distribution and return instance of pyplot
    figure and axis.

    Parameters
    ----------
    sigmap : list of floats
        Significance distribution
    ax : pyplot
        Subplot
    nbins : int
        Number of bins in histogram
    sigma_min : float
        Lower limit of the x axis in the plot
    sigma_max : float
        Upper limit of the x axis in the plot
    """
    # Compute bin edges, hence use nbins+1
    binwidth  = (sigma_max - sigma_min) / float(nbins)
    bin_edges = [sigma_min + binwidth*i for i in range(nbins+1)]

    # Draw significance histogram
    y, _, _ = ax.hist(sigmap, bins=bin_edges, histtype='step', color='k')

    # Set initial Gaussian parameters
    y_max = float(y.max())
    if y_max < 1.0:
        y_max = 1.0
    par1  = gammalib.GOptimizerPar('Norm',  y_max)
    par2  = gammalib.GOptimizerPar('Mean',  0.0)
    par3  = gammalib.GOptimizerPar('Sigma', 1.0)
    par1.min(0.1)
    par3.min(0.1)
    pars  = gammalib.GOptimizerPars()
    pars.append(par1)
    pars.append(par2)
    pars.append(par3)

    # Set fit function
    x   = [0.5*(bin_edges[i]+bin_edges[i+1]) for i in range(nbins)]
    fct = gaussian(x, y)

    # Optimize function and compute errors
    opt = gammalib.GOptimizerLM()
    opt.optimize(fct, pars)
    opt.errors(fct, pars)

    # Recover parameters and errors
    norm    = pars[0].value()
    e_norm  = pars[0].error()
    mean    = pars[1].value()
    e_mean  = pars[1].error()
    sigma   = pars[2].value()
    e_sigma = pars[2].error()

    # Draw text box
    msg = 'mean: $%.3f\pm%.3f$\nwidth: $%.3f\pm%.3f$' % \
           (mean, e_mean, sigma, e_sigma)
    ax.text(0.97, 0.95, msg, ha='right', va='top',
            bbox=dict(edgecolor='red', facecolor='white'),
            transform=ax.transAxes)

    # Plot the normal
    yvals = [norm * math.exp(-0.5*(x-mean)**2/(sigma**2)) for x in bin_edges]
    ax.plot(bin_edges, yvals, 'r-')

    # Configure the plot
    ax.set_ylim(0.5, ax.get_ylim()[1]*2.5)
    ax.set_xlim(sigma_min, sigma_max)
    ax.set_xlabel('Significance', fontsize=12)
    ax.set_ylabel('Entries', fontsize=12)
    ax.grid()
    ax.set_yscale('log')

    # Return
    return


# ====================== #
# Plot residual spectrum #
# ====================== #
def plot_resspec(resspec, ax1, ax2, fontsize=11, markersize=4, legend=True):
    """
    Plot residual spectrum

    Parameters
    ----------
    resspec : dict
        Dictionary
    ax1 : pyplot
        First subplot
    ax2 : pyplot
        First subplot
    fontsize : int, optional
        Font size
    markersize : int, optional
        Marker size
    legend : bool, optional
        Plot legend
    """
    # Initialise vectors
    emean      = []
    ebounds    = []
    counts     = []
    e_counts   = []
    model      = []
    residual   = []
    e_residual = []

    # Build vectors
    ebounds.append(resspec[0]['emin'])
    for result in resspec:
        emean.append(result['elogmean'])
        ebounds.append(result['emax'])
        counts.append(result['counts'])
        e_counts.append(math.sqrt(result['counts']))
        model.append(result['model'])
        if result['model'] > 0.0:
            residual.append((result['counts']-result['model'])/math.sqrt(result['model']))
        else:
            residual.append(0.0)
        e_residual.append(1.0)

    # Add model value to be compatible with plt.step
    model = [model[0]] + model

    # Get energy range
    emin = resspec[0]['emin']
    emax = resspec[len(resspec)-1]['emax']

    # Set axes scales
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax1.set_yscale('log')

    # Set energy range
    ax1.set_xlim([emin,emax])
    ax2.set_xlim([emin,emax])

    # Counts and model
    ax1.errorbar(emean, counts, yerr=e_counts, markersize=markersize,
                 fmt='ko', capsize=0, linewidth=1, zorder=2, label='Data')
    ax1.step(ebounds, model, color='r', linewidth=1, zorder=1, label='Model')
    ax1.set_xlabel('Energy (TeV)', fontsize=fontsize)
    ax1.set_ylabel('Counts', fontsize=fontsize)
    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax1.xaxis.set_tick_params(labelsize=fontsize)
    ax1.yaxis.set_tick_params(labelsize=fontsize)

    # Residuals
    ax2.errorbar(emean, residual, yerr=e_residual, markersize=markersize,
                 fmt='ko', capsize=0, linewidth=1, zorder=2)
    ax2.axhline(0, color='r', linestyle='--')
    ax2.set_xlabel('Energy (TeV)', fontsize=fontsize)
    ax2.set_ylabel(r'Residuals ($\sigma$)', fontsize=fontsize)
    ax2.xaxis.set_tick_params(labelsize=fontsize)
    ax2.yaxis.set_tick_params(labelsize=fontsize)

    # Resort labels
    handles, labels = ax1.get_legend_handles_labels()
    order           = [1,0]

    # Add legend
    if legend:
        ax1.legend([handles[idx] for idx in order],
                   [labels[idx] for idx in order],
                   loc='upper right', fontsize=10)

    # Return
    return


# ================= #
# Plot residual map #
# ================= #
def plot_resmap(map, plt1, plt2, steps=10):
    """
    Plot pull histogram

    Parameters
    ----------
    map : str
        Pull filename
    plt1 : pyplot
        First subplot
    plt2 : pyplot
        First subplot
    steps : int
        Number of steps for rebinning
    """
    # Create longitude array from skymap
    x_lon = []
    y_lon = []
    for ix in range(0, map.nx(), steps):
        value = 0.0
        for istep in range(steps):
            if istep+ix < map.nx():
                for iy in range(map.ny()):
                    index  = istep+ix+iy*map.nx()
                    value += map[index]
        x_lon.append(ix)
        y_lon.append(value)

    # Create latitude array from skymap
    x_lat = []
    y_lat = []
    for iy in range(0, map.ny(), steps):
        value = 0.0
        for istep in range(steps):
            if istep+iy < map.ny():
                for ix in range(map.nx()):
                    index  = ix+(istep+iy)*map.nx()
                    value += map[index]
        x_lat.append(iy)
        y_lat.append(value)

    # Plot longitude distribution
    plt1.step(x_lon, y_lon, 'r-', linewidth=2)
    plt1.set_xlabel('Right Ascension (pixels)')
    plt1.set_ylabel('Residual counts')
    plt1.grid(True)

    # Plot latitude distribution
    plt2.step(x_lat, y_lat, 'r-', linewidth=2)
    plt2.set_xlabel('Declination (pixels)')
    plt2.set_ylabel('Residual counts')
    plt2.grid(True)

    # Return
    return


# ==================== #
# Plot sector location #
# ==================== #
def plot_sector_location(ax, index):
    """
    Plot sector location

    Parameters
    ----------
    ax : pyplot
        Subplot
    index : int
        Sector index
    """
    # Set rectangle width and height
    x0     = 0.87
    y0     = 0.88
    width  = 0.04
    height = 0.1

    # Set rectangle locations
    locs = []
    for iy in range(3):
        for ix in range(3):
            x   = x0 + ix * width
            y   = y0 - iy * height
            loc = (x,y)
            locs.append(loc)

    # Plot empty rectangles
    for loc in locs:
        ax.add_patch(Rectangle(loc, width, height, alpha=1, transform=ax.transAxes,
                               fill=None))

    # Plot filled rectangle
    ax.add_patch(Rectangle(locs[index], width, height, alpha=1,
                           transform=ax.transAxes, facecolor='red'))

    # Return
    return


# ======== #
# Plot map #
# ======== #
def plot_map(inmap, sub, smooth=0.0):
    """
    Plot Aitoff map

    Parameters
    ----------
    inmap : `~gammalib.GSkyMap()`
        Input sky map
    sub : pyplot
        Frame for map
    smooth : float, optional
        Map smoothing parameter (degrees)
    """
    # Optionally smooth map
    if smooth > 0.0:
        map = inmap.copy()
        map.smooth('DISK', smooth)
    else:
        map = inmap

    # Create array from skymap
    array = []
    v_max = 0.0
    for iy in range(map.ny()):
        row = []
        for ix in range(map.nx()):
            index = ix+iy*map.nx()
            value = map[index]
            row.append(value)
        array.append(row)

    # Get skymap boundaries
    ra_min  = map.pix2dir(gammalib.GSkyPixel(0.0,map.ny()/2)).ra_deg()
    ra_max  = map.pix2dir(gammalib.GSkyPixel(map.nx(),map.ny()/2)).ra_deg()
    dec_min = map.pix2dir(gammalib.GSkyPixel(map.nx()/2,0)).dec_deg()
    dec_max = map.pix2dir(gammalib.GSkyPixel(map.nx()/2,map.ny())).dec_deg()
    aspect  = abs((ra_max-ra_min)/(dec_max-dec_min))

    # Show Aitoff projection
    c = sub.imshow(array, extent=(ra_min,ra_max,dec_min,dec_max),
                   cmap=plt.get_cmap('jet'), aspect='auto')
    cbar = plt.colorbar(c, orientation='vertical', shrink=0.7)
    sub.grid(True)
    sub.set_xlim([ra_min,ra_max])
    sub.set_ylim([dec_min,dec_max])
    sub.set_xlabel('Right Ascension (deg)', fontsize=12)
    sub.set_ylabel('Declination (deg)', fontsize=12)

    # Return
    return


# =================== #
# Plot radial profile #
# =================== #
def plot_profile(profile, ax1, ax2, fontsize=11, markersize=4, legend=True,
                 rmin=0.0, rmax=2.0):
    """
    Plot radial profile

    Parameters
    ----------
    profile : dict
        Profile dictionary
    ax1 : pyplot
        First subplot
    ax2 : pyplot
        First subplot
    fontsize : int, optional
        Font size
    markersize : int, optional
        Marker size
    legend : bool, optional
        Plot legend
    rmin : float, optional
        Minimum radius (deg)
    rmax : float, optional
        Maximum radius (deg)
    """
    # Quit if nothing to plot
    if len(profile['radius']) == 0:
        return

    # Build vector for residual error
    e_residual = [1.0 for r in profile['radius']]
    
    # Build model vector for step
    dr      = profile['radius'][0] / 0.5
    rbounds = [r-0.5*dr for r in profile['radius']]
    model   = [m for m in profile['model']]
    rbounds = rbounds + [rbounds[len(rbounds)-1]+dr]
    model   = [model[0]] + model

    # Set axes scales
    ax1.set_xscale('linear')
    ax2.set_xscale('linear')
    ax1.set_yscale('linear')

    # Set radius range
    ax1.set_xlim([rmin,rmax])
    ax2.set_xlim([rmin,rmax])

    # Counts and model
    ax1.errorbar(profile['radius'], profile['counts'], yerr=profile['error'],
                 markersize=markersize, fmt='ko', capsize=0, linewidth=1,
                 zorder=2, label='Data')
    ax1.step(rbounds, model, color='r', linewidth=1, zorder=1, label='Model')
    ax1.set_xlabel('Radius (deg)', fontsize=fontsize)
    ax1.set_ylabel('Counts/pixel', fontsize=fontsize)
    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax1.xaxis.set_tick_params(labelsize=fontsize)
    ax1.yaxis.set_tick_params(labelsize=fontsize)

    # Residuals
    ax2.errorbar(profile['radius'], profile['residual'], yerr=e_residual,
                 markersize=markersize, fmt='ko', capsize=0, linewidth=1, zorder=2)
    ax2.axhline(0, color='r', linestyle='--')
    ax2.set_xlabel('Radius (deg)', fontsize=fontsize)
    ax2.set_ylabel(r'Residuals ($\sigma$)', fontsize=fontsize)
    ax2.xaxis.set_tick_params(labelsize=fontsize)
    ax2.yaxis.set_tick_params(labelsize=fontsize)

    # Resort labels
    handles, labels = ax1.get_legend_handles_labels()
    order           = [1,0]

    # Add legend
    if legend:
        ax1.legend([handles[idx] for idx in order],
                   [labels[idx] for idx in order],
                   loc='upper right', fontsize=10)

    # Draw text box with energy range
    msg = '$%.2f - %.2f$ TeV' % (profile['emin'], profile['emax'])
    ax1.text(0.97, 0.95, msg, ha='right', va='top',
             #bbox=dict(edgecolor='red', facecolor='white'),
             transform=ax1.transAxes, fontsize=10)

    # Return
    return


# ====================== #
# Select relevant events #
# ====================== #
def select_events(obs, emin='UNDEF', emax='UNDEF', rad=2.0):
    """
    Select events within a given RoI radius

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    rad : float, optional
        Radial selection radius (deg)

    Returns
    -------
    obs : `~gammalib.GObservations`
        Observation container
    """
    # Setup task parameters
    select = ctools.ctselect(obs)
    select['ra']       = 'UNDEF'
    select['dec']      = 'UNDEF'
    select['rad']      = rad
    select['tmin']     = 'UNDEF'
    select['tmax']     = 'UNDEF'
    select['emin']     = emin
    select['emax']     = emax
    select['usethres'] = 'DEFAULT'

    # Select events
    select.run()

    # Extract observation
    obs = select.obs().copy()

    # Return observation
    return obs


# ==================== #
# Get background model #
# ==================== #
def get_bkg_model(obs, rad=2.0):
    """
    Get background model

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    rad : float, optional
        RoI radius (deg)
    suffix : str, optional
        Background lookup file suffix

    Returns
    -------
    models : `~gammalib.GModels`
        Model container
    """
    # Get energy range for observation
    emin = obs[0].ebounds().emin().TeV()
    emax = obs[0].ebounds().emax().TeV()

    # Setup task parameters
    bkg = cscripts.csbkgmodel(obs)
    bkg['instrument'] = 'HESS'
    bkg['spatial']    = 'LOOKUP'
    bkg['slufile']    = 'bkg_lookup.fits'
    bkg['gradient']   = True
    bkg['spectral']   = 'NODES'
    bkg['ebinalg']    = 'LOG'
    bkg['emin']       = emin
    bkg['emax']       = emax
    bkg['enumbins']   = 8
    bkg['runwise']    = True
    bkg['rad']        = rad
    bkg['debug']      = True

    # Generate background model
    bkg.run()

    # Extract models
    models = bkg.models().copy()

    # Return models
    return models


# ================== #
# Create counts cube #
# ================== #
def create_cntcube(obs, emin=None, emax=None, ebins=20, npix=200, binsz=0.02):
    """
    Create counts cube

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    ebins : int, optional
        Number of energy bins
    npix : int, optional
        Number of spatial pixels in x and y
    binsz : float, optional
        Spatial pixels bin size (deg)

    Returns
    -------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    """
    # Get energy range for observation
    if emin == None:
        emin = obs[0].ebounds().emin().TeV()
    if emax == None:
        emax = obs[0].ebounds().emax().TeV()

    # Setup task parameters
    ctbin = ctools.ctbin(obs)
    ctbin['ebinalg']  = 'LOG'
    ctbin['emin']     = emin
    ctbin['emax']     = emax
    ctbin['enumbins'] = ebins
    ctbin['coordsys'] = 'CEL'
    ctbin['proj']     = 'TAN'
    ctbin['usepnt']   = True
    ctbin['nxpix']    = npix
    ctbin['nypix']    = npix
    ctbin['binsz']    = binsz

    # Generate counts cube
    ctbin.run()

    # Extract counts cube
    cntcube = ctbin.cube().copy()

    # Return counts cube
    return cntcube


# ================= #
# Create model cube #
# ================= #
def create_modcube(obs, cntcube):
    """
    Create model cube

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube

    Returns
    -------
    modcube : `~gammalib.GCTAEventCube`
        Model cube
    """
    # Setup task parameters
    ctmodel = ctools.ctmodel(obs)
    ctmodel.cube(cntcube)

    # Generate model cube
    ctmodel.run()

    # Extract model cube
    modcube = ctmodel.cube().copy()

    # Return model cube
    return modcube


# ============================ #
# Determine spectral residuals #
# ============================ #
def residual_spectral(cntcube, modcube):
    """
    Determine spectral residuals

    Parameters
    ----------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    modcube : `~gammalib.GCTAEventCube`
        Model cube

    Returns
    -------
    residual : dict
        Residual dictionary
    """
    # Extract energy boundaries
    ebounds = cntcube.ebounds().copy()

    # Extract counts and model cubes sky maps
    cnt = cntcube.counts()
    mod = modcube.counts()

    # Initialise result array
    residual = []

    # Loop over energies
    for i in range(ebounds.size()):

        # Determine counts and model sums
        counts = 0.0
        model  = 0.0
        for k in range(cnt.npix()):
            counts += cnt[k,i]
            model  += mod[k,i]

        # Build result dictionary
        result = {'emin':     ebounds.emin(i).TeV(),
                  'emax':     ebounds.emax(i).TeV(),
                  'elogmean': ebounds.elogmean(i).TeV(),
                  'ewidth':   ebounds.ewidth(i).TeV(),
                  'counts':   counts,
                  'model':    model}

        # Append dictionary
        residual.append(result)

    # Return residual
    return residual


# =========================== #
# Determine spatial residuals #
# =========================== #
def residual_spatial(cntcube, modcube, nsample=10, smooth=0.0):
    """
    Determine spatial residuals

    Parameters
    ----------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    modcube : `~gammalib.GCTAEventCube`
        Model cube
    nsample : int, optional
        Sampling for spatial residuals
    smooth : float, optional
        Smoothing for spatial residuals

    Returns
    -------
    map, sigmap : `~gammalib.GSkyMap()`, list
        Residual sky map and array of residuals
    """
    # Extract counts and model cubes sky maps
    cnt = cntcube.counts()
    mod = modcube.counts()

    # Generate stacked maps
    cnt_stacked = cnt.copy()
    mod_stacked = mod.copy()
    cnt_stacked.stack_maps()
    mod_stacked.stack_maps()

    # Optionally smooth stacked maps
    if smooth > 0.0:
        cnt_smoothed = cnt_stacked.copy()
        mod_smoothed = mod_stacked.copy()
        cnt_smoothed.smooth('DISK', smooth)
        mod_smoothed.smooth('DISK', smooth)
        dx            = abs(cnt_smoothed.projection().cdelt(0))
        dy            = abs(cnt_smoothed.projection().cdelt(1))
        renorm        = gammalib.pi * smooth * smooth / (dx*dy)
        cnt_smoothed *= renorm
        mod_smoothed *= renorm
    else:
        cnt_smoothed = cnt_stacked
        mod_smoothed = mod_stacked

    # Generate residual significance map
    map  = cnt_smoothed.copy()
    sign = (cnt_smoothed - mod_smoothed).sign()
    for i in range(map.npix()):
        model_val = mod_smoothed[i]
        if model_val > 1.0e-12:     # Accomodate for smoothed zero's
            data_val = cnt_smoothed[i]
            if data_val > 0.0:
                log_val = math.log(data_val / model_val)
                map[i]  = (data_val * log_val) + model_val - data_val
            else:
                map[i] = model_val
        else:
            map[i] = 0.0
    map *= 2.0
    map  = map.sqrt()
    map *= sign

    # Generate significance histogram
    sigmap = []
    nx     = int(cnt_stacked.nx() / nsample)
    ny     = int(cnt_stacked.ny() / nsample)
    for ix in range(nx):
        for iy in range(ny):
            counts = 0.0
            model  = 0.0
            for kx in range(nsample):
                inx_x = ix*nsample + kx
                for ky in range(nsample):
                    inx_y   = iy*nsample + ky
                    inx     = inx_x + inx_y*map.nx()
                    counts += cnt_stacked[inx]
                    model  += mod_stacked[inx]
            if model >= 5.0:
                sigma = (counts - model) / math.sqrt(model)
                sigmap.append(sigma)

    # Return residual map and sigmap array
    return map, sigmap


# ======================================== #
# Determine spectral residuals for sectors #
# ======================================== #
def residual_spectral_sectors(cntcube, modcube):
    """
    Determine spectral residuals for sectors

    Parameters
    ----------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    modcube : `~gammalib.GCTAEventCube`
        Model cube

    Returns
    -------
    resspec : list of list of dict
        Residual spectra for sectors
    """
    # Log action
    print('Determine spectral residuals for sectors')

    # Extract energy boundaries
    ebounds = cntcube.ebounds().copy()

    # Extract counts and model cubes sky maps
    cnt = cntcube.counts()
    mod = modcube.counts()

    # Initialise result array
    resspec = []

    #
    nx     = 3
    ny     = 3
    nx_pix = int(cnt.nx() / nx)
    ny_pix = int(cnt.ny() / ny)

    # Loop over sectors
    for iy in range(ny):
        iy_min = iy     * ny_pix
        iy_max = iy_min + ny_pix
        for ix in range(nx):
            ix_min = ix     * nx_pix
            ix_max = ix_min + nx_pix

            # Initialise result array
            res = []

            # Loop over energies
            for i in range(ebounds.size()):

                # Determine counts and model sums
                counts = 0.0
                model  = 0.0
                for ix_sub in range(ix_min,ix_max):
                    for iy_sub in range(iy_min,iy_max):
                        inx     = ix_sub + iy_sub * cnt.nx()
                        counts += cnt[inx,i]
                        model  += mod[inx,i]

                # Build result dictionary
                result = {'emin':     ebounds.emin(i).TeV(),
                          'emax':     ebounds.emax(i).TeV(),
                          'elogmean': ebounds.elogmean(i).TeV(),
                          'ewidth':   ebounds.ewidth(i).TeV(),
                          'counts':   counts,
                          'model':    model}

                # Append dictionary
                res.append(result)

            # Append result list to list of residuals
            resspec.append(res)

    # Return residual spectra
    return resspec


# ========================== #
# Determine radial residuals #
# ========================== #
def residual_radial(obs, cntcube, modcube, ebins=4, rad=2.0, dr=0.1):
    """
    Determine radial residuals

    Parameters
    ----------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    modcube : `~gammalib.GCTAEventCube`
        Model cube
    ebins : int, optional
        Number of energy bins
    rad : float, optional
        Maximum radius
    dr : float, optional
        Radius binsize

    Returns
    -------
    residual : dict
        Residual dictionary
    """
    # Extract energy boundaries
    ebounds = cntcube.ebounds().copy()

    # Extract counts and model cubes sky maps
    cnt = cntcube.counts()
    mod = modcube.counts()

    # Define radius axis
    nr     = int(rad/dr+0.5)
    radius = [(r+0.5)*dr for r in range(nr)]

    # Define energy binning
    nebins = ebounds.size()
    debins = float(ebins)/float(nebins)

    # Initialise result arrays
    counts = [[0.0 for r in range(nr)] for e in range(ebins)]
    model  = [[0.0 for r in range(nr)] for e in range(ebins)]
    pixels = [[0.0 for r in range(nr)] for e in range(ebins)]

    # Initialise centre direction
    centre = obs[0].pointing().dir()

    # Loop over cube pixels
    for k in range(cnt.npix()):

        # Compute theta angle (deg)
        theta = centre.dist_deg(cnt.inx2dir(k))

        # Get radius index
        ir = int(theta / dr)

        # Skip pixel if it is not inside radius axis
        if ir < 0 or ir >= nr:
            continue

        # Loop over energies
        for i in range(nebins):

            # Get energy index
            ie = int(float(i)*debins)

            # Skip layer if it is not inside energy bins
            if ie < 0 or ie >= ebins:
                continue

            # Add cube content to array
            counts[ie][ir] += cnt[k,i]
            model[ie][ir]  += mod[k,i]
            pixels[ie][ir] += 1.0

    # Build result dictionary
    residual = []
    for ie in range(ebins):

        # Compute normalised arrays
        rad_array      = []
        counts_array   = []
        error_array    = []
        model_array    = []
        residual_array = []
        for ir in range(nr):
            if pixels[ie][ir] > 0.0:
                if model[ie][ir] > 0.0:
                    res = (counts[ie][ir]-model[ie][ir])/math.sqrt(model[ie][ir])
                else:
                    res = 0.0
                rad_array.append(radius[ir])
                counts_array.append(counts[ie][ir]/pixels[ie][ir])
                error_array.append(math.sqrt(counts[ie][ir])/pixels[ie][ir])
                model_array.append(model[ie][ir]/pixels[ie][ir])
                residual_array.append(res)

        # Build result dictionary
        iemin  = int(ie/debins+0.5)
        iemax  = int((ie+1)/debins-0.5)
        result = {'emin':      ebounds.emin(iemin).TeV(),
                  'emax':      ebounds.emax(iemax).TeV(),
                  'radius':    rad_array,
                  'counts':    counts_array,
                  'error':     error_array,
                  'model':     model_array,
                  'residual':  residual_array}

        # Append dictionary
        residual.append(result)

    # Return residual
    return residual


# ================ #
# Plot observation #
# ================ #
def plot(obs, cntcube, modcube):
    """
    Plot residuals for observation

    Parameters
    ----------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    modcube : `~gammalib.GCTAEventCube`
        Model cube
    """
    # Get spectral residuals
    resspec = residual_spectral(cntcube, modcube)

    # Get spatial residuals
    resmap, sigmap = residual_spatial(cntcube, modcube, smooth=0.2)

    # Initialise figure
    fig = plt.figure(figsize=(14,4))
    fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.15,
                        wspace=0.3)

    # Create new figure for residual spectra
    gs1 = gridspec.GridSpec(3,3)
    gs1.update(hspace=0)
    plt1 = fig.add_subplot(gs1[0:2,0])
    plt2 = fig.add_subplot(gs1[2,0])
    plot_resspec(resspec, plt1, plt2)

    # Create new figure for residual profiles
    plt3 = fig.add_subplot(132)
    plot_sigmap(sigmap, plt3)

    # Create new figure for residual map
    plt4 = fig.add_subplot(133, aspect='equal')
    plot_map(resmap, plt4, smooth=0.0)

    # Build title string
    title = 'Observation ID: %s' % obs[0].id()
    fig.suptitle(title, fontsize=16)

    # Show figure
    #plt.show()

    # Save figure
    filename = 'off_%s_residuals.png' % obs[0].id()
    fig.savefig(filename, dpi=300)

    # Return results
    return


# ============================ #
# Plot observation for sectors #
# ============================ #
def plot_sectors(obs, cntcube, modcube):
    """
    Analysis and plot observation for nine sectors

    Parameters
    ----------
    obs : `~gammalib.GObservation`
        Observation container
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    modcube : `~gammalib.GCTAEventCube`
        Model cube
    """
    # Get spectral residuals
    resspec = residual_spectral_sectors(cntcube, modcube)

    # Initialise figure
    fig = plt.figure(figsize=(10,6))
    fig.subplots_adjust(wspace=0.31)

    # Create panel for first sector
    gs1 = gridspec.GridSpec(9,3)
    gs1.update(left=0.08, right=0.98, top=0.98, bottom=0.28, hspace=0)
    plt1 = plt.subplot(gs1[0:2,0])
    plt2 = plt.subplot(gs1[2,0])
    plot_resspec(resspec[0], plt1, plt2, fontsize=9, legend=False)
    plot_sector_location(plt1, 0)

    # Create panel for second sector
    gs2 = gridspec.GridSpec(9,3)
    gs2.update(left=0.08, right=0.98, top=0.98, bottom=0.28, hspace=0)
    plt3 = fig.add_subplot(gs2[0:2,1])
    plt4 = fig.add_subplot(gs2[2,1])
    plot_resspec(resspec[1], plt3, plt4, fontsize=9, legend=False)
    plot_sector_location(plt3, 1)

    # Create panel for third sector
    gs3 = gridspec.GridSpec(9,3)
    gs3.update(left=0.08, right=0.98, top=0.98, bottom=0.28, hspace=0)
    plt5 = fig.add_subplot(gs3[0:2,2])
    plt6 = fig.add_subplot(gs3[2,2])
    plot_resspec(resspec[2], plt5, plt6, fontsize=9, legend=False)
    plot_sector_location(plt5, 2)

    # Create panel for forth sector
    gs4 = gridspec.GridSpec(9,3)
    gs4.update(left=0.08, right=0.98, top=0.875, bottom=0.175, hspace=0)
    plt7 = plt.subplot(gs4[3:5,0])
    plt8 = plt.subplot(gs4[5,0])
    plot_resspec(resspec[3], plt7, plt8, fontsize=9, legend=False)
    plot_sector_location(plt7, 3)

    # Create panel for fifth sector
    gs5 = gridspec.GridSpec(9,3)
    gs5.update(left=0.08, right=0.98, top=0.875, bottom=0.175, hspace=0)
    plt9  = fig.add_subplot(gs5[3:5,1])
    plt10 = fig.add_subplot(gs5[5,1])
    plot_resspec(resspec[4], plt9, plt10, fontsize=9, legend=False)
    plot_sector_location(plt9, 4)

    # Create panel for sixths sector
    gs6 = gridspec.GridSpec(9,3)
    gs6.update(left=0.08, right=0.98, top=0.875, bottom=0.175, hspace=0)
    plt11 = fig.add_subplot(gs6[3:5,2])
    plt12 = fig.add_subplot(gs6[5,2])
    plot_resspec(resspec[5], plt11, plt12, fontsize=9, legend=False)
    plot_sector_location(plt11, 5)

    # Create panel for seventh sector
    gs7 = gridspec.GridSpec(9,3)
    gs7.update(left=0.08, right=0.98, top=0.77, bottom=0.07, hspace=0)
    plt13 = plt.subplot(gs7[6:8,0])
    plt14 = plt.subplot(gs7[8,0])
    plot_resspec(resspec[6], plt13, plt14, fontsize=9, legend=False)
    plot_sector_location(plt13, 6)

    # Create panel for eight sector
    gs8 = gridspec.GridSpec(9,3)
    gs8.update(left=0.08, right=0.98, top=0.77, bottom=0.07, hspace=0)
    plt15 = plt.subplot(gs8[6:8,1])
    plt16 = plt.subplot(gs8[8,1])
    plot_resspec(resspec[7], plt15, plt16, fontsize=9, legend=False)
    plot_sector_location(plt15, 7)

    # Create panel for ninth sector
    gs9 = gridspec.GridSpec(9,3)
    gs9.update(left=0.08, right=0.98, top=0.77, bottom=0.07, hspace=0)
    plt17 = plt.subplot(gs7[6:8,2])
    plt18 = plt.subplot(gs7[8,2])
    plot_resspec(resspec[8], plt17, plt18, fontsize=9, legend=False)
    plot_sector_location(plt17, 8)

    # Show figure
    #plt.show()

    # Save figure
    filename = 'off_%s_residuals_sectors.png' % obs[0].id()
    fig.savefig(filename, dpi=300)

    # Return results
    return


# ==================== #
# Plot radial profiles #
# ==================== #
def plot_radial_profiles(obs, cntcube, modcube):
    """
    Plot radial profiles

    Parameters
    ----------
    obs : `~gammalib.GObservation`
        Observation container
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    modcube : `~gammalib.GCTAEventCube`
        Model cube
    """
    # Get radial profiles
    profiles = residual_radial(obs, cntcube, modcube, ebins=6)

    # Initialise figure
    fig = plt.figure(figsize=(10,6))
    fig.subplots_adjust(wspace=0.31)

    # Create panel for first energy bin
    gs1 = gridspec.GridSpec(6,3)
    gs1.update(left=0.08, right=0.98, top=0.98, bottom=0.23, hspace=0)
    plt1 = plt.subplot(gs1[0:2,0])
    plt2 = plt.subplot(gs1[2,0])
    plot_profile(profiles[0], plt1, plt2, fontsize=9, legend=False)

    # Create panel for second energy bin
    gs2 = gridspec.GridSpec(6,3)
    gs2.update(left=0.08, right=0.98, top=0.98, bottom=0.23, hspace=0)
    plt3 = fig.add_subplot(gs2[0:2,1])
    plt4 = fig.add_subplot(gs2[2,1])
    plot_profile(profiles[1], plt3, plt4, fontsize=9, legend=False)

    # Create panel for third energy bin
    gs3 = gridspec.GridSpec(6,3)
    gs3.update(left=0.08, right=0.98, top=0.98, bottom=0.23, hspace=0)
    plt5 = fig.add_subplot(gs3[0:2,2])
    plt6 = fig.add_subplot(gs3[2,2])
    plot_profile(profiles[2], plt5, plt6, fontsize=9, legend=False)

    # Create panel for forth energy bin
    gs4 = gridspec.GridSpec(6,3)
    gs4.update(left=0.08, right=0.98, top=0.85, bottom=0.1, hspace=0)
    plt7 = plt.subplot(gs4[3:5,0])
    plt8 = plt.subplot(gs4[5,0])
    plot_profile(profiles[3], plt7, plt8, fontsize=9, legend=False)

    # Create panel for fifth energy bin
    gs5 = gridspec.GridSpec(6,3)
    gs5.update(left=0.08, right=0.98, top=0.85, bottom=0.1, hspace=0)
    plt9  = fig.add_subplot(gs5[3:5,1])
    plt10 = fig.add_subplot(gs5[5,1])
    plot_profile(profiles[4], plt9, plt10, fontsize=9, legend=False)

    # Create panel for sixths energy bin
    gs6 = gridspec.GridSpec(6,3)
    gs6.update(left=0.08, right=0.98, top=0.85, bottom=0.1, hspace=0)
    plt11 = fig.add_subplot(gs6[3:5,2])
    plt12 = fig.add_subplot(gs6[5,2])
    plot_profile(profiles[5], plt11, plt12, fontsize=9, legend=False)

    # Show figure
    #plt.show()

    # Save figure
    filename = 'off_%s_residuals_profiles.png' % obs[0].id()
    fig.savefig(filename, dpi=300)

    # Return results
    return


# ========================== #
# Generate background lookup #
# ========================== #
def generate_background_lookup():
    """
    Generate background lookup from empty field observations
    """
    # Set filenames
    obsname  = '$HESSDATA/obs/obs_off.xml'
    filename = 'bkg_lookup.fits'

    # Continue only if lookup table does not yet exist
    if not os.path.isfile(filename):

        # Initialise lookup table
        emin   = gammalib.GEnergy(0.2,  'TeV')
        emax   = gammalib.GEnergy(50.0, 'TeV')
        ebds   = gammalib.GEbounds(10, emin, emax)
        lookup = gammalib.GCTAModelSpatialLookup(2.0, 0.2, ebds)

        # Load empty field observations
        obs = gammalib.GObservations(obsname)

        # Fill lookup table
        lookup.fill(obs)

        # Save lookup table
        lookup.save(filename, True)

    # Return
    return


# ======================================= #
# Show residuals for all OFF observations #
# ======================================= #
def show_residuals(emin=0.2, emax=20.0):
    """
    Show residuals for all OFF observations

    Parameters
    ----------
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    """
    # Set XML filenames
    obsname = '$HESSDATA/obs/obs_off.xml'

    # Generate background lookup
    generate_background_lookup()

    # Load observations
    inobs = gammalib.GObservations(obsname)

    # Loop over runs in observations
    for run in inobs:

        # Build observation container with single run
        obs = gammalib.GObservations()
        obs.append(run)

        # Select events
        obs = select_events(obs, emin=emin, emax=emax)

        # Generate background model
        models = get_bkg_model(obs)

        # Attach models to observation container
        obs.models(models)

        # Create counts cube
        cntcube = create_cntcube(obs, emin=emin, emax=emax)

        # Create model cube
        modcube = create_modcube(obs, cntcube)

        # Plot observation
        plot(obs, cntcube, modcube)

        # Plot sectors
        plot_sectors(obs, cntcube, modcube)

        # Plot radial profiles
        plot_radial_profiles(obs, cntcube, modcube)

        # Stop after first run
        #break

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show residuals
    show_residuals()
