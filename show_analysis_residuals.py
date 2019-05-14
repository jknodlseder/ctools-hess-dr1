#! /usr/bin/env python
# ==========================================================================
# Show analysis residuals
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
import os
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


# =========================== #
# Plot significance histogram #
# =========================== #
def plot_sighist(sighist, ax, nbins=80, sigma_min=-4.0, sigma_max=4.0, fontsize=11, title=None):
    """
    Plot the significance distribution

    Parameters
    ----------
    sighist : list of floats
        Significance values
    ax : pyplot
        Pyplot panel
    nbins : int, optional
        Number of histogram bins
    sigma_min : float, optional
        Lower limit of the x axis in the plot
    sigma_max : float, optional
        Upper limit of the x axis in the plot
    fontsize : int, optional
        Font size for axes
    title : str, optional
        Title string
    """
    # Compute bin edges, hence use nbins+1
    binwidth  = (sigma_max - sigma_min) / float(nbins)
    bin_edges = [sigma_min + binwidth*i for i in range(nbins+1)]

    # Draw significance histogram
    y, _, _ = ax.hist(sighist, bins=bin_edges, histtype='step', color='k')

    # Set initial Gaussian parameters
    y_max = float(y.max())
    par1  = gammalib.GOptimizerPar('Norm',  y_max)
    par2  = gammalib.GOptimizerPar('Mean',  0.0)
    par3  = gammalib.GOptimizerPar('Sigma', 1.0)
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
            transform=ax.transAxes, fontsize=8)

    # Plot the normal
    yvals = [norm * math.exp(-0.5*(x-mean)**2/(sigma**2)) for x in bin_edges]
    ax.plot(bin_edges, yvals, 'r-')

    # Configure the plot
    ax.set_ylim(0.5, ax.get_ylim()[1]*2.5)
    ax.set_xlim(sigma_min, sigma_max)
    ax.set_xlabel('Significance', fontsize=fontsize)
    ax.set_ylabel('Entries', fontsize=fontsize)
    ax.grid()
    ax.set_yscale('log')

    # Add title
    if title != None:
        ax.set_title(title)

    # Return
    return


# ====================== #
# Plot residual spectrum #
# ====================== #
def plot_resspec(resspec, ax1, ax2, label=None, title=None, legend=True,
                 source=None, markersize=4, fontsize=11):
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
    label : str, optional
        Label string
    title : str, optional
        Title string
    legend : bool, optional
        Plot legend
    source : str, optional
        Source name
    markersize : int, optional
        Marker size
    fontsize : int, optional
        Font size
    """
    # Initialise vectors
    emean      = []
    ebounds    = []
    counts     = []
    e_counts   = []
    model      = []
    bkg        = []
    src        = []
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
        bkg.append(result['bkg'])
        src.append(result['model']-result['bkg'])
        if result['model'] > 0.0:
            if result['counts'] > 0.0:
                res = math.sqrt(2.0*(result['counts']*math.log(result['counts']/result['model'])+
                                result['model']-result['counts']))
            else:
                res = result['model']
            if result['counts'] < result['model']:
                res *= -1.0
            residual.append(res)
        else:
            residual.append(0.0)
        e_residual.append(1.0)

    # Add model value to be compatible with plt.step
    model = [model[0]] + model
    bkg   = [bkg[0]]   + bkg
    src   = [src[0]]   + src

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

    # Set source name
    if source != None:
        srcname = source
    else:
        srcname = 'Source'

    # Counts and model
    ax1.errorbar(emean, counts, yerr=e_counts, markersize=markersize,
                 fmt='ko', capsize=0, linewidth=1, zorder=2, label='Data')
    ax1.step(ebounds, src, color='b', linewidth=1, zorder=1, label=srcname)
    ax1.step(ebounds, bkg, color='g', linewidth=1, zorder=1, label='Background')
    ax1.step(ebounds, model, color='r', linewidth=1, zorder=1, label='Total')
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
    order           = [3,0,1,2]

    # Add legend
    if legend:
        ax1.legend([handles[idx] for idx in order],
                   [labels[idx] for idx in order],
                   loc='lower left', fontsize=7)

    # Add title
    if title != None:
        ax1.set_title(title)

    # Optionally draw label
    if label != None:
        ax1.text(0.03, 0.10, label, ha='left', va='bottom',
                 bbox=dict(edgecolor='red', facecolor='white'),
                 transform=ax1.transAxes, fontsize=10)

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
def plot_map(inmap, sub, smooth=0.0, fontsize=11, title=None):
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
    fontsize : int, optional
        Font size
    title : str, optional
        Title string
    """
    # Optionally smooth map
    if smooth > 0.0:
        map = inmap.copy()
        map.smooth('DISK', smooth)
    else:
        map = inmap

    # Create array from skymap
    array = []
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
                   #cmap=plt.get_cmap('jet'), aspect=aspect)
    cbar = plt.colorbar(c, orientation='vertical', shrink=0.7)
    sub.grid(True)
    sub.set_xlim([ra_min,ra_max])
    sub.set_ylim([dec_min,dec_max])
    sub.set_xlabel('Right Ascension (deg)', fontsize=fontsize)
    sub.set_ylabel('Declination (deg)', fontsize=fontsize)

    # Add title
    if title != None:
        sub.set_title(title)

    # Return
    return


# =================== #
# Plot radial profile #
# =================== #
def plot_profile(profile, ax1, ax2, fontsize=11, markersize=4, legend=True,
                 rmin=0.0, rmax=2.0, textpos='', source=None):
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
    textpos : str, optional
        Text position (e.g. 'bottom')
    source : str, optional
        Source name
    """
    # Build vector for residual error
    e_residual = [1.0 for r in profile['radius']]
    
    # Build model and background vectors for step
    dr         = profile['radius'][0] / 0.5
    rbounds    = [r-0.5*dr for r in profile['radius']]
    model      = [m for m in profile['model']]
    background = [b for b in profile['background']]
    rbounds    = rbounds + [rbounds[len(rbounds)-1]+dr]
    model      = [model[0]] + model
    background = [background[0]] + background

    # Build source vector
    src = [m-profile['background'][i] for i, m in enumerate(profile['model'])]
    src = [src[0]] + src

    # Set source name
    if source != None:
        srcname = source
    else:
        srcname = 'Source'

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
    ax1.step(rbounds, src,        color='b', linewidth=1, zorder=1, label=srcname)
    ax1.step(rbounds, background, color='g', linewidth=1, zorder=1, label='Background')
    ax1.step(rbounds, model,      color='r', linewidth=1, zorder=1, label='Total')
    ax1.set_xlabel('Offset angle (deg)', fontsize=fontsize)
    ax1.set_ylabel('Cts./pix.', fontsize=fontsize)
    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax1.xaxis.set_tick_params(labelsize=fontsize)
    ax1.yaxis.set_tick_params(labelsize=fontsize)
    ax1.get_yaxis().set_label_coords(-0.2,0.5)

    # Residuals
    ax2.errorbar(profile['radius'], profile['residual'], yerr=e_residual,
                 markersize=markersize, fmt='ko', capsize=0, linewidth=1, zorder=2)
    ax2.axhline(0, color='r', linestyle='--')
    ax2.set_xlabel('Offset angle (deg)', fontsize=fontsize)
    ax2.set_ylabel(r'Res. ($\sigma$)', fontsize=fontsize)
    ax2.xaxis.set_tick_params(labelsize=fontsize)
    ax2.yaxis.set_tick_params(labelsize=fontsize)
    ax2.get_yaxis().set_label_coords(-0.2,0.5)

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
    if textpos is 'bottom':
        ax1.text(0.97, 0.18, msg, ha='right', va='top',
                 transform=ax1.transAxes, fontsize=9)
    else:
        ax1.text(0.97, 0.95, msg, ha='right', va='top',
                 transform=ax1.transAxes, fontsize=9)

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
    # Log action
    print('Select events')

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


# ================== #
# Create counts cube #
# ================== #
def create_cntcube(obs, source, emin, emax, xref, yref, nxpix, nypix, ebins=20, binsz=0.02):
    """
    Create counts cube

    Parameters
    ----------
    obs : `~gammalib.GObservation`
        Observation container
    source : str
        Source name
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    xref : float
        Reference Right Ascension value (deg)
    yref : float
        Reference Declination value (deg)
    nxpix : int
        Number of pixels in Right Ascension
    nypix : int
        Number of pixels in Declination
    ebins : int, optional
        Number of energy bins
    binsz : float, optional
        Spatial bin size (deg)

    Returns
    -------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    """
    # Log action
    print('Create counts cube')

    # Set counts cube name
    srcname = source.replace(' ', '_')
    cntname = 'residual_cntcube_%s.fits' % srcname

    # Continue only if selected counts cube does not exist
    if not os.path.isfile(cntname):

        # Setup task parameters
        ctbin = ctools.ctbin(obs)
        ctbin['ebinalg']  = 'LOG'
        ctbin['emin']     = emin
        ctbin['emax']     = emax
        ctbin['enumbins'] = ebins
        if xref == None or xref == None:
            ctbin['usepnt'] = True
        else:
            ctbin['usepnt'] = False
            ctbin['xref']   = xref
            ctbin['yref']   = yref
        ctbin['coordsys'] = 'CEL'
        ctbin['proj']     = 'TAN'
        ctbin['nxpix']    = nxpix
        ctbin['nypix']    = nypix
        ctbin['binsz']    = binsz
        ctbin['outobs']   = cntname

        # Generate counts cube
        ctbin.execute()

    # Load model cube
    cntcube = gammalib.GCTAEventCube(cntname)

    # Return counts cube
    return cntcube


# ================= #
# Create model cube #
# ================= #
def create_modcube(obs, source, cntcube, edisp=True):
    """
    Create model cube

    Parameters
    ----------
    obs : `~gammalib.GObservation`
        Observation container
    source : str
        Source name
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    edisp : bool, optional
        Use energy dispersion

    Returns
    -------
    modcube : `~gammalib.GCTAEventCube`
        Model cube
    """
    # Log action
    print('Create model cube')

    # Set model cube name
    srcname = source.replace(' ', '_')
    modname = 'residual_modcube_%s.fits' % srcname

    # Continue only if selected model cube does not exist
    if not os.path.isfile(modname):

        # Setup task parameters
        ctmodel = ctools.ctmodel(obs)
        ctmodel.cube(cntcube)
        ctmodel['edisp']   = edisp
        ctmodel['outcube'] = modname
        ctmodel['debug']   = True

        # Generate model cube
        ctmodel.execute()

    # Load model cube
    modcube = gammalib.GCTAEventCube(modname)

    # Return model cube
    return modcube


# ====================== #
# Create background cube #
# ====================== #
def create_bkgcube(obs, source, cntcube):
    """
    Create background cube

    Parameters
    ----------
    obs : `~gammalib.GObservation`
        Observation container
    source : str
        Source name
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube

    Returns
    -------
    bkgcube : `~gammalib.GCTAEventCube`
        Background cube
    """
    # Log action
    print('Create background cube')

    # Set model cube name
    srcname = source.replace(' ', '_')
    modname = 'residual_bkgcube_%s.fits' % srcname

    # Continue only if selected model cube does not exist
    if not os.path.isfile(modname):

        # Keep only background model
        models     = obs.models()
        bkg_models = gammalib.GModels()
        for model in models:
            if model.type() == 'CTABackground':
                bkg_models.append(model)
        obs.models(bkg_models)

        # Setup task parameters
        ctmodel = ctools.ctmodel(obs)
        ctmodel.cube(cntcube)
        ctmodel['edisp']   = True
        ctmodel['outcube'] = modname

        # Generate model cube
        ctmodel.execute()

    # Load background cube
    bkgcube = gammalib.GCTAEventCube(modname)

    # Return background cube
    return bkgcube


# ===================================== #
# Create stacked counts and model cubes #
# ===================================== #
def create_stacked_cubes(inobs, source, emin, emax, xref, yref, npix,
                         ebins=20, binsz=0.02, edisp=True):
    """
    Create stacked counts and model cubes

    Parameters
    ----------
    inobs : `~gammalib.GObservations`
        Observation container
    source : str
        Source name
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    xref : float
        Reference Right Ascension value (deg)
    yref : float
        Reference Declination value (deg)
    npix : int
        Number of pixels in Right Ascension and Declination
    ebins : int, optional
        Number of energy bins
    binsz : float, optional
        Spatial bin size (deg)
    edisp : bool, optional
        Use energy dispersion

    Returns
    -------
    cntcube, modcube, bkgcube : tupe of `~gammalib.GCTAEventCube`
        Stacked counts cube, model cube and background cube
    """
    # Log action
    print('Create stacked counts and model cubes')

    # Set counts cube name
    srcname = source.replace(' ', '_')
    cntname = 'residual_stacked_cntcube_%s.fits' % srcname
    modname = 'residual_stacked_modcube_%s.fits' % srcname
    bkgname = 'residual_stacked_bkgcube_%s.fits' % srcname

    # If stacked cubes exist then load them
    if os.path.isfile(cntname) and os.path.isfile(modname) and os.path.isfile(bkgname):
        cntcube_stacked = gammalib.GCTAEventCube(cntname)
        modcube_stacked = gammalib.GCTAEventCube(modname)
        bkgcube_stacked = gammalib.GCTAEventCube(bkgname)

    # ... otherwise generate them
    else:

        # Define counts and model cubes for stacking
        map             = gammalib.GSkyMap('TAN','CEL',0.0,0.0,-binsz,binsz,npix,npix,ebins)
        ebds            = gammalib.GEbounds(ebins,gammalib.GEnergy(emin,'TeV'),gammalib.GEnergy(emax,'TeV'))
        gti             = gammalib.GGti(gammalib.GTime(0.0,'s'),gammalib.GTime(1.0,'s'))
        cntcube_stacked = gammalib.GCTAEventCube(map, ebds, gti)
        modcube_stacked = gammalib.GCTAEventCube(map, ebds, gti)
        bkgcube_stacked = gammalib.GCTAEventCube(map, ebds, gti)

        # Loop over runs in observations
        for run in inobs:

            # Set source name
            source_run = source + '_%s' % run.id()

            # Build observation container with single run
            obs = gammalib.GObservations()
            obs.append(run)
            obs.models(inobs.models())

            # Select events
            obs = select_events(obs, emin=emin, emax=emax)

            # Create counts cube
            cntcube = create_cntcube(obs, source_run, emin, emax, None, None,
                                     npix, npix, ebins=ebins, binsz=binsz)

            # Create model cube
            modcube = create_modcube(obs, source_run, cntcube, edisp=edisp)

            # Create background cube
            bkgcube = create_bkgcube(obs, source_run, cntcube)

            # Stack cubes
            cntcube_stacked = stack_cube(cntcube_stacked, cntcube)
            modcube_stacked = stack_cube(modcube_stacked, modcube)
            bkgcube_stacked = stack_cube(bkgcube_stacked, bkgcube)

            # Stop after first run
            #break

        # Save cubes
        cntcube_stacked.save(cntname, True)
        modcube_stacked.save(modname, True)
        bkgcube_stacked.save(bkgname, True)

    # Return stacked cubes
    return cntcube_stacked, modcube_stacked, bkgcube_stacked


# ============================ #
# Determine spectral residuals #
# ============================ #
def residual_spectral(cntcube, modcube, bkgcube, ra, dec, rad):
    """
    Determine spectral residuals

    Parameters
    ----------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    modcube : `~gammalib.GCTAEventCube`
        Model cube
    bkgcube : `~gammalib.GCTAEventCube`
        Background cube
    ra : float
        Right Ascension of On region centre (deg)
    dec : float
        Declination of On region centre (deg)
    rad : float
        Radius of On region centre (deg)

    Returns
    -------
    res_on : list of dict
        Residual spectrum for On region
    res_fov : list of dict
        Residual spectrum for full FoV region
    """
    # Log action
    print('Determine spectral residuals')

    # Extract energy boundaries
    ebounds = cntcube.ebounds().copy()

    # Extract counts and model cubes sky maps
    cnt = cntcube.counts()
    mod = modcube.counts()
    bkg = bkgcube.counts()

    # Initialise result array
    res_on  = []
    res_fov = []

    # Set
    on = gammalib.GSkyRegionCircle(ra, dec, rad)

    # Loop over energies
    for i in range(ebounds.size()):

        # Determine counts and model sums
        counts_on  = 0.0
        model_on   = 0.0
        bkg_on     = 0.0
        counts_fov = 0.0
        model_fov  = 0.0
        bkg_fov    = 0.0
        for k in range(cnt.npix()):
            counts_fov += cnt[k,i]
            model_fov  += mod[k,i]
            bkg_fov    += bkg[k,i]
            if on.contains(cnt.inx2dir(k)):
                counts_on += cnt[k,i]
                model_on  += mod[k,i]
                bkg_on    += bkg[k,i]

        # Build result dictionaries
        result_on  = {'emin':     ebounds.emin(i).TeV(),
                      'emax':     ebounds.emax(i).TeV(),
                      'elogmean': ebounds.elogmean(i).TeV(),
                      'ewidth':   ebounds.ewidth(i).TeV(),
                      'counts':   counts_on,
                      'model':    model_on,
                      'bkg':      bkg_on}
        result_fov = {'emin':     ebounds.emin(i).TeV(),
                      'emax':     ebounds.emax(i).TeV(),
                      'elogmean': ebounds.elogmean(i).TeV(),
                      'ewidth':   ebounds.ewidth(i).TeV(),
                      'counts':   counts_fov,
                      'model':    model_fov,
                      'bkg':      bkg_fov}

        # Append dictionaries
        res_on.append(result_on)
        res_fov.append(result_fov)

    # Return residual arrays
    return res_on, res_fov


# ======================================== #
# Determine spectral residuals for sectors #
# ======================================== #
def residual_spectral_sectors(cntcube, modcube, bkgcube):
    """
    Determine spectral residuals for sectors

    Parameters
    ----------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    modcube : `~gammalib.GCTAEventCube`
        Model cube
    bkgcube : `~gammalib.GCTAEventCube`
        Background cube

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
    bkg = bkgcube.counts()

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
                back   = 0.0
                for ix_sub in range(ix_min,ix_max):
                    for iy_sub in range(iy_min,iy_max):
                        inx     = ix_sub + iy_sub * cnt.nx()
                        counts += cnt[inx,i]
                        model  += mod[inx,i]
                        back   += bkg[inx,i]

                # Build result dictionary
                result = {'emin':     ebounds.emin(i).TeV(),
                          'emax':     ebounds.emax(i).TeV(),
                          'elogmean': ebounds.elogmean(i).TeV(),
                          'ewidth':   ebounds.ewidth(i).TeV(),
                          'counts':   counts,
                          'model':    model,
                          'bkg':      back}

                # Append dictionary
                res.append(result)

            # Append result list to list of residuals
            resspec.append(res)

    # Return residual spectra
    return resspec


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
    map : `~gammalib.GSkyMap()`
        Residual sky map
    """
    # Log action
    print('Determine spatial residuals')

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


# ========================== #
# Determine radial residuals #
# ========================== #
def residual_radial(cntcube, modcube, bkgcube, ebins=4, rad=2.0, dr=0.1):
    """
    Determine radial residuals

    Parameters
    ----------
    cntcube : `~gammalib.GCTAEventCube`
        Counts cube
    modcube : `~gammalib.GCTAEventCube`
        Model cube
    bkgcube : `~gammalib.GCTAEventCube`
        Background cube
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
    bkg = bkgcube.counts()

    # Define radius axis
    nr     = int(rad/dr+0.5)
    radius = [(r+0.5)*dr for r in range(nr)]

    # Define energy binning
    nebins = ebounds.size()
    debins = float(ebins)/float(nebins)

    # Initialise result arrays
    counts      = [[0.0 for r in range(nr)] for e in range(ebins)]
    model       = [[0.0 for r in range(nr)] for e in range(ebins)]
    background  = [[0.0 for r in range(nr)] for e in range(ebins)]
    pixels      = [[0.0 for r in range(nr)] for e in range(ebins)]

    # Initialise centre direction
    centre = gammalib.GSkyDir()
    centre.radec_deg(0.0, 0.0)

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
            counts[ie][ir]      += cnt[k,i]
            model[ie][ir]       += mod[k,i]
            background[ie][ir]  += bkg[k,i]
            pixels[ie][ir]      += 1.0

    # Build result dictionary
    residual = []
    for ie in range(ebins):

        # Compute normalised arrays
        rad_array        = []
        counts_array     = []
        error_array      = []
        model_array      = []
        background_array = []
        residual_array   = []
        for ir in range(nr):
            if pixels[ie][ir] > 0.0:
                if model[ie][ir] > 0.0:
                    if counts[ie][ir] > 0.0:
                        res = math.sqrt(2.0*(counts[ie][ir]*math.log(counts[ie][ir]/model[ie][ir])+
                                model[ie][ir]-counts[ie][ir]))
                    else:
                        res = model[ie][ir]
                    if counts[ie][ir] < model[ie][ir]:
                        res *= -1.0
                else:
                    res = 0.0
                rad_array.append(radius[ir])
                counts_array.append(counts[ie][ir]/pixels[ie][ir])
                error_array.append(math.sqrt(counts[ie][ir])/pixels[ie][ir])
                model_array.append(model[ie][ir]/pixels[ie][ir])
                background_array.append(background[ie][ir]/pixels[ie][ir])
                residual_array.append(res)

        # Build result dictionary
        iemin  = int(ie/debins+0.5)
        iemax  = int((ie+1)/debins-0.5)
        result = {'emin':       ebounds.emin(iemin).TeV(),
                  'emax':       ebounds.emax(iemax).TeV(),
                  'radius':     rad_array,
                  'counts':     counts_array,
                  'error':      error_array,
                  'model':      model_array,
                  'background': background_array,
                  'residual':   residual_array}

        # Append dictionary
        residual.append(result)

    # Return residual
    return residual


# ========== #
# Stack cube #
# ========== #
def stack_cube(cube_stacked, cube):
    """
    Stack cube

    Parameters
    ----------
    cube_stacked : `~gammalib.GCTAEventCube`
        Stacked cube
    cube : `~gammalib.GCTAEventCube`
        Cube

    Returns
    -------
    cube_stacked : `~gammalib.GCTAEventCube`
        Stacked cube
    """
    # Extract sky maps
    map     = cube.counts()
    stacked = cube_stacked.counts()

    # Copy content
    for k in range(map.nmaps()):
        for i in range(map.npix()):
            stacked[i,k] += map[i,k]

    # Put back skymap
    cube_stacked.counts(stacked)

    # Return stacked map
    return cube_stacked


# ============================= #
# Analyse and plot observations #
# ============================= #
def analyse_and_plot(source, obs, emin, emax, xref, yref, nxpix, nypix, edisp=True):
    """
    Analysis and plot observations

    Parameters
    ----------
    source : str
        Source name
    obs : `~gammalib.GObservation`
        Observation container
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    xref : float
        Reference Right Ascension value (deg)
    yref : float
        Reference Declination value (deg)
    nxpix : int
        Number of pixels in Right Ascension
    nypix : int
        Number of pixels in Declination
    edisp : bool, optional
        Use energy dispersion
    """
    # Select events
    obs = select_events(obs)

    # Create counts cube
    cntcube = create_cntcube(obs,  source, emin, emax, xref, yref, nxpix, nypix)

    # Create model cube
    modcube = create_modcube(obs, source, cntcube, edisp=edisp)

    # Get spectral residuals
    resspec = residual_spectral(cntcube, modcube)

    # Get spatial residuals
    resmap, sigmap = residual_spatial(cntcube, modcube)

    # Initialise figure
    fig = plt.figure(figsize=(14,4))
    fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.15,
                        wspace=0.3)

    # Create new figure for residual spectra
    plt1 = fig.add_subplot(231)
    plt2 = fig.add_subplot(234)
    fig.subplots_adjust(hspace=0)
    fig.add_axes(sharex=True)
    plot_resspec(resspec, plt1, plt2)

    # Create new figure for residual histogram
    plt3 = fig.add_subplot(132)
    plot_sighist(sigmap, plt3)

    # Create new figure for residual map
    plt4 = fig.add_subplot(133, aspect='equal')
    plot_map(resmap, plt4, smooth=0.2)

    # Build title string
    title = '%s observations' % source
    fig.suptitle(title, fontsize=16)

    # Show figure
    plt.show()

    # Save figure
    filename = 'residuals_%s.png' % source
    fig.savefig(filename, dpi=300)

    # Return results
    return


# ============================= #
# Analyse and plot observations #
# ============================= #
def analyse_and_plot_large(source, obs, emin, emax, xref, yref, nxpix, nypix,
                           ra, dec, rad, edisp=True):
    """
    Analysis and plot observations (2nd variant)

    Parameters
    ----------
    source : str
        Source name
    obs : `~gammalib.GObservation`
        Observation container
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    xref : float
        Reference Right Ascension value (deg)
    yref : float
        Reference Declination value (deg)
    nxpix : int
        Number of pixels in Right Ascension
    nypix : int
        Number of pixels in Declination
    ra : float
        Right Ascension of On region centre (deg)
    dec : float
        Declination of On region centre (deg)
    rad : float
        Radius of On region centre (deg)
    edisp : bool, optional
        Use energy dispersion
    """
    # Select events
    obs = select_events(obs)

    # Create counts cube
    cntcube = create_cntcube(obs,  source, emin, emax, xref, yref, nxpix, nypix)

    # Create model cube
    modcube = create_modcube(obs, source, cntcube, edisp=edisp)

    # Create background cube
    bkgcube = create_bkgcube(obs, source, cntcube)

    # Get spectral residuals
    res_on, res_fov = residual_spectral(cntcube, modcube, bkgcube, ra, dec, rad)

    # Get spatial residuals
    resmap, sigmap = residual_spatial(cntcube, modcube, smooth=0.2)

    # Initialise figure
    fig = plt.figure(figsize=(9,6))

    # Create new figure for residual spectra
    gs1 = gridspec.GridSpec(4,2)
    gs1.update(left=0.08, right=0.95, hspace=0, top=0.95, bottom=0.20)
    plt1 = fig.add_subplot(gs1[0,0])
    plt2 = fig.add_subplot(gs1[1,0])
    fig.add_axes(sharex=True)
    plot_resspec(res_on, plt1, plt2, title='On region')

    # Create new figure for residual spectra
    gs2 = gridspec.GridSpec(4,2)
    gs2.update(left=0.11, right=0.98, hspace=0, top=0.95, bottom=0.20)
    plt3 = fig.add_subplot(gs2[0,1])
    plt4 = fig.add_subplot(gs2[1,1])
    fig.add_axes(sharex=True)
    plot_resspec(res_fov, plt3, plt4, title='Full field of view', legend=False)

    # Create new figure for residual histogram
    gs3  = gridspec.GridSpec(2,2)
    gs3.update(left=0.08, right=1.10, top=0.90, bottom=0.10)
    plt5 = fig.add_subplot(gs3[1,0])
    plot_sighist(sigmap, plt5)

    # Create new figure for residual map
    gs4  = gridspec.GridSpec(2,2)
    gs4.update(left=0.30, right=0.98, top=0.90, bottom=0.10)
    plt6 = fig.add_subplot(gs4[1,1], aspect='equal')
    plot_map(resmap, plt6, smooth=0.0)

    # Build title string
    #title = '%s observations' % source
    #fig.suptitle(title, fontsize=16)

    # Show figure
    plt.show()

    # Save figure
    filename = 'residuals_%s.png' % source
    fig.savefig(filename, dpi=300)

    # Return results
    return


# ============================= #
# Analyse and plot observations #
# ============================= #
def analyse_and_plot_small(source, obs, emin, emax, xref, yref, nxpix, nypix,
                           ra, dec, rad, edisp=True):
    """
    Analysis and plot observations (small variant)

    Parameters
    ----------
    source : str
        Source name
    obs : `~gammalib.GObservation`
        Observation container
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    xref : float
        Reference Right Ascension value (deg)
    yref : float
        Reference Declination value (deg)
    nxpix : int
        Number of pixels in Right Ascension
    nypix : int
        Number of pixels in Declination
    ra : float
        On region centre Right Ascension (deg)
    dec : float
        On region centre Declination (deg)
    rad : float
        On region radius (deg)
    edisp : bool, optional
        Use energy dispersion
    """
    # Select events
    obs = select_events(obs)

    # Create counts cube
    cntcube = create_cntcube(obs,  source, emin, emax, xref, yref, nxpix, nypix)

    # Create model cube
    modcube = create_modcube(obs, source, cntcube, edisp=edisp)

    # Create background cube
    bkgcube = create_bkgcube(obs, source, cntcube)

    # Get spectral residuals
    res_on, res_fov = residual_spectral(cntcube, modcube, bkgcube, ra, dec, rad)

    # Get spatial residuals
    resmap, sigmap = residual_spatial(cntcube, modcube, smooth=0.2)

    # Initialise figure
    fig = plt.figure(figsize=(14,3.2))
    fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.15,
                        wspace=0.3)

    # Create panel for On region residual spectra
    gs1 = gridspec.GridSpec(3,4)
    gs1.update(hspace=0)
    plt1 = fig.add_subplot(gs1[0:2,0])
    plt2 = fig.add_subplot(gs1[2,0])
    fig.add_axes(sharex=True)
    plot_resspec(res_on, plt1, plt2, source=source, title='On region spectral residuals')

    # Create panel for FoV region residual spectra
    gs2 = gridspec.GridSpec(3,4)
    gs2.update(hspace=0)
    plt3 = fig.add_subplot(gs2[0:2,1])
    plt4 = fig.add_subplot(gs2[2,1])
    fig.add_axes(sharex=True)
    plot_resspec(res_fov, plt3, plt4, source=source, title='FoV spectral residuals')

    # Create panel for residual histogram
    plt5 = fig.add_subplot(143)
    plot_sighist(sigmap, plt5, title='Spatial residuals')

    # Create new figure for residual map
    gs4  = gridspec.GridSpec(1,4)
    gs4.update(left=0.0, right=1.0, top=0.90, bottom=0.15)
    plt6 = fig.add_subplot(gs4[0,3], aspect='equal')
    plot_map(resmap, plt6, smooth=0.0, title='Residual significance map')

    # Show figure
    plt.show()

    # Save figure
    filename = 'residuals_%s.png' % source
    fig.savefig(filename, dpi=300)

    # Return results
    return


# ========================== #
# Analyse and plot residuals #
# ========================== #
def analyse_and_plot_sectors(source, obs, emin, emax, xref, yref, nxpix, nypix,
                             edisp=True):
    """
    Analysis and plot residuals (nine sectors)

    Parameters
    ----------
    source : str
        Source name
    obs : `~gammalib.GObservations`
        Observation container
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    xref : float
        Reference Right Ascension value (deg)
    yref : float
        Reference Declination value (deg)
    nxpix : int
        Number of pixels in Right Ascension
    nypix : int
        Number of pixels in Declination
    edisp : bool, optional
        Use energy dispersion
    """
    # Select events
    obs = select_events(obs)

    # Create counts cube
    cntcube = create_cntcube(obs,  source, emin, emax, xref, yref, nxpix, nypix)

    # Create model cube
    modcube = create_modcube(obs, source, cntcube, edisp=edisp)

    # Create background cube
    bkgcube = create_bkgcube(obs, source, cntcube)

    # Get spectral residuals
    resspec = residual_spectral_sectors(cntcube, modcube, bkgcube)

    # Initialise figure
    fig = plt.figure(figsize=(10,6))
    fig.subplots_adjust(wspace=0.31)

    # Create panel for first sector
    gs1 = gridspec.GridSpec(9,3)
    gs1.update(left=0.08, right=0.98, top=0.98, bottom=0.28, hspace=0)
    plt1 = plt.subplot(gs1[0:2,0])
    plt2 = plt.subplot(gs1[2,0])
    plot_resspec(resspec[0], plt1, plt2, source=source, fontsize=9, legend=False)
    plot_sector_location(plt1, 0)

    # Create panel for second sector
    gs2 = gridspec.GridSpec(9,3)
    gs2.update(left=0.08, right=0.98, top=0.98, bottom=0.28, hspace=0)
    plt3 = fig.add_subplot(gs2[0:2,1])
    plt4 = fig.add_subplot(gs2[2,1])
    plot_resspec(resspec[1], plt3, plt4, source=source, fontsize=9, legend=False)
    plot_sector_location(plt3, 1)

    # Create panel for third sector
    gs3 = gridspec.GridSpec(9,3)
    gs3.update(left=0.08, right=0.98, top=0.98, bottom=0.28, hspace=0)
    plt5 = fig.add_subplot(gs3[0:2,2])
    plt6 = fig.add_subplot(gs3[2,2])
    plot_resspec(resspec[2], plt5, plt6, source=source, fontsize=9, legend=False)
    plot_sector_location(plt5, 2)

    # Create panel for forth sector
    gs4 = gridspec.GridSpec(9,3)
    gs4.update(left=0.08, right=0.98, top=0.875, bottom=0.175, hspace=0)
    plt7 = plt.subplot(gs4[3:5,0])
    plt8 = plt.subplot(gs4[5,0])
    plot_resspec(resspec[3], plt7, plt8, source=source, fontsize=9, legend=False)
    plot_sector_location(plt7, 3)

    # Create panel for fifth sector
    gs5 = gridspec.GridSpec(9,3)
    gs5.update(left=0.08, right=0.98, top=0.875, bottom=0.175, hspace=0)
    plt9  = fig.add_subplot(gs5[3:5,1])
    plt10 = fig.add_subplot(gs5[5,1])
    plot_resspec(resspec[4], plt9, plt10, source=source, fontsize=9, legend=False)
    plot_sector_location(plt9, 4)

    # Create panel for sixths sector
    gs6 = gridspec.GridSpec(9,3)
    gs6.update(left=0.08, right=0.98, top=0.875, bottom=0.175, hspace=0)
    plt11 = fig.add_subplot(gs6[3:5,2])
    plt12 = fig.add_subplot(gs6[5,2])
    plot_resspec(resspec[5], plt11, plt12, source=source, fontsize=9, legend=False)
    plot_sector_location(plt11, 5)

    # Create panel for seventh sector
    gs7 = gridspec.GridSpec(9,3)
    gs7.update(left=0.08, right=0.98, top=0.77, bottom=0.07, hspace=0)
    plt13 = plt.subplot(gs7[6:8,0])
    plt14 = plt.subplot(gs7[8,0])
    plot_resspec(resspec[6], plt13, plt14, source=source, fontsize=9, legend=False)
    plot_sector_location(plt13, 6)

    # Create panel for eight sector
    gs8 = gridspec.GridSpec(9,3)
    gs8.update(left=0.08, right=0.98, top=0.77, bottom=0.07, hspace=0)
    plt15 = plt.subplot(gs8[6:8,1])
    plt16 = plt.subplot(gs8[8,1])
    plot_resspec(resspec[7], plt15, plt16, source=source, fontsize=9, legend=False)
    plot_sector_location(plt15, 7)

    # Create panel for ninth sector
    gs9 = gridspec.GridSpec(9,3)
    gs9.update(left=0.08, right=0.98, top=0.77, bottom=0.07, hspace=0)
    plt17 = plt.subplot(gs7[6:8,2])
    plt18 = plt.subplot(gs7[8,2])
    plot_resspec(resspec[8], plt17, plt18, source=source, fontsize=9, legend=False)
    plot_sector_location(plt17, 8)

    # Show figure
    plt.show()

    # Save figure
    filename = 'residuals_sectors_%s.png' % source
    fig.savefig(filename, dpi=300)

    # Return results
    return


# ========================== #
# Analyse and plot residuals #
# ========================== #
def analyse_and_plot_profiles(source, obs, emin, emax, xref, yref, npix,
                              edisp=True):
    """
    Analysis and plot radial profiles

    Parameters
    ----------
    source : str
        Source name
    obs : `~gammalib.GObservation`
        Observation container
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    xref : float
        Reference Right Ascension value (deg)
    yref : float
        Reference Declination value (deg)
    npix : int
        Number of pixels in Right Ascension and Declination (deg)
    edisp : bool, optional
        Use energy dispersion
    """
    # Get stacked cubes
    cntcube, modcube, bkgcube = create_stacked_cubes(obs, source, emin, emax,
                                                     xref, yref, npix,
                                                     ebins=20, binsz=0.02,
                                                     edisp=edisp)

    # Get radial profiles
    profiles = residual_radial(cntcube, modcube, bkgcube, ebins=6)

    # Initialise figure
    fig = plt.figure(figsize=(10,3.5))
    fig.subplots_adjust(wspace=0.31)

    # Create panel for first energy bin
    gs1 = gridspec.GridSpec(6,3)
    gs1.update(left=0.08, right=0.98, top=0.98, bottom=0.26, hspace=0)
    plt1 = plt.subplot(gs1[0:2,0])
    plt2 = plt.subplot(gs1[2,0])
    plot_profile(profiles[0], plt1, plt2, fontsize=9, legend=False, source=source)

    # Create panel for second energy bin
    gs2 = gridspec.GridSpec(6,3)
    gs2.update(left=0.08, right=0.98, top=0.98, bottom=0.26, hspace=0)
    plt3 = fig.add_subplot(gs2[0:2,1])
    plt4 = fig.add_subplot(gs2[2,1])
    plot_profile(profiles[1], plt3, plt4, fontsize=9, legend=False, source=source)

    # Create panel for third energy bin
    gs3 = gridspec.GridSpec(6,3)
    gs3.update(left=0.08, right=0.98, top=0.98, bottom=0.26, hspace=0)
    plt5 = fig.add_subplot(gs3[0:2,2])
    plt6 = fig.add_subplot(gs3[2,2])
    plot_profile(profiles[2], plt5, plt6, fontsize=9, legend=False, source=source)

    # Create panel for forth energy bin
    gs4 = gridspec.GridSpec(6,3)
    gs4.update(left=0.08, right=0.98, top=0.84, bottom=0.12, hspace=0)
    plt7 = plt.subplot(gs4[3:5,0])
    plt8 = plt.subplot(gs4[5,0])
    plot_profile(profiles[3], plt7, plt8, fontsize=9, legend=False, source=source)

    # Create panel for fifth energy bin
    gs5 = gridspec.GridSpec(6,3)
    gs5.update(left=0.08, right=0.98, top=0.84, bottom=0.12, hspace=0)
    plt9  = fig.add_subplot(gs5[3:5,1])
    plt10 = fig.add_subplot(gs5[5,1])
    plot_profile(profiles[4], plt9, plt10, fontsize=9, legend=False, textpos='bottom', source=source)

    # Create panel for sixths energy bin
    gs6 = gridspec.GridSpec(6,3)
    gs6.update(left=0.08, right=0.98, top=0.84, bottom=0.12, hspace=0)
    plt11 = fig.add_subplot(gs6[3:5,2])
    plt12 = fig.add_subplot(gs6[5,2])
    plot_profile(profiles[5], plt11, plt12, fontsize=9, legend=False, textpos='bottom', source=source)

    # Show figure
    plt.show()

    # Save figure
    filename = 'residuals_profiles_%s.png' % source
    fig.savefig(filename, dpi=300)

    # Return
    return


# ======================= #
# Show analysis residuals #
# ======================= #
def show_analysis_residuals(source, obsfile, resfile, emin, emax,
                            xref, yref, nxpix, nypix,
                            ra, dec, rad, edisp=True):
    """
    Show analysis residuals

    Parameters
    ----------
    source : str
        Source name
    obsfile : str
        Result observation definition XML file
    resfile : str
        Result model definition XML file
    emin : float
        Minimum energy (TeV)
    emax : float
        Maximum energy (TeV)
    xref : float
        Reference Right Ascension value (deg)
    yref : float
        Reference Declination value (deg)
    nxpix : int
        Number of pixels in Right Ascension
    nypix : int
        Number of pixels in Declination
    ra : float
        Right Ascension of On region centre (deg)
    dec : float
        Declination of On region centre (deg)
    rad : float
        Radius of On region centre (deg)
    edisp : bool, optional
        Use energy dispersion
    """
    # Load observations and models
    obs    = gammalib.GObservations(obsfile)
    models = gammalib.GModels(resfile)

    # Attach models to observations
    obs.models(models)

    # Analyse and plot observation
    analyse_and_plot_small(source, obs, emin, emax, xref, yref, nxpix, nypix,
                           ra, dec, rad, edisp=edisp)

    # Analyse and plot sectors
    analyse_and_plot_sectors(source, obs, emin, emax, xref, yref, nxpix, nypix,
                             edisp=edisp)

    # Analyse and plot radial profiles
    analyse_and_plot_profiles(source, obs, emin, emax, xref, yref, 200,
                              edisp=edisp)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Set usage string
    usage = 'show_analysis_residuals.py [source]'

    # Set default options
    options = []

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Set source dependent parameters
    if args[0] == 'Crab':
        source  = 'Crab'
        obsfile = 'obs_crab_selected.xml'
        resfile = 'crab_results_gauss_eplaw_lookup_grad_hess_edisp.xml'
        emin    =  0.67
        emax    = 30.0
        xref    = 83.63
        yref    = 22.01
        nxpix   = 350
        nypix   = 250
        ra      = 83.6330
        dec     = 22.0145
        rad     = 0.2
        edisp   = True
    elif args[0] == 'MSH-PL':
        source  = 'MSH 15-52 plaw'
        obsfile = 'obs_msh_selected.xml'
        resfile = 'msh_results_egauss_plaw_lookup_grad_hess_edisp.xml'
        emin    =  0.381
        emax    = 40.0
        xref    = 228.4817
        yref    = -59.1358
        nxpix   = 200
        nypix   = 250
        ra      = 228.547
        dec     = -59.174
        rad     = 0.2
        edisp   = False
    elif args[0] == 'MSH':
        source  = 'MSH 15-52'
        obsfile = 'obs_msh_selected.xml'
        resfile = 'msh_results_egauss_eplaw_lookup_grad_hess_edisp.xml'
        emin    =  0.381
        emax    = 40.0
        xref    = 228.4817
        yref    = -59.1358
        nxpix   = 200
        nypix   = 250
        ra      = 228.547
        dec     = -59.174
        rad     = 0.2
        edisp   = False
    elif args[0] == 'RX':
        source  = 'RX J1713.7-3946'
        obsfile = 'obs_rx_selected.xml'
        resfile = 'rx_results_map0.48_eplaw_lookup_grad_hess_edisp.xml'
        emin    =  0.3
        emax    = 50.0
        xref    = 258.1125
        yref    = -39.6867
        nxpix   = 300
        nypix   = 350
        ra      = 258.1125
        dec     = -39.6867
        rad     = 0.6
        edisp   = False
    elif args[0] == 'PKS':
        source  = 'PKS 2155-304'
        obsfile = 'obs_pks_t200_selected.xml'
        resfile = 'pks_t200_results_ptsrc_logparabola_lookup_grad_hess_edisp.xml'
        emin    =  0.2
        emax    = 10.0
        xref    = 329.71694
        yref    = -30.22559
        nxpix   = 250
        nypix   = 250
        ra      = 329.71694
        dec     = -30.22559
        rad     = 0.2
        edisp   = True
    else:
        obsfile = None
        resfile = None
        print('Source "'+args[0]+'" not valid. Select one of "Crab", "MSH", "RX" or "PKS"')

    # Show analysis residuals
    if obsfile != None:
        show_analysis_residuals(source, obsfile, resfile, emin, emax,
                                xref, yref, nxpix, nypix,
                                ra, dec, rad, edisp=edisp)
