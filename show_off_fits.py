#! /usr/bin/env python
# ==========================================================================
# Show OFF observations model fit results
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
import pickle
import gammalib
import ctools
import cscripts
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# ======================= #
# Get observation results #
# ======================= #
def get_obs_results(obsname, resname):
    """
    Get observation results

    Parameters
    ----------
    obsname : str
        Observation definition XML file
    resname : str
        Model result XML file
    """
    # Initialise results
    results = []

    # Load observations and models
    inobs  = gammalib.GObservations(obsname)
    models = gammalib.GModels(resname)

    # Loop over runs in observations
    for run in inobs:

        # Build observation container with single run
        obs = gammalib.GObservations()
        obs.append(run)

        # Get background model
        has_model = False
        for model in models:
            if model.classname() == 'GCTAModelBackground' and model.is_valid('HESS', run.id()):
                has_model = True
                break

        # Fall through if we have no model
        if not has_model:
            continue

        # Set results
        result = {'id': run.id(),
                  'ontime': run.ontime(), 'livetime': run.livetime(),
                  'deadc': run.deadc(),
                  'ra': run.pointing().dir().ra_deg(),
                  'dec': run.pointing().dir().dec_deg(),
                  'zenith': run.pointing().zenith(),
                  'azimuth': run.pointing().azimuth(),
                  'detx': model.spatial()[1].value(),
                  'e_detx': model.spatial()[1].error(),
                  'dety': model.spatial()[2].value(),
                  'e_dety': model.spatial()[2].error()}

        # Append results
        results.append(result)

    # Return
    return results


# ======================== #
# Plot observation results #
# ======================== #
def plot_obs_results(ax, obsname, resname, color, srcname):
    """
    Plot observation results

    Parameters
    ----------
    ax : pyplot
        Subplot
    obsname : str
        Observation definition XML file
    resname : str
        Model result XML file
    color : str
        Color string
    srcname : str
        Source name
    """
    # Get results
    results = get_obs_results(obsname, resname)

    # Generate arrays
    zenith  = [res['zenith'] for res in results]
    detx    = [abs(res['detx']) for res in results]
    e_detx  = [res['e_detx'] for res in results]
    dety    = [abs(res['dety']) for res in results]
    e_dety  = [res['e_dety'] for res in results]
    detxy   = [math.sqrt(res['detx']*res['detx']+res['dety']*res['dety']) for res in results]
    e_detxy = [math.sqrt(res['e_detx']*res['e_detx']*res['detx']*res['detx'] /
                         (res['detx']*res['detx']+res['dety']*res['dety']) +
                         res['e_dety']*res['e_dety']*res['dety']*res['dety'] /
                         (res['detx']*res['detx']+res['dety']*res['dety']))
               for res in results]

    # Plot data
    ax.errorbar(zenith, detxy, yerr=e_detxy, color=color, fmt='o', capsize=0,
                zorder=1, label=srcname)

    # Return
    return


# ======================================= #
# Plot zenith angle band for observations #
# ======================================= #
def plot_obs_zeniths(ax, obsname, color, srcname, y0=-0.015):
    """
    Plot observation results

    Parameters
    ----------
    ax : pyplot
        Subplot
    obsname : str
        Observation definition XML file
    color : str
        Color string
    srcname : str
        Source name
    y0 : float, optional
        Vertical axis location of lower edge of bar
    """
    # Get zenith angle range
    zenith_min = 90.0
    zenith_max =  0.0
    obs        = gammalib.GObservations(obsname)
    for run in obs:
        zenith = run.pointing().zenith()
        if zenith < zenith_min:
            zenith_min = zenith
        if zenith > zenith_max:
            zenith_max = zenith

    # Plot filled rectangle
    x0     = zenith_min
    width  = zenith_max - zenith_min
    height = 0.006
    ax.add_patch(Rectangle((x0,y0), width, height, alpha=1,
                           facecolor=color))
    ax.text(x0+0.5*width,y0+0.5*height, srcname, fontsize=7,
            horizontalalignment='center',
            verticalalignment='center')

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
        ROI radius (deg)

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
    Generate background model

    Parameters
    ----------
    obs : `~gammalib.GObservations`
        Observation container
    rad : float, optional
        ROI radius (deg)

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


# ============ #
# Plot results #
# ============ #
def plot_results(results, fontsize=12, show_source_obs=False):
    """
    Plot fit results

    Parameters
    ----------
    results : list of dict
        Fit results
    fontsize : int, optional
        Font size
    show_source_obs : bool, optional
        Show results for source observations
    """
    # Initialise figure
    fig = plt.figure(figsize=(6,4))
    fig.subplots_adjust(left=0.14, right=0.99, top=0.99, bottom=0.12,
                        wspace=0.3)

    # Get subplot
    ax = plt.subplot()

    # Generate arrays
    zenith  = [res['zenith'] for res in results]
    detx    = [abs(res['detx']) for res in results]
    e_detx  = [res['e_detx'] for res in results]
    dety    = [abs(res['dety']) for res in results]
    e_dety  = [res['e_dety'] for res in results]
    detxy   = [math.sqrt(res['detx']*res['detx']+res['dety']*res['dety']) for res in results]
    e_detxy = [math.sqrt(res['e_detx']*res['e_detx']*res['detx']*res['detx'] /
                         (res['detx']*res['detx']+res['dety']*res['dety']) +
                         res['e_dety']*res['e_dety']*res['dety']*res['dety'] /
                         (res['detx']*res['detx']+res['dety']*res['dety']))
               for res in results]
    #angle = [math.atan2(res['dety'],res['detx'])*gammalib.rad2deg for res in results]
    #print(angle)

    # Plot data
    ax.errorbar(zenith, detxy, yerr=e_detxy, fmt='ro', capsize=0, linewidth=1,
                zorder=3, label=r'$\sqrt{G_\mathrm{x}^2+G_\mathrm{y}^2}$')
    if not show_source_obs:
        ax.errorbar(zenith, detx, yerr=e_detx, color='gray', fmt='o', capsize=0, linewidth=1,
                    zorder=2, label=r'$|G_\mathrm{x}|$')
        ax.errorbar(zenith, dety, yerr=e_dety, color='silver', fmt='o', capsize=0, linewidth=1,
                    zorder=1, label=r'$|G_\mathrm{y}|$')

    # Plot source data
    if show_source_obs:
        plot_obs_results(ax,
                         '$HESSDATA/obs/obs_crab.xml',
                         '../crab/2019-04-26/crab_results_ptsrc_eplaw_lookup_grad_hess_edisp.xml',
                         'skyblue', 'Crab')
        plot_obs_results(ax,
                         '$HESSDATA/obs/obs_msh.xml',
                         '../msh15-52/2019-04-26/msh_results_egauss_eplaw_lookup_grad_hess.xml',
                         'darkseagreen', 'MSH 15-52')
        plot_obs_results(ax,
                         '$HESSDATA/obs/obs_rx.xml',
                         '../rxj1713/2019-04-26/rx_results_map0.40_eplaw_lookup_grad_hess_edisp.xml',
                         'plum', 'RX J1713-3946')
        plot_obs_results(ax,
                         '$HESSDATA/obs/obs_pks_flare.xml',
                         '../pks/2019-04-26/pks_t200_results_ptsrc_logparabola_lookup_grad_hess_edisp.xml',
                         'tan', 'PKS 2155-304')
    else:
        plot_obs_zeniths(ax, '$HESSDATA/obs/obs_crab.xml', 'skyblue', 'Crab', y0=-0.01)
        plot_obs_zeniths(ax, '$HESSDATA/obs/obs_msh.xml', 'darkseagreen', 'MSH 15-52', y0=-0.017)
        plot_obs_zeniths(ax, '$HESSDATA/obs/obs_rx.xml', 'plum', 'RX J1713-3946', y0=-0.024)
        plot_obs_zeniths(ax, '$HESSDATA/obs/obs_pks_flare.xml', 'tan', 'PKS 2155-304', y0=-0.031)

    # Plot blue line at zero gradient
    ax.plot([0.0,55.0], [0.0,0.0], color='b', linewidth=1, zorder=4)
    ax.set_xlim([0.0, 55.0])

    # Set labels
    ax.set_xlabel('Zenith angle (deg)', fontsize=fontsize)
    ax.set_ylabel('Gradient (deg$^{-1}$)', fontsize=fontsize)

    # Show legend
    if show_source_obs:
        ax.legend(loc='upper left', fontsize=8)
    else:
        ax.legend(loc='upper left')

    # Show figure
    plt.show()

    # Save figure
    fig.savefig('show_off_fits.png',dpi=300)

    # Return
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


# ================================ #
# Show off observations model fits #
# ================================ #
def show_off_fits(emin=0.2, emax=20.0):
    """
    Show off observations model fits

    Parameters
    ----------
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    """
    # Set filenames
    obsname = '$HESSDATA/obs/obs_off.xml'
    resname = 'show_off_fits.dmp'

    # If results do not yet exist then generate them
    if not os.path.isfile(resname):

        # Initialise results
        results = []

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

            # Set results
            result = {'id': run.id(),
                      'ontime': run.ontime(), 'livetime': run.livetime(),
                      'deadc': run.deadc(),
                      'ra': run.pointing().dir().ra_deg(),
                      'dec': run.pointing().dir().dec_deg(),
                      'zenith': run.pointing().zenith(),
                      'azimuth': run.pointing().azimuth(),
                      'detx': models[0].spatial()[1].value(),
                      'e_detx': models[0].spatial()[1].error(),
                      'dety': models[0].spatial()[2].value(),
                      'e_dety': models[0].spatial()[2].error()}

            # Append results
            results.append(result)

            # Stop after first run
            #break

        # Pickle results
        pickle.dump(results, open(resname, 'wb'))

    # Load results
    results = pickle.load(open(resname, 'rb'))

    # Plot analysis results
    plot_results(results)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show off observations model fits
    show_off_fits()
