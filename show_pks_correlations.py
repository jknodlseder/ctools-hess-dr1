#! /usr/bin/env python
# ==========================================================================
# Display correlation between ctools and H.E.S.S. lightcurve data
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
import gammalib
import cscripts
import matplotlib.pyplot as plt


# =================== #
# Get HESS lightcurve #
# =================== #
def get_hess_lightcurve(filename='pks-hess-data.csv'):
    """
    Get HESS lightcurve

    Parameters
    ----------
    filename : str, optional
        Light curve data file name

    Returns
    -------
    x, y, y_err : tuple of floats
        Energy, flux and flux error
    """
    # Open data file
    file1 = gammalib.GFilename('$CTAGITROOT/analysis/hess_dr1/%s' % filename)
    file2 = gammalib.GFilename('%s' % filename)
    if file1.exists():
        csv = gammalib.GCsv(file1,'\t')
    else:
        csv = gammalib.GCsv(file2,'\t')

    # Initialise arrays
    x     = []
    y     = []
    y_err = []

    # Loop over data
    for row in range(csv.nrows()):
        x.append(float(csv[row,0])+1.0/(24.0*60.0))
        y.append(float(csv[row,1])*1.0e-9)
        y_err.append((float(csv[row,2])-float(csv[row,1]))*1.0e-9)

    # Return
    return x,y,y_err


# =================== #
# Get HESS lightcurve #
# =================== #
def get_hess_lightcurve2(filename='pks-hess-data2.csv'):
    """
    Get HESS lightcurve

    Parameters
    ----------
    filename : str, optional
        Light curve data file name

    Returns
    -------
    x, y, y_err : tuple of floats
        Energy, flux and flux error
    """
    # Open data file
    file1 = gammalib.GFilename('$CTAGITROOT/analysis/hess_dr1/%s' % filename)
    file2 = gammalib.GFilename('%s' % filename)
    if file1.exists():
        csv = gammalib.GCsv(file1,' ')
    else:
        csv = gammalib.GCsv(file2,' ')

    # Initialise arrays
    x         = []
    y         = []
    y_err_min = []
    y_err_max = []

    # Loop over data
    for row in range(csv.nrows()):
        x.append(float(csv[row,0])-53900.0+1.0/(24.0*60.0))
        y.append(float(csv[row,2]))
        y_err_min.append(float(csv[row,3]))
        y_err_max.append(float(csv[row,4]))

    # Return
    return x,y,[y_err_min,y_err_max]


# =============== #
# Plot lightcurve #
# =============== #
def plot_lightcurve(filename, label='ctools (On-Off - wstat)'):
    """
    Plot lightcurve

    Parameters
    ----------
    filename : str
        Name of lightcurve FITS file
    label : str, optional
        Label string
    """
    # Read spectrum file    
    fits    = gammalib.GFits(filename)
    table   = fits.table(1)
    
    # Extract standard columns
    c_mjd   = table['MJD']
    c_emjd  = table['e_MJD']
    c_ts    = table['TS']

    # Extract columns dependent on flux type
    if table.contains('EnergyFlux'):
        c_flux  = table['EnergyFlux']
        c_eflux = table['e_EnergyFlux']
        c_upper = table['EFluxUpperLimit']
        ylabel  = r'E $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)'
    elif table.contains('PhotonFlux'):
        c_flux  = table['PhotonFlux']
        c_eflux = table['e_PhotonFlux']
        c_upper = table['FluxUpperLimit']
        ylabel  = r'flux(>700 GeV) (ph cm$^{-2}$ s$^{-1}$)'
    else:
        c_flux  = table['Prefactor']
        c_eflux = table['e_Prefactor']
        c_upper = table['DiffUpperLimit']
        ylabel  = r'dN/dE (cm$^{-2}$ s$^{-1}$ MeV$^{-1}$)'

    # Initialise arrays to be filled
    flux             = []
    e_flux           = []
    ul_flux          = []
    ul_e_flux        = []
    hess_flux        = []
    hess_flux_min    = []
    hess_flux_max    = []
    ul_hess_flux     = []
    ul_hess_flux_min = []
    ul_hess_flux_max = []

    # Get H.E.S.S. lightcurve
    #hess_x, hess_y, hess_y_err = get_hess_lightcurve()
    hess_x, hess_y, hess_y_err = get_hess_lightcurve2()

    # Loop over rows of the file
    nrows = table.nrows()
    for row in range(nrows):

        # Get MJD, Test Statistic, flux and flux error
        mjd   = c_mjd.real(row)-53900.0
        ts    = c_ts.real(row)
        flx   = c_flux.real(row)
        e_flx = c_eflux.real(row)

        # Find corresponding H.E.S.S. flux
        max_delta_mjd = 1.0e30
        for i, hess_mjd in enumerate(hess_x):
            delta_mjd = abs(hess_mjd-mjd)
            if delta_mjd < max_delta_mjd:
                max_delta_mjd      = delta_mjd
                best_hess_mjd      = hess_mjd
                best_hess_flux     = hess_y[i]
                best_hess_flux_min = hess_y_err[0][i]
                best_hess_flux_max = hess_y_err[0][i]
        if (max_delta_mjd*24.0*60.0) > 0.5:
            print('No H.E.S.S. flux found. Skip data point.')
            continue

        # If Test Statistic is larger than 9 and twice the flux error is
        # smaller than the flux, then append flux point ...
        if ts > 0.0 and e_flx > 0.0:
            flux.append(c_flux.real(row))
            e_flux.append(c_eflux.real(row))
            hess_flux.append(best_hess_flux)
            hess_flux_min.append(best_hess_flux_min)
            hess_flux_max.append(best_hess_flux_max)

        # ... otherwise append upper limit
        else:
            #if c_upper.real(row) > 0.6e-10:
            #    continue
            ul_flux.append(c_upper.real(row))
            ul_e_flux.append(0.5*c_upper.real(row))
            ul_hess_flux.append(best_hess_flux)
            ul_hess_flux_min.append(best_hess_flux_min)
            ul_hess_flux_max.append(best_hess_flux_max)

    # Create figure
    fig = plt.figure('Flux correlation', (8,6))
    ax  = plt.subplot()
    ax.set_xlabel('H.E.S.S. '+ylabel, fontsize=14)
    ax.set_ylabel('ctools '+ylabel, fontsize=14)
    ax.set_xlim([-0.2e-10, 3.4e-10])
    ax.set_ylim([-0.2e-10, 3.4e-10])
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.grid()

    # Plot the light curve 
    ax.plot([-0.2e-10,3.4e-10], [-0.2e-10,3.4e-10], color='b', linewidth=2, zorder=3)
    ax.errorbar(hess_flux, flux, xerr=[hess_flux_min, hess_flux_max],
                yerr=e_flux, fmt='ro', zorder=2)
    ax.errorbar(ul_hess_flux, ul_flux, xerr=[ul_hess_flux_min, ul_hess_flux_max], yerr=ul_e_flux,
                 uplims=True, fmt='ro', zorder=1)

    # Show title
    fig.suptitle(label, fontsize=16)

    # Show figure
    plt.show()
    if '_3D' in filename:
        fig.savefig('pks_correlation_unbinned.png',dpi=300)
    else:
        fig.savefig('pks_correlation_onoff.png',dpi=300)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Plot 3D lightcurve
    plot_lightcurve('pks_lightcrv_ptsrc_fplaw700_lookup_grad_hess_edisp_3D.fits',
                    label='Unbinned')

    # Plot OnOff lightcurve
    plot_lightcurve('pks_lightcrv_ptsrc_fplaw700_lookup_grad_hess_edisp_ONOFF_10bins.fits',
                     label='On-Off - wstat')
