#! /usr/bin/env python
# ==========================================================================
# Display butterfly generated by ctbutterfly and spectrum generated by csspec
#
# Copyright (C) 2018-2019 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation. either version 3 of the License, or
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
import cscripts
import matplotlib.pyplot as plt


# ============================================= #
# Return RX J1713 data (Aharonian et al. (2006) #
# ============================================= #
def rx_sed_2006():
    """
    Return SED data of Aharonian et al. (2006)
    """
    engs = [0.209,     0.253,     0.309,     0.366,     0.456,
            0.548,     0.649,     0.793,     0.951,     1.159,
            1.399,     1.705,     2.037,     2.480,     3.049,
            3.624,     4.460,     5.353,     6.460,     7.761,
            9.478,    11.552,    15.263,    22.765,    32.941]
    ergs = [2.163e-11, 2.417e-11, 3.065e-11, 3.352e-11, 2.908e-11,
            3.477e-11, 2.919e-11, 3.333e-11, 3.333e-11, 3.027e-11,
            2.432e-11, 2.713e-11, 2.622e-11, 2.888e-11, 2.279e-11,
            1.985e-11, 2.189e-11, 2.630e-11, 2.079e-11, 1.806e-11,
            1.616e-11, 1.041e-11, 5.372e-12, 4.526e-12, 6.388e-12]
    emin = [1.515e-11, 1.895e-11, 2.632e-11, 2.947e-11, 2.579e-11,
            3.189e-11, 2.606e-11, 3.067e-11, 3.046e-11, 2.752e-11,
            2.127e-11, 2.389e-11, 2.308e-11, 2.581e-11, 1.966e-11,
            1.671e-11, 1.869e-11, 2.308e-11, 1.781e-11, 1.517e-11,
            1.290e-11, 7.331e-12, 3.084e-12, 2.555e-12, 4.492e-12]
    yerr = [erg - emin[i] for i, erg in enumerate(ergs)]

    # Return data
    return engs, ergs, yerr


# ============================================= #
# Return RX J1713 data (Aharonian et al. (2011) #
# ============================================= #
def rx_sed_2011():
    """
    Return SED data of Aharonian et al. (2011)
    """
    # Set data
    engs   = [    0.33,     0.40,     0.49,     0.59,     0.71,     0.86,
                  1.04,     1.26,     1.53,     1.85,     2.24,     2.71,
                  3.28,     3.98,     4.81,     5.82,     7.05,     8.53,
                 10.33,    12.51,    15.41,    18.32,    22.18,    26.85,
                 32.50,    47.19,    81.26]
    ebins  = [    0.07,     0.07,     0.10,     0.11,     0.14,     0.16,
                  0.20,     0.24,     0.30,     0.35,     0.43,     0.52,
                  0.64,     0.76,     0.92,     1.12,     1.36,     1.64,
                  1.98,     2.40,     2.91,     3.52,     4.26,     5.16,
                  6.25,    27.80,    49.31]
    flux   = [2.29e-10, 1.25e-10, 9.46e-11, 6.06e-11, 4.37e-11, 2.15e-11,
              1.82e-11, 1.17e-11, 8.87e-12, 5.63e-12, 3.78e-12, 2.49e-12,
              1.64e-12, 1.04e-12, 7.48e-13, 4.34e-13, 2.32e-13, 1.25e-13,
              1.07e-13, 5.61e-14, 2.17e-14, 1.84e-14, 9.24e-15, 7.40e-15,
              6.46e-15, 9.63e-16, 1.98e-16]
    e_flux = [0.32e-10, 0.16e-10, 0.90e-11, 0.52e-11, 0.31e-11, 0.18e-11,
              0.11e-11, 0.07e-11, 0.50e-12, 0.33e-12, 0.23e-12, 0.16e-12,
              0.11e-12, 0.08e-12, 0.57e-13, 0.40e-13, 0.30e-13, 0.20e-13,
              0.14e-13, 0.97e-14, 0.69e-14, 0.46e-14, 3.16e-15, 2.06e-15,
              1.50e-15, 3.93e-16, 1.29e-16]

    # Convert per TeV in per ergs
    for i, _ in enumerate(flux):
        flux[i]   = flux[i]   * engs[i] * engs[i] * 1.60218
        e_flux[i] = e_flux[i] * engs[i] * engs[i] * 1.60218

    # Return data
    return engs, flux, e_flux


# =========================================== #
# Return RX J1713 data (Abdalla et al. (2018) #
# =========================================== #
def rx_sed_2018():
    """
    Return SED data of Abdalla et al. (2018)
    """
    # Set data
    engs   = [    0.23,     0.27,     0.32,     0.38,     0.45,     0.53,
                  0.63,     0.75,     0.89,     1.05,     1.25,     1.48,
                  1.76,     2.08,     2.47,     2.93,     3.47,     4.12,
                  4.88,     5.79,     6.86,     8.14,     9.65,    11.44,
                 13.56,    16.08,    19.99,    24.62,    29.19,    34.61]
    flux   = [4.44e-11, 3.99e-11, 4.28e-11, 3.92e-11, 3.96e-11, 3.38e-11,
              3.43e-11, 3.41e-11, 3.21e-11, 2.94e-11, 3.02e-11, 3.25e-11,
              2.99e-11, 3.09e-11, 2.62e-11, 2.54e-11, 2.59e-11, 2.39e-11,
              2.50e-11, 2.34e-11, 2.01e-11, 1.74e-11, 1.42e-11, 1.58e-11,
              8.34e-12, 9.67e-12, 3.82e-12, 5.68e-12, 2.79e-12, 5.08e-12]
    e_flux = [8.35e-12, 5.17e-12, 3.52e-12, 2.41e-12, 2.04e-12, 1.73e-12,
              1.41e-12, 1.34e-12, 1.22e-12, 1.12e-12, 1.09e-12, 1.07e-12,
              1.01e-12, 1.19e-12, 1.24e-12, 1.17e-12, 1.18e-12, 1.28e-12,
              1.23e-12, 1.31e-12, 1.38e-12, 1.51e-12, 1.53e-12, 1.59e-12,
              1.70e-12, 1.67e-12, 1.41e-12, 1.75e-12, 1.57e-12, 1.27e-12]

    # Convert per TeV in per ergs
    for i, _ in enumerate(flux):
        flux[i]   = flux[i]
        e_flux[i] = e_flux[i]

    # Return data
    return engs, flux, e_flux


# =========================== #
# Plot spectrum and butterfly #
# =========================== #
def plot_spectrum(specfile, butterfile, ax, fmt='ro'):
    """
    Plot sensitivity data

    Parameters
    ----------
    specfile : str
        Name of spectrum FITS file
    butterfile : str
        Name of butterfly FITS file
    ax : pyplot
        Subplot
    fmt : str
        Format string
    """
    # Read spectrum file    
    fits     = gammalib.GFits(specfile)
    table    = fits.table(1)
    c_energy = table['Energy']
    c_ed     = table['ed_Energy']
    c_eu     = table['eu_Energy']
    c_flux   = table['Flux']
    c_eflux  = table['e_Flux']
    c_ts     = table['TS']
    c_upper  = table['UpperLimit']

    # Read butterfly file
    csv = gammalib.GCsv(butterfile)

    # Initialise arrays to be filled
    energies    = []
    flux        = []
    ed_engs     = []
    eu_engs     = []
    e_flux      = []
    ul_energies = []
    ul_ed_engs  = []
    ul_eu_engs  = []
    ul_flux     = []
    butterfly_x = []
    butterfly_y = []
    line_x      = []
    line_y      = []

    # Loop over rows of the file
    nrows = table.nrows()
    for row in range(nrows):

        # Get Test Statistic, flux and flux error
        ts    = c_ts.real(row)
        flx   = c_flux.real(row)
        e_flx = c_eflux.real(row)

        # If Test Statistic is larger than 9 and flux error is smaller than
        # flux then append flux plots ...
        if ts > -99999.0 and e_flx < flx and e_flx > 0.0:
            energies.append(c_energy.real(row))
            flux.append(c_flux.real(row))
            ed_engs.append(c_ed.real(row))
            eu_engs.append(c_eu.real(row))
            e_flux.append(c_eflux.real(row))

        # ... otherwise append upper limit
        else:
            ul_energies.append(c_energy.real(row))
            ul_flux.append(c_upper.real(row))
            ul_ed_engs.append(c_ed.real(row))
            ul_eu_engs.append(c_eu.real(row))

    # Set upper limit errors
    yerr = [0.6 * x for x in ul_flux]

    # Loop over rows of the file
    nrows = csv.nrows()
    for row in range(nrows):

        # Limit values to positive
        value = csv.real(row,2)
        if value < 1e-40:
            value = 1e-40

        # Compute upper edge of confidence band
        butterfly_x.append(csv.real(row,0)/1.0e6) # TeV
        butterfly_y.append(value*csv.real(row,0)*csv.real(row,0)*gammalib.MeV2erg)

        # Set line values
        line_x.append(csv.real(row,0)/1.0e6) # TeV
        line_y.append(csv.real(row,1)*csv.real(row,0)*csv.real(row,0)*gammalib.MeV2erg)

    # Loop over the rows backwards to compute the lower edge of the
    # confidence band
    for row in range(nrows):
        index = nrows - 1 - row
        butterfly_x.append(csv.real(index,0)/1.0e6)
        low_error = max(csv.real(index,3)*csv.real(index,0)*csv.real(index,0)*gammalib.MeV2erg, 1e-26)
        butterfly_y.append(low_error)

    # Plot the spectrum
    ax.errorbar(energies, flux, yerr=e_flux, xerr=[ed_engs, eu_engs],
                fmt=fmt, label='ctools')
    ax.errorbar(ul_energies, ul_flux, xerr=[ul_ed_engs, ul_eu_engs],
                yerr=yerr, uplims=True, fmt=fmt, label='_nolegend_')

    # Plot RX J1713 SED
    #engs, ergs, yerrs = rx_sed_2006()
    #engs, ergs, yerrs = rx_sed_2011()
    engs, ergs, yerrs = rx_sed_2018()
    ax.errorbar(engs, ergs, yerr=yerrs, fmt='bo', label='H.E.S.S. collaboration et al. (2018)')

	# Plot spectral law from publication
    abdalla_x = [i*0.1+0.2 for i in range(400)]
    abdalla_y = []
    for x in abdalla_x:
        y = 2.3e-11 * math.pow(x, -2.06) * math.exp(-x/12.9)
        abdalla_y.append(y*x*x*gammalib.MeV2erg*1.0e6)
    ax.plot(abdalla_x, abdalla_y, color='b', ls='-', label='_nolegend_')

    # Plot the butterfly
    ax.plot(line_x, line_y, color='red', ls='-')
    ax.fill(butterfly_x, butterfly_y, color='salmon', alpha=0.5)

    # Add legend
    ax.legend()

    # Return
    return


# ============= #
# Show spectrum #
# ============= #
def show_spectrum():
    """
    Show spectrum
    """
    # Create figure
    fig = plt.figure('Spectra', (7, 4))
    fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.15)
    ax = plt.subplot()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (TeV)', fontsize=14)
    ax.set_ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)', fontsize=14)
    ax.set_xlim([0.3,40.0])
    ax.set_ylim([1e-12,1e-10])
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.grid()

    # Set usage string
    usage = 'show_rx_spectrum.py [-p plotfile] [specfile1] [butterfly1]'

    # Set default options
    options = [{'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Show spectrum
    plot_spectrum(args[0], args[1], ax, fmt='ro')

    # Show and save figure
    plt.show()
    fig.savefig('rx_sed.png',dpi=300)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show spectrum
    show_spectrum()
