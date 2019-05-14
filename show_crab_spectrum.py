#! /usr/bin/env python
# ==========================================================================
# Show Crab spectrum
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


# =========================== #
# Plot spectrum and butterfly #
# =========================== #
def plot_spectrum(specfile, butterfile, ax, fmt='ro', color='green'):
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
    color : str
        Color string
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
        low_error = max(csv.real(index,3)*csv.real(index,0)*csv.real(index,0)*gammalib.MeV2erg, 1e-20)
        butterfly_y.append(low_error)

    # Plot the spectrum
    ax.errorbar(energies, flux, yerr=e_flux, xerr=[ed_engs, eu_engs],
                fmt=fmt, label='ctools')
    ax.errorbar(ul_energies, ul_flux, xerr=[ul_ed_engs, ul_eu_engs],
                yerr=yerr, uplims=True, fmt=fmt, label='_nolegend_')

    # Plot Crab SED from Aharonian et al. (2006a)
    hess_x = [i*0.1+0.4 for i in range(400)]
    hess_y = []
    for x in hess_x:
        y = 3.84e-11 * math.pow(x, -2.41) * math.exp(-x/15.1)
        hess_y.append(y*x*x*gammalib.MeV2erg*1.0e6)
    ax.plot(hess_x, hess_y, color='blue', ls='-', label='Aharonian et al. (2006a)')

    # Plot Crab SED from Nigro et al. (2019)
    nigro_x = [i*0.1+0.4 for i in range(400)]
    nigro_y = []
    for x in nigro_x:
        y = 4.47e-11 * math.pow(x, -2.39 - 0.16 * math.log(x))
        nigro_y.append(y*x*x*gammalib.MeV2erg*1.0e6)
    ax.plot(nigro_x, nigro_y, color='green', ls='-', label='Nigro et al. (2019)')

    # Plot the butterfly
    ax.plot(line_x, line_y, color='red', ls='-')
    ax.fill(butterfly_x, butterfly_y, color='salmon', alpha=0.5)

    # Resort labels
    handles, labels = ax.get_legend_handles_labels()
    order           = [2,0,1]

    # Add legend
    ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

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
    ax.set_xlim([0.4,40.0])
    ax.set_ylim([1e-12,1e-10])
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.grid()

    # Set usage string
    usage = 'show_crab_spectrum.py [specfile1] [butterfly1]'

    # Set default options
    options = [{'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Show spectrum
    plot_spectrum(args[0], args[1], ax, fmt='ro')

    # Save figure
    plt.show()
    fig.savefig('crab_sed.png',dpi=300)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show spectrum
    show_spectrum()
