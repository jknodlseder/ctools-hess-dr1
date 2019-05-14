#! /usr/bin/env python
# ==========================================================================
# Display MSH 15-52 spectrum and butterfly
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


# ===================== #
# Return MSH 15-52 data #
# ===================== #
def msh_sed():
    """
    Return SED data of Aharonian et al. (2005)
    """
    ebins = [[0.28,0.34], [0.34,0.49], [0.49,0.73], [0.73,1.1],
             [1.1,1.6],   [1.6,2.3],   [2.3,3.6],   [3.6,4.8],
             [4.8,7.6],   [7.6,10.4],  [10.4,15.7], [15.7,23.2],
             [23.2,34.8], [34.8,46.2]]
    engs = [0.31, 0.41, 0.60, 0.89,
            1.31, 1.90, 2.87, 4.15,
            6.09, 8.90, 12.8, 19.1,
            28.4, 40.1]
    flux = [1.3e-10, 3.6e-11, 1.9e-11, 7.2e-12,
            3.2e-12, 1.6e-12, 5.5e-13, 2.2e-13,
            9.5e-14, 4.7e-14, 1.2e-14, 6.2e-15,
            2.4e-15, 7.8e-16]
    e_flux = [3.2e-11, 3.8e-12, 1.6e-12, 6.5e-13,
              6.5e-13, 1.5e-13, 7.0e-14, 3.4e-14,
              1.7e-14, 9.5e-15, 4.1e-15, 2.0e-15,
              1.2e-15, 5.5e-16]

    # Convert per TeV in per ergs
    for i, _ in enumerate(flux):
        flux[i]   = flux[i]   * engs[i] * engs[i] * 1.60218
        e_flux[i] = e_flux[i] * engs[i] * engs[i] * 1.60218

    # Return data
    return engs, ebins, flux, e_flux


# =========================== #
# Plot spectrum and butterfly #
# =========================== #
def plot_spectrum(specfile, butterfile1, butterfile2, ax, fmt='ro'):
    """
    Plot sensitivity data

    Parameters
    ----------
    specfile : str
        Name of spectrum FITS file
    butterfile1 : str
        Name of first butterfly text file
    butterfile2 : str
        Name of second butterfly text file
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

    # Read butterfly files
    csv1 = gammalib.GCsv(butterfile1)
    csv2 = gammalib.GCsv(butterfile2)

    # Initialise arrays to be filled
    energies     = []
    flux         = []
    ed_engs      = []
    eu_engs      = []
    e_flux       = []
    ul_energies   = []
    ul_ed_engs   = []
    ul_eu_engs   = []
    ul_flux      = []
    butterfly1_x = []
    butterfly1_y = []
    butterfly2_x = []
    butterfly2_y = []
    line1_x      = []
    line1_y      = []
    line2_x      = []
    line2_y      = []

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

    # Loop over rows of the first butterfly file
    nrows = csv1.nrows()
    for row in range(nrows):

        # Limit values to positive
        value = csv1.real(row,2)
        if value < 1e-40:
            value = 1e-40

        # Compute upper edge of confidence band
        butterfly1_x.append(csv1.real(row,0)/1.0e6) # TeV
        butterfly1_y.append(value*csv1.real(row,0)*csv1.real(row,0)*gammalib.MeV2erg)

        # Set line values
        line1_x.append(csv1.real(row,0)/1.0e6) # TeV
        line1_y.append(csv1.real(row,1)*csv1.real(row,0)*csv1.real(row,0)*gammalib.MeV2erg)

    # Loop over the rows backwards to compute the lower edge of the
    # confidence band
    for row in range(nrows):
        index = nrows - 1 - row
        butterfly1_x.append(csv1.real(index,0)/1.0e6)
        low_error = max(csv1.real(index,3)*csv1.real(index,0)*csv1.real(index,0)*gammalib.MeV2erg, 1e-26)
        butterfly1_y.append(low_error)

    # Loop over rows of the second butterfly file
    nrows = csv2.nrows()
    for row in range(nrows):

        # Limit values to positive
        value = csv2.real(row,2)
        if value < 1e-40:
            value = 1e-40

        # Compute upper edge of confidence band
        butterfly2_x.append(csv2.real(row,0)/1.0e6) # TeV
        butterfly2_y.append(value*csv2.real(row,0)*csv2.real(row,0)*gammalib.MeV2erg)

        # Set line values
        line2_x.append(csv2.real(row,0)/1.0e6) # TeV
        line2_y.append(csv2.real(row,1)*csv2.real(row,0)*csv2.real(row,0)*gammalib.MeV2erg)

    # Loop over the rows backwards to compute the lower edge of the
    # confidence band
    for row in range(nrows):
        index = nrows - 1 - row
        butterfly2_x.append(csv2.real(index,0)/1.0e6)
        low_error = max(csv2.real(index,3)*csv2.real(index,0)*csv2.real(index,0)*gammalib.MeV2erg, 1e-26)
        butterfly2_y.append(low_error)

    # Plot the spectrum 
    ax.errorbar(energies, flux, yerr=e_flux, xerr=[ed_engs, eu_engs],
                 fmt=fmt, label='ctools')
    ax.errorbar(ul_energies, ul_flux, xerr=[ul_ed_engs, ul_eu_engs],
                yerr=yerr, uplims=True, fmt=fmt, label='_nolegend_')

    # Plot MSH 15-52 SED
    engs, ebins, flux, e_flux = msh_sed()
    ax.errorbar(engs, flux, yerr=e_flux, fmt='bo', label='Aharonian et al. (2005)')

    # Plot MSH 15-52 spectrum from Dubois (2009)
    dubois_x = [i*0.1+0.2 for i in range(1000)]
    dubois_y = []
    for x in dubois_x:
        y = 7.13e-12 * math.pow(x, -2.06) * math.exp(-x/11.9)
        dubois_y.append(y*x*x*gammalib.MeV2erg*1.0e6)
    ax.plot(dubois_x, dubois_y, color='blue', ls='-', label='Dubois (2009)')

    # Plot MSH 15-52 spectrum from HGPS
    #hgps_x = [i*0.1+0.2 for i in range(1000)]
    #hgps_y = []
    #for x in hgps_x:
    #    y = 6.859766e-12 * math.pow(x, -2.0507) * math.exp(-x/19.2)
    #    hgps_y.append(y*x*x*gammalib.MeV2erg*1.0e6)
    #ax.plot(hgps_x, hgps_y, color='blue', ls='-', label='HGPS')

    # Plot the first butterfly
    ax.plot(line1_x, line1_y, color='red',ls='--')
    ax.fill(butterfly1_x, butterfly1_y, color='lightsalmon', alpha=0.5)

    # Plot the second butterfly
    ax.plot(line2_x, line2_y, color='red',ls='-')
    ax.fill(butterfly2_x, butterfly2_y, color='salmon', alpha=0.5)

    # Resort labels
    handles, labels = ax.get_legend_handles_labels()
    order           = [1,2,0]

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
    ax.set_xlim([0.2,50.0])
    ax.set_ylim([1e-13,5e-11])
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.grid()

    # Set usage string
    usage = 'show_msh_spectrum.py [-p plotfile] [specfile1] [butterfly1] [butterfly2]'

    # Set default options
    options = [{'option': '-p', 'value': ''}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    plotfile = options[0]['value']

    # Show spectrum
    plot_spectrum(args[0], args[1], args[2], ax, fmt='ro')

    # Save figure
    plt.show()
    fig.savefig('msh_sed.png',dpi=300)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show spectrum
    show_spectrum()
