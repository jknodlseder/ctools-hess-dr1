#! /usr/bin/env python
# ==========================================================================
# Display PKS 2155-304 spectrum and butterfly
#
# Copyright (C) 2019 Juergen Knoedlseder
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
import gammalib
import cscripts
import matplotlib.pyplot as plt


# ======================== #
# Return PKS 2155-304 data #
# ======================== #
def pks_sed():
    """
    Return SED data of Aharonian et al. (2009)
    """
    # Set data points
    engs = [0.22013, 0.26655, 0.32227, 0.38989, 0.47216,
            0.57465, 0.69447, 0.83883, 1.01583, 1.22748,
            1.49510, 1.81097, 2.18739, 2.64572, 3.22640,
            3.90068, 4.71999, 5.70876]
    flux = [7.19578e-10, 5.94289e-10, 5.40850e-10, 4.55817e-10, 3.52754e-10,
            2.71342e-10, 2.01465e-10, 1.49872e-10, 1.09441e-10, 8.36712e-11,
            4.80279e-11, 4.86552e-11, 3.42818e-11, 2.13132e-11, 1.35004e-11,
            1.16156e-11, 3.65306e-12, 3.68269e-12]
    flux_max = [7.19578e-10, 5.94289e-10, 5.40850e-10, 4.55817e-10, 3.52754e-10,
                2.71342e-10, 2.01465e-10, 1.49872e-10, 1.09441e-10, 8.36712e-11,
                4.80279e-11, 4.86552e-11, 4.01484e-11, 2.59826e-11, 1.79978e-11,
                1.55800e-11, 6.99166e-12, 7.27051e-12]
    flux_min = [7.19578e-10, 5.94289e-10, 5.40850e-10, 4.55817e-10, 3.52754e-10,
                2.71342e-10, 2.01465e-10, 1.49872e-10, 1.09441e-10, 8.36712e-11,
                4.80279e-11, 4.86552e-11, 2.90317e-11, 1.68057e-11, 9.65650e-12,
                7.71089e-12, 1.03840e-12, 1.16634e-12]

    # Compute flux errors
    for i, _ in enumerate(flux):
         flux_max[i] = flux_max[i] - flux[i]
         flux_min[i] = flux[i] - flux_min[i]

    # Return data
    return engs, flux, [flux_max, flux_min]


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
    ax : plt
        Frame
    fmt : str, optional
        Format for ctools points
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

    # Plot PKS 2155-304 SED
    engs, pks, e_pks = pks_sed()
    ax.errorbar(engs, pks, yerr=e_pks, fmt='bo', label='Aharonian et al. (2009)')
    ax.errorbar([6.90492], [9.76439e-12], yerr=0.6*9.76439e-12, uplims=True,
                fmt='bo', label='_nolegend_')

    # Plot the spectrum
    ax.errorbar(energies, flux, yerr=e_flux, xerr=[ed_engs, eu_engs],
                 fmt=fmt, label='ctools')
    ax.errorbar(ul_energies, ul_flux, xerr=[ul_ed_engs, ul_eu_engs],
                yerr=yerr, uplims=True, fmt=fmt, label='_nolegend_')

    # Plot the butterfly
    ax.plot(line_x, line_y, color='red',ls='-')
    ax.fill(butterfly_x, butterfly_y, color='salmon', alpha=0.5)

    # Resort labels
    handles, labels = ax.get_legend_handles_labels()
    order           = [1,0]

    # Add legend
    ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper right')

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
    fig = plt.figure('Spectra', (7,4))
    fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.15)
    ax = plt.subplot()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (TeV)', fontsize=14)
    ax.set_ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)', fontsize=14)
    ax.set_xlim([0.15,15.0])
    ax.set_ylim([8e-13,4e-9])
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.grid()

    # Set usage string
    usage = 'show_spectrum.py [-p plotfile] [specfile] [butterfly]'

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
    fig.savefig('pks_sed.png',dpi=300)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show spectrum
    show_spectrum()
