#! /usr/bin/env python
# ==========================================================================
# Show background lookup table
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
import cscripts
import matplotlib.pyplot as plt


# ================= #
# Plot lookup table #
# ================= #
def plot_lookup(filename, ax, nengs=300, nthetas=300):
    """
    Plot lookup table

    Parameters
    ----------
    filename : str
        Lookup table filename
    ax : figure
        Subplot
    nengs : int, optional
        Number of energies
    nthetas : int, optional
        Number of offset angles
    """
    # Load lookup table
    lookup = gammalib.GCTAModelSpatialLookup(filename)
    table  = lookup.table()

    # Get boundaries
    axis_eng   = table.axis('ENERG')
    axis_theta = table.axis('THETA')
    emin       = table.axis_lo(axis_eng, 0)
    emax       = table.axis_hi(axis_eng, table.axis_bins(axis_eng)-1)
    tmin       = table.axis_lo(axis_theta, 0)
    tmax       = table.axis_hi(axis_theta, table.axis_bins(axis_theta)-1)

    # Use log energies
    emin = math.log10(emin)
    emax = math.log10(emax)

    # Set axes
    denergy     = (emax - emin)/(nengs-1)
    dtheta      = (tmax - tmin)/(nthetas-1)
    logenergies = [emin+i*denergy for i in range(nengs)]
    thetas      = [tmax-i*dtheta  for i in range(nthetas)]

    # Initialise image
    image = []

    # Set event time
    time = gammalib.GTime()

    # Loop over offset angles
    for theta in thetas:

        # Set event direction
        dir = gammalib.GCTAInstDir(theta*gammalib.deg2rad,0.0)

        # Initialise row
        row = []

        # Loop over energies
        for logenergy in logenergies:

            # Set event energy
            energy = gammalib.GEnergy()
            energy.log10TeV(logenergy)

            # Get lookup table value
            value = lookup.eval(dir, energy, time)

            # Append value
            row.append(value)

        # Append row
        image.append(row)

    # Plot image
    c    = ax.imshow(image, extent=[emin,emax,tmin,tmax],
                     cmap=plt.get_cmap('jet'), aspect=0.8)
    cbar = plt.colorbar(c, orientation='vertical', shrink=0.7)
    cbar.set_label('sr$^{-1}$', labelpad=-22, y=-0.05, rotation=0)

    # Plot title and axis
    ax.set_xlabel('log10(E/TeV)', fontsize=12)
    ax.set_ylabel('Offset angle (deg)', fontsize=12)

    # Return
    return


# ========================== #
# Generate background lookup #
# ========================== #
def generate_background_lookup(filename):
    """
    Generate background lookup from empty field observations

    Parameters
    ----------
    filename : str
        Lookup table filename
    """
    # Set filenames
    obsname  = '$HESSDATA/obs/obs_off.xml'

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


# ============================ #
# Show background lookup table #
# ============================ #
def show_bkg_lookup():
    """
    Show background lookup table
    """
    # Set usage string
    usage = 'show_bkg_lookup.py file'

    # Set default options
    options = []

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Set filename
    if len(args) > 0:
        filename = args[0]
    else:
        filename = 'bkg_lookup.fits'
        generate_background_lookup(filename)

    # Create figure
    fig = plt.figure('Background lookup table', (6,3.6))
    fig.subplots_adjust(left=0.12, right=1.05, top=1.00, bottom=0.10)
    ax = plt.subplot()

    # Plot background lookup table
    plot_lookup(filename, ax)

    # Show map
    plt.show()

    # Save figure
    pltname = os.path.splitext(filename)[0] + '.png'
    fig.savefig(pltname,dpi=300)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show background lookup table
    show_bkg_lookup()
