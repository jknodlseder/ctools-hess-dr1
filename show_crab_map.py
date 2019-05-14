#! /usr/bin/env python
# ==========================================================================
# Show Crab sky map
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
import glob
import math
import gammalib
import matplotlib.pyplot as plt


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
        Map smoothing parameter (deg)
    """
    # Optionally smooth map
    if smooth > 0.0:
        map = inmap.copy()
        map.smooth('GAUSSIAN', smooth)
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
    aspect  = 1.0

    # Show Aitoff projection
    c    = sub.imshow(array, extent=(ra_min,ra_max,dec_min,dec_max),
                      cmap=plt.get_cmap('jet'), aspect=aspect)
    cbar = plt.colorbar(c, orientation='vertical', shrink=0.9)
    sub.set_xlabel('Right Ascension (deg)', fontsize=14)
    sub.set_ylabel('Declination (deg)', fontsize=14)

    # Set range
    sub.set_xlim([ra_min, ra_max])
    sub.set_ylim([dec_min, dec_max])

    # Return
    return


# =========== #
# Plot circle #
# =========== #
def plot_circle(ra, dec, rad, color='w', fmt='-', linewidth=1):
    """
    Plot circle

    Parameters
    ----------
    ra : float
        Right Ascension (deg)
    dec : float
        Declination (deg)
    rad : float
        Radius (deg)
    color : str, optional
        Color string
    fmt : str, optional
        Format string
    linewidth : int, optional
        Line width
    """
    # Create arrays
    x = [ra +rad*math.sin(angle*gammalib.deg2rad) for angle in range(360)]
    y = [dec+rad*math.cos(angle*gammalib.deg2rad) for angle in range(360)]

    # Plot circle
    plt.plot(x, y, fmt, color=color, linewidth=linewidth)

    # Return
    return


# ============ #
# Plot circles #
# ============ #
def plot_circles():
    """
    Plot Off region circles
    """
    # Set colors
    colors = ['cyan', 'lime', 'yellow', 'pink']
    
    # Get on and off circles
    ons  = glob.glob('onoff-joint40_cstat_ptsrc_on.reg')
    offs = glob.glob('onoff-joint40_cstat_ptsrc_*_off.reg')

    # Loop over on files
    for file in ons:

        # Load regions
        regions = gammalib.GSkyRegions(file)

        # Loop over On regions
        for region in regions:

            # Get RA, Dec and radius
            ra     = region.ra()
            dec    = region.dec()
            radius = region.radius()

            # Plot circle
            plot_circle(ra, dec, radius, linewidth=2)

    # Loop over off files
    for i, file in enumerate(offs):

        # Load regions
        regions = gammalib.GSkyRegions(file)

        # Loop over Off regions
        for region in regions:

            # Get RA, Dec and radius
            ra     = region.ra()
            dec    = region.dec()
            radius = region.radius()

            # Plot circle
            plot_circle(ra, dec, radius, color=colors[i], fmt='--', linewidth=2)

    # Return
    return


# ============== #
# Plot pointings #
# ============== #
def plot_pointings():
    """
    Plot pointings
    """
    # Set colors
    colors = ['cyan', 'lime', 'yellow', 'pink']

    # Load observations
    obs = gammalib.GObservations('obs_crab_selected.xml')

    # Loop over observations
    for i, run in enumerate(obs):

        # Get pointing
        pnt = run.pointing().dir()
        ra  = pnt.ra_deg()
        dec = pnt.dec_deg()

        # Plot pointing
        plt.plot([ra], [dec], 'P', color=colors[i])

    # Return
    return


# ======== #
# Plot all #
# ======== #
def plot_all(filename, smooth=0.02):
    """
    Plot residuals

    Parameters
    ----------
    filename : str
        Sky map file name
    smooth : float, optional
        Map smoothing parameter (deg)
    """
    # Load Crab sky map
    resmap = gammalib.GSkyMap(filename)
    resmap.stack_maps()

    # Create figure
    fig = plt.figure('Map', (10,5))
    fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.10)
    ax = plt.subplot()

    # Plot map
    plot_map(resmap, ax, smooth=smooth)

    # Plot circles
    plot_circles()

    # Plot pointings
    plot_pointings()

    # Show map
    plt.show()

    # Save figure
    filename = os.path.splitext(filename)[0] + '.png'
    fig.savefig(filename,dpi=300)

    # Return
    return


# ======== #
# Show map #
# ======== #
def show_map():
    """
    Show map
    """
    # Plot
    plot_all('crab_resmap.fits', smooth=0.02)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Set filenames
    show_map()
