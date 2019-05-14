#! /usr/bin/env python
# ==========================================================================
# Show Crab extension
#
# Copyright (C) 2018-2019 Juergen Knoedlseder
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


# ======== #
# Plot map #
# ======== #
def plot_map(inmap, sub, smooth=0.0):
    """
    Plot sky map

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
        map.smooth('GAUSSIAN', smooth)
    else:
        map = inmap
    
    # Create array from skymap
    array = []
    v_max = 0.0
    for iy in range(map.ny()):
        row = []
        for ix in range(map.nx()):
            index = ix+(map.ny()-iy-1)*map.nx()
            value = map[index] - 10
            if value < 0:
                value = 0
            row.append(value)
        array.append(row)

    # Get skymap boundaries
    ra_min  = map.pix2dir(gammalib.GSkyPixel(0.0,map.ny()/2)).ra_deg()
    ra_max  = map.pix2dir(gammalib.GSkyPixel(map.nx(),map.ny()/2)).ra_deg()
    dec_min = map.pix2dir(gammalib.GSkyPixel(map.nx()/2,0)).dec_deg()
    dec_max = map.pix2dir(gammalib.GSkyPixel(map.nx()/2,map.ny())).dec_deg()
    ra_min  -= 83.63
    ra_max  -= 83.63
    dec_min -= 22.015
    dec_max -= 22.015
    aspect  = abs((ra_max-ra_min)/(dec_max-dec_min))

    # Show Aitoff projection
    c = sub.imshow(array, extent=(ra_min,ra_max,dec_min,dec_max),
                   cmap=plt.get_cmap('jet'), aspect=aspect)
    sub.set_xlim([ra_min,ra_max])
    sub.set_ylim([dec_min,dec_max])
    sub.set_xlabel('Right Ascension + 83.63 (deg)', fontsize=14)
    sub.set_ylabel('Declination + 22.015 (deg)', fontsize=14)

    # Set ticks
    for tick in sub.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in sub.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    sub.xaxis.set_ticks_position('both')
    sub.yaxis.set_ticks_position('both')

    # Return
    return


# =========== #
# Plot circle #
# =========== #
def plot_circle(sub, ra, dec, rad, e_ra=0, e_dec=0, color='w', fmt='-', linewidth=3):
    """
    Plot circle

    Parameters
    ----------
    sub : pyplot
        Frame for map
    ra : float
        Right Ascension (deg)
    dec : float
        Declination (deg)
    rad : float
        Radius (deg)
    e_ra : float, optional
        Uncertainty in Right Ascension (deg)
    e_dec : float, optional
        Uncertainty in Declination (deg)
    color : str, optional
        Color string
    fmt : str, optional
        Format string
    linewidth : int, optional
        Line width
    """
    # Create arrays
    x = [ra-83.63  +rad*math.sin(angle*gammalib.deg2rad) for angle in range(360)]
    y = [dec-22.015+rad*math.cos(angle*gammalib.deg2rad) for angle in range(360)]

    # Plot circle
    sub.plot(x, y, fmt, color=color, linewidth=linewidth)

    # Plot RA error bar
    if e_ra > 0.0:
        x = [ra-83.63-e_ra, ra-83.63+e_ra]
        y = [dec-22.015,    dec-22.015]
        sub.plot(x, y, fmt, color=color, linewidth=linewidth)

    # Plot DEC error bar
    if e_dec > 0.0:
        x = [ra-83.63,         ra-83.63]
        y = [dec-22.015-e_dec, dec-22.015+e_dec]
        sub.plot(x, y, fmt, color=color, linewidth=linewidth)

    # Return
    return


# ======== #
# Plot all #
# ======== #
def plot_all(crabmap, smooth=0.0):
    """
    Plot residuals

    Parameters
    ----------
    crabmap : str
        Crab map file name
    smooth : float, optional
        Smoothing width (deg)
    """
    # Set parameters
    size = 1000

    # Load Crab map
    chandra = gammalib.GSkyMap(crabmap)

    # Create a slighly larger map
    crab = gammalib.GSkyMap('TAN', 'CEL', 83.63, 22.015, -6.833e-05, 6.833e-05, size, size)
    crab += chandra

    # Create figure
    fig = plt.figure('Map', (7, 7))
    fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.10)
    ax = plt.subplot()

    # Plot map
    plot_map(crab, ax, smooth=smooth)

    # Plot Holler et al. 2017 circle
    plot_circle(ax, 83.62875, 22.012361, 0.0145, e_ra=0.00033, e_dec=0.0003055,
                color='white', linewidth=4)
    plot_circle(ax, 83.62875, 22.012361, 0.0145, e_ra=0.00033, e_dec=0.0003055,
                color='deepskyblue')
    ax.text(-0.005, -0.02, 'Holler et al. (2017)', fontsize=14, fontweight='bold',
            color='deepskyblue')

    # Plot ctools circle
    # RA .......................: 83.6206711191164 +/- 0.00293905226757571 deg
    # DEC ......................: 22.0249161487514 +/- 0.00271599100158073 deg
    # Sigma ....................: 0.0131191472448878 +/- 0.00513796531189789 deg
    plot_circle(ax, 83.620671, 22.024916, 0.013119, e_ra=0.002939, e_dec=0.002716,
                color='white', linewidth=4)
    plot_circle(ax, 83.620671, 22.024916, 0.013119, e_ra=0.002939, e_dec=0.002716,
                color='red')
    ax.text(-0.023, 0.021, 'ctools', fontsize=14, fontweight='bold',
            color='salmon')

    # Plot systematic uncertainty
    plot_circle(ax, 83.63+0.028, 22.015+0.025, 0.005555, color='white')
    ax.text(0.023, 0.029, 'systematic uncertainty', fontsize=12, color='white')

    # Save figure
    plt.show()
    fig.savefig('crab_extension.png',dpi=300)

    # Return
    return


# =================== #
# Show Crab extension #
# =================== #
def show_crab_extension():
    """
    Show Crab extension
    """
    # Set usage string
    usage = 'show_map.py [-s smooth] [-crab map]'

    # Set default options
    options = [{'option': '-s', 'value': 0.0},
               {'option': '-crab', 'value': '../chandra_crab.fits'}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    smooth  = float(options[0]['value'])
    crabmap = options[1]['value']

    # Plot
    plot_all(crabmap, smooth=smooth)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show Crab extension
    show_crab_extension()
