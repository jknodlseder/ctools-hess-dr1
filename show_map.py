#! /usr/bin/env python
# ==========================================================================
# Show sky map
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
import os
import sys
import math
import gammalib
import ctools
import cscripts
import matplotlib.pyplot as plt


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
    # Set radius size in longitude and latitude
    dx = rad/math.cos(dec*gammalib.deg2rad)
    dy = rad

    # Create arrays
    x = [ra +dx*math.sin(angle*gammalib.deg2rad) for angle in range(360)]
    y = [dec+dy*math.cos(angle*gammalib.deg2rad) for angle in range(360)]

    # Plot circle
    plt.plot(x, y, fmt, color=color, linewidth=linewidth)

    # Return
    return


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
    c    = sub.imshow(array, extent=(ra_min,ra_max,dec_min,dec_max),
                      cmap=plt.get_cmap('jet'), aspect=aspect)
    cbar = plt.colorbar(c, orientation='vertical', shrink=0.7)
    sub.set_xlabel('Right Ascension (deg)', fontsize=14)
    sub.set_ylabel('Declination (deg)', fontsize=14)

    # Return
    return


# =========================== #
# Determine spatial residuals #
# =========================== #
def residual_spatial(cntcube, modcube):
    """
    Determine spatial residuals

    Parameters
    ----------
    cntcube : str
        Counts cube name
    modcube : str
        Model cube name

    Returns
    -------
    map, error : `~gammalib.GSkyMap()`
        Residual and error sky maps
    """
    # Load counts and model cubes
    cnt = gammalib.GSkyMap(cntcube)
    mod = gammalib.GSkyMap(modcube)

    # Generate residual map
    map = cnt - mod
    map.stack_maps()

    # Generate residual error bars
    error = mod.copy()
    error.stack_maps()
    error = error.sqrt()

    # Return residual map and error map
    return map, error


# ======== #
# Plot all #
# ======== #
def plot_all(cntcube, modcube, smooth=0.02):
    """
    Plot residuals

    Parameters
    ----------
    cntcube : str
        Counts cube name
    modcube : str
        Model cube name
    smooth : float, optional
        Smoothing width (deg)
    """
    # Get map
    try:
        resmap, _ = residual_spatial(cntcube, modcube)
    except:
        resmap = gammalib.GSkyMap(modcube)
        resmap.stack_maps()

    # Create figure
    fig = plt.figure('Map', (7, 7))
    fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.10)
    ax = plt.subplot()

    # Plot map
    plot_map(resmap, ax, smooth=smooth)

    # Plot circles
    #plot_circle(258.1125, -39.6867, 0.6, color='w', fmt='-', linewidth=1)
    #plot_circle(258.1125, -39.6867, 0.8, color='w', fmt='--', linewidth=1)
    #plot_circle(258.1125, -39.6867, 1.0, color='w', fmt='--', linewidth=1)

    # Show map
    plt.show()

    # Save figure
    filename = os.path.splitext(modcube)[0] + '.png'
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
    # Set usage string
    usage = 'show_map.py [-s smooth] [-c cntcube] file'

    # Set default options
    options = [{'option': '-s', 'value': 0.0},
               {'option': '-c', 'value': 'cntcube.fits'}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    smooth  = float(options[0]['value'])
    cntcube = options[1]['value']

    # Plot
    plot_all(cntcube, args[0], smooth=smooth)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Set filenames
    show_map()
