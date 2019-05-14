#! /usr/bin/env python
# ==========================================================================
# Show PKS 2155-304 position
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
            value = map[index]
            if value < 0:
                value = 0.0
            row.append(math.sqrt(value))
        array.append(row)

    # Get skymap boundaries
    ra_min  = map.pix2dir(gammalib.GSkyPixel(0.0,map.ny()/2)).ra_deg()
    ra_max  = map.pix2dir(gammalib.GSkyPixel(map.nx(),map.ny()/2)).ra_deg()
    dec_min = map.pix2dir(gammalib.GSkyPixel(map.nx()/2,0)).dec_deg()
    dec_max = map.pix2dir(gammalib.GSkyPixel(map.nx()/2,map.ny())).dec_deg()
    ra_min  -= 329.71694
    ra_max  -= 329.71694
    dec_min -= -30.22559
    dec_max -= -30.22559
    aspect  = abs((ra_max-ra_min)/(dec_max-dec_min))

    # Show Aitoff projection
    c = sub.imshow(array, extent=(ra_min,ra_max,dec_min,dec_max),
                   cmap=plt.get_cmap('jet'), aspect=aspect)
    sub.set_xlim([ra_min,ra_max])
    sub.set_ylim([dec_min,dec_max])
    sub.set_xlabel('Right Ascension + 329.71694 (deg)', fontsize=14)
    sub.set_ylabel('Declination - 30.22559 (deg)', fontsize=14)

    # Set ticks
    for tick in sub.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in sub.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    sub.xaxis.set_ticks_position('both')
    sub.yaxis.set_ticks_position('both')

    # Return
    return


# ================== #
# Plot error ellipse #
# ================== #
def plot_error_ellipse(sub, ra, dec, e_ra, e_dec, label='_nolabel_',
                       color='w', fmt='-', linewidth=2):
    """
    Plot error ellipse

    Parameters
    ----------
    sub : pyplot
        Frame for error ellipse
    ra : float
        Right Ascension of centre (deg)
    dec : float
        Declination of centre (deg)
    e_ra : float
        Error in Right Ascension of centre (deg)
    e_dec : float
        Error in Declination of centre (deg)
    label : string, optional
        Ellipse label
    color : string, optional
        Ellipse color
    fmt : string, optional
        Ellipse format
    linewidth : int, optional
        Ellipse line width
    """
    # Get RA normalization (stretch)
    ra_norm = 1.0 / math.cos(30.22559*gammalib.deg2rad)

    # Set ellipse parameters
    a  = e_ra
    b  = e_dec
    pa = 90.0

    # Create arrays
    sin_rad = [a*b/(math.sqrt(a*a*math.sin(angle*gammalib.deg2rad)**2 +
                              b*b*math.cos(angle*gammalib.deg2rad)**2)) *
               math.sin((angle+pa)*gammalib.deg2rad) for angle in range(360)]
    cos_rad = [a*b/(math.sqrt(a*a*math.sin(angle*gammalib.deg2rad)**2 +
                              b*b*math.cos(angle*gammalib.deg2rad)**2)) *
               math.cos((angle+pa)*gammalib.deg2rad) for angle in range(360)]
    x       = [ra  - 329.71694 + sin_r for sin_r in sin_rad]
    y       = [dec +  30.22559 + cos_r for cos_r in cos_rad]

    # Plot ellipse
    sub.plot(x, y, fmt, color=color, linewidth=linewidth, label=label)

    # Compute distance
    fitted   = gammalib.GSkyDir()
    position = gammalib.GSkyDir()
    hess = gammalib.GSkyDir()
    fitted.radec_deg(ra, dec)
    position.radec_deg(329.7169, -30.2256)
    hess.radec_deg(329.7192, -30.2249)
    dist  = fitted.dist_deg(position)
    dhess = fitted.dist_deg(hess)
    print('%s distance from true position: %.3f' % (label, dist*3600.0))
    print('%s distance from H.E.S.S. position: %.3f' % (label, dhess*3600.0))

    # Return
    return


# ======== #
# Plot all #
# ======== #
def plot_all(pksmap, size=500, smooth=0.0):
    """
    Plot everything

    Parameters
    ----------
    pksmap : str
        PKS sky map
    size : int, optional
        Number of bins of sky map
    smooth : float, optional
        Width of smoothing kernel (deg)
    """
    # Load PKS 2155-304 map
    panstarrs = gammalib.GSkyMap(pksmap)

    # Create map and merge in PKS 2155-304 map
    pks = gammalib.GSkyMap('TAN', 'CEL', 329.71694, -30.22559, -1.5e-5, 1.5e-5,
                           size, size)
    pks += panstarrs

    # Create figure
    fig = plt.figure('Map', (8,7))
    fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.10)
    ax = plt.subplot()

    # Plot map
    plot_map(pks, ax, smooth=smooth)

    # Plot ctools error ellipse (curved power law)
    plot_error_ellipse(ax, 329.718468529345, -30.2236388733347,
                       0.000757258427185772, 0.00066055030359906,
                       color='red', label='ctools (T200)')
    plot_error_ellipse(ax, 329.718438975006, -30.2239335046493,
                       0.00071405381958194, 0.000619571439060395,
                       color='red', fmt='--', label='ctools (T300)')
    plot_error_ellipse(ax, 329.717902296568, -30.2241135104247,
                       0.00106581586903552, 0.000924063456354341,
                       color='red', fmt=':', label='ctools (T700)')

    # Plot Aharonian et al. 2009 error ellipse
    plot_error_ellipse(ax, 329.71917, -30.224944, 0.0004166, 0.0005277,
                       color='deepskyblue', label='Aharonian et al. (2009)')

    # Add 20" line
    x = [0.0018, 0.0018 - 20.0/3600.0]
    y = [0.0032, 0.0032]
    ax.plot(x, y, '-', color='white', linewidth=2)
    ax.text(-0.0008, 0.0028, '20"', fontsize=14, color='white')

    # Add legend
    ax.legend(loc='lower right')

    # Save figure
    plt.show()
    fig.savefig('pks_position.png',dpi=300)

    # Return
    return


# ========================== #
# Show PKS 2155-304 position #
# ========================== #
def show_pks_positions():
    """
    Show PKS 2155-304 position
    """
    # Set usage string
    usage = 'show_pks_position.py [-s smooth] [-pks map]'

    # Set default options
    options = [{'option': '-s', 'value': 0.0},
               {'option': '-pks', 'value': '$CTAGITROOT/analysis/hess_dr1/pks_PanSTARRS_y.fits'}]

    # Get arguments and options from command line arguments
    args, options = cscripts.ioutils.get_args_options(options, usage)

    # Extract script parameters from options
    smooth = float(options[0]['value'])
    pksmap = options[1]['value']

    # Plot
    plot_all(pksmap, smooth=smooth)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show PKS 2155-304 position
    show_pks_positions()
