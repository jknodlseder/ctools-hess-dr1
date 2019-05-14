#! /usr/bin/env python
# ==========================================================================
# Show MSH 15-52 extension
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
import math
import gammalib
import cscripts
import matplotlib.pyplot as plt


# ====================== #
# Generate MSH 15-52 map #
# ====================== #
def generate_map():
    """
    Generate MSH 15-52 map
    """
    # Define files
    inmodel = 'msh_results_egauss_plaw_lookup_grad_hess_edisp.xml'
    logfile = 'msh_skymap.log'
    outfile = 'msh_skymap.fits'
    modfile = 'msh_skymap.xml'

    # Continue only if map or result does not yet exist
    if not os.path.isfile(outfile) and os.path.isfile(inmodel):

        # Read input model
        models = gammalib.GModels(inmodel)

        # Remove MSH 15-52 component
        models.remove('MSH 15-52')

        # Save input model
        models.save(modfile)

        # Generate residual map
        resmap = cscripts.csresmap()
        resmap['inobs']     = 'obs_msh_selected.xml'
        resmap['inmodel']   = modfile
        resmap['outmap']    = outfile
        resmap['emin']      = 0.381
        resmap['emax']      = 40.0
        resmap['enumbins']  = 40
        resmap['nxpix']     = 200
        resmap['nypix']     = 200
        resmap['binsz']     = 0.005
        resmap['coordsys']  = 'CEL'
        resmap['proj']      = 'TAN'
        resmap['xref']      = 228.55
        resmap['yref']      = -59.17
        resmap['algorithm'] = 'SUB'
        resmap['logfile']   = logfile
        resmap['chatter']   = 4
        resmap.execute()

    # Load map if it exists
    if os.path.isfile(outfile):
        map = gammalib.GSkyMap(outfile)
    else:
        map = None

    # Return map
    return map


# ======== #
# Plot map #
# ======== #
def plot_map(inmap, sub, smooth=0.02, ra_offset=228.55, dec_offset=-59.17):
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
    ra_offset : float, optional
        Right Ascension offset (deg)
    dec_offset : float, optional
        Declination offset (deg)
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
                value = 0
            row.append(value)
        array.append(row)

    # Get skymap boundaries
    ra_min   = map.pix2dir(gammalib.GSkyPixel(0.0,map.ny()/2)).ra_deg()
    ra_max   = map.pix2dir(gammalib.GSkyPixel(map.nx(),map.ny()/2)).ra_deg()
    dec_min  = map.pix2dir(gammalib.GSkyPixel(map.nx()/2,0)).dec_deg()
    dec_max  = map.pix2dir(gammalib.GSkyPixel(map.nx()/2,map.ny())).dec_deg()
    ra_min  -= ra_offset
    ra_max  -= ra_offset
    dec_min -= dec_offset
    dec_max -= dec_offset
    aspect   = abs((ra_max-ra_min)/(dec_max-dec_min))

    # Show map
    c = sub.imshow(array, extent=(ra_min,ra_max,dec_min,dec_max),
                   cmap=plt.get_cmap('jet'), aspect=aspect)
    sub.set_xlim([ra_min,ra_max])
    sub.set_ylim([dec_min,dec_max])
    sub.set_xlabel('Right Ascension + %.2f (deg)' % ra_offset, fontsize=14)
    sub.set_ylabel('Declination - %.2f (deg)' % abs(dec_offset), fontsize=14)

    # Set ticks
    for tick in sub.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in sub.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    sub.xaxis.set_ticks_position('both')
    sub.yaxis.set_ticks_position('both')

    # Return
    return


# ============ #
# Plot ellipse #
# ============ #
def plot_ellipse(sub, ra, dec, a, b, pa, label='_nolabel_', e_ra=0, e_dec=0,
                 color='w', fmt='-', linewidth=2,
                 ra_offset=228.55, dec_offset=-59.17):
    """
    Plot ellipse

    Parameters
    ----------
    sub : pyplot
        Frame for map
    ra : float
        Right Ascension of centre (deg)
    dec : float
        Declination of centre (deg)
    a : float
        Major axis (deg)
    b : float
        Minor axis (deg)
    pa : float
        Position angle (deg)
    label : str, optional
        Label string
    e_ra : float, optional
        Error of Right Ascension of centre (deg)
    e_dec : float, optional
        Error of Declination of centre (deg)
    color : str, optional
        Color string
    fmt : str, optional
        Format string
    linewidth : int, optional
        Line width
    ra_offset : float, optional
        Right Ascension offset (deg)
    dec_offset : float, optional
        Declination offset (deg)
    """
    # Get RA normalization (stretch)
    ra_norm = 1.0 / math.cos(59.14*gammalib.deg2rad)

    # Create arrays
    sin_rad = [a*b/(math.sqrt(a*a*math.sin(angle*gammalib.deg2rad)**2 +
                              b*b*math.cos(angle*gammalib.deg2rad)**2)) *
               math.sin((angle+pa)*gammalib.deg2rad) for angle in range(360)]
    cos_rad = [a*b/(math.sqrt(a*a*math.sin(angle*gammalib.deg2rad)**2 +
                              b*b*math.cos(angle*gammalib.deg2rad)**2)) *
               math.cos((angle+pa)*gammalib.deg2rad) for angle in range(360)]
    x       = [ra  - ra_offset  + sin_r*ra_norm for sin_r in sin_rad]
    y       = [dec - dec_offset + cos_r for cos_r in cos_rad]

    # Plot circle
    sub.plot(x, y, fmt, color=color, linewidth=linewidth, label=label)

    # Plot RA error bar
    if e_ra > 0.0:
        x = [ra-ra_offset - e_ra, ra-ra_offset +e_ra]
        y = [dec-dec_offset,      dec-dec_offset]
        sub.plot(x, y, fmt, color=color, linewidth=1)

    # Plot DEC error bar
    if e_dec > 0.0:
        x = [ra-ra_offset,         ra-ra_offset]
        y = [dec-dec_offset-e_dec, dec-dec_offset+e_dec]
        sub.plot(x, y, fmt, color=color, linewidth=1)

    # Return
    return


# ======================== #
# Show MSH 15-52 extension #
# ======================== #
def show_map_extension():
    """
    Show MSH 15-52 extension
    """
    # Generate map
    map = generate_map()

    # Continue only if we have a map
    if map != None:

        # Create figure
        fig = plt.figure('Map', (7, 7))
        fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.10)
        ax = plt.subplot()

        # Plot map
        plot_map(map, ax)

        # Plot ctools ellipse
        plot_ellipse(ax, 228.547875537934, -59.1744150877659,
                     0.105545452579442, 0.0617228294639243, 57.2950204130012+90.0,
                     e_ra=0.0116431938530672, e_dec=0.00689124693198501,
                     color='red', label='ctools')

        # Plot Aharonian et al. (2005) ellipse
        plot_ellipse(ax, 228.52917, -59.1575,
                     0.106667, 0.038333, 139.0,
                     e_ra=0.005833, e_dec=0.0030556, color='blue',
                     label='Aharonian et al. (2005)')

        # Plot On region
        plot_ellipse(ax, 228.547, -59.174,
                     0.2, 0.2, 0.0,
                     color='white',
                     label='_nolegend_')

        # Add legend
        ax.legend()

        # Plot figure
        plt.show()

        # Save figure
        fig.savefig('msh_skymap.png',dpi=300)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Show MSH 15-52 extension
    show_map_extension()
