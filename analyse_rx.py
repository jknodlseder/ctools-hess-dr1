#! /usr/bin/env python
# ==========================================================================
# Analyse RX J1713.7-3946 data
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
import math
import gammalib
import ctools
import cscripts


# ================== #
# Create region maps #
# ================== #
def create_region_maps(rad_on=0.6, rad_off=[0.8,1.0]):
    """
    Create region maps

    Parameters
    ----------
    rad_on : float, optional
        On region radius (deg)
    rad_off : list of float, optional
        Off region minimum and maximum radius (deg)
    """
    # Set centre
    centre = gammalib.GSkyDir()
    centre.radec_deg(258.1125, -39.6867)

    # Create sky maps
    srcreg = gammalib.GSkyMap('TAN', 'CEL', 258.1125, -39.6867, 0.01, 0.01, 300, 300)
    bkgreg = gammalib.GSkyMap('TAN', 'CEL', 258.1125, -39.6867, 0.01, 0.01, 300, 300)

    # Fill sky map
    for i in range(srcreg.npix()):
        rad = centre.dist_deg(srcreg.inx2dir(i))
        if rad <= rad_on:
            srcreg[i] = 1
        if rad > rad_off[0] and rad < rad_off[1]:
            bkgreg[i] = 1

    # Save maps
    srcreg.save('rx_srcreg_map.fits', True)
    bkgreg.save('rx_bkgreg_map.fits', True)

    # Return
    return


# ============================= #
# Add attribute to XML filename #
# ============================= #
def add_attribute(filename, attribute):
    """
    Add attribute to XML filename

    Parameters
    ----------
    filename : str
        XML filename
    attribute : str
        Attribute

    Returns
    -------
    filename : str
        XML filename with attribute
    """
    # Add attribute to XML filename
    fname = os.path.splitext(filename)[0] + attribute + os.path.splitext(filename)[1]

    # Return filename
    return fname


# =================== #
# Scaled template map #
# =================== #
def scale_map(map, alpha):
    """
    Create scaled template map

    Parameters
    ----------
    map : `~gammalib.GSkyMap()`
        Original sky map
    alpha : float
        Scaling parameter

    Returns
    -------
    map : `~gammalib.GSkyMap()`
        Scaled template map
    """
    # Create copy of map
    scaled_map = map.copy()

    # Scale all pixels
    for i in range(scaled_map.npix()):
        if scaled_map[i] < 0.0:
            scaled_map[i] = 0.0
        scaled_map[i] = math.pow(scaled_map[i], alpha)

    # Return scale map
    return scaled_map


# ================ #
# Set source model #
# ================ #
def set_source_model(srcmodel, spec='plaw', alpha=1.0, mapname='../map_RXJ1713.fits'):
    """
    Set source model

    Parameters
    ----------
    srcmodel : str
        Source model name
    spec : str, optional
        Spectral model
    alpha : float, optional
        Map scaling factor
    mapname : str, optional
        Sky map name

    Returns
    -------
    source : `~gammalib.GModelSky()`
        Source model
    """
    # Set spectral component
    if spec == 'plaw':
        spectral = gammalib.GModelSpectralPlaw(1.0e-17, -2.0,
                                               gammalib.GEnergy(1.0,'TeV'))
        spectral['Prefactor'].min(1.0e-25)
    elif spec == 'eplaw':
        spectral = gammalib.GModelSpectralExpPlaw(1.0e-17, -2.0,
                                                  gammalib.GEnergy(1.0,'TeV'),
                                                  gammalib.GEnergy(10.0,'TeV'))
        spectral['CutoffEnergy'].min(1.0e6)
        spectral['CutoffEnergy'].max(1.0e8)
        spectral['Prefactor'].min(1.0e-25)
    elif spec == 'inveplaw':
        spectral = gammalib.GModelSpectralExpInvPlaw(1.0e-17, -2.0,
                                                     gammalib.GEnergy(1.0,'TeV'),
                                                     gammalib.GEnergy(10.0,'TeV'))
        spectral['Prefactor'].min(1.0e-25)
    elif spec == 'logparabola':
        spectral = gammalib.GModelSpectralLogParabola(1.0e-17, -2.0,
                                                      gammalib.GEnergy(1.0,'TeV'),
                                                      -0.3)
        spectral['Prefactor'].min(1.0e-25)
    elif spec == 'abdalla2018':
        spectral = gammalib.GModelSpectralExpPlaw(2.3e-17, -2.06,
                                                  gammalib.GEnergy(1.0,'TeV'),
                                                  gammalib.GEnergy(12.9,'TeV'))
        spectral['Prefactor'].fix()
        spectral['Index'].fix()
        spectral['CutoffEnergy'].fix()
    elif spec == 'abdalla2018Ec':
        spectral = gammalib.GModelSpectralExpPlaw(2.3e-17, -2.06,
                                                  gammalib.GEnergy(1.0,'TeV'),
                                                  gammalib.GEnergy(12.9,'TeV'))
        spectral['CutoffEnergy'].fix()
        spectral['Prefactor'].min(1.0e-25)
    elif spec == 'aharonian2007':
        spectral = gammalib.GModelSpectralLogParabola(2.06e-17, -2.02,
                                                      gammalib.GEnergy(1.0,'TeV'),
                                                      -0.29)
        spectral['Prefactor'].fix()
        spectral['Index'].fix()
        spectral['Curvature'].fix()

    # Set spatial component
    if 'map' in srcmodel:
        filename = '%s' % mapname
        map      = gammalib.GSkyMap(filename)
        map      = scale_map(map, alpha)
        filename = 'map_RXJ1713_%.2f.fits' % alpha
        map.save(filename, True)
        spatial = gammalib.GModelSpatialDiffuseMap(filename)
    elif 'disk' in srcmodel:
        dir     = gammalib.GSkyDir()
        dir.radec_deg(258.3, -39.7)
        spatial = gammalib.GModelSpatialRadialDisk(dir, 0.5)
        spatial['RA'].free()
        spatial['DEC'].free()
    elif 'gauss' in srcmodel:
        dir     = gammalib.GSkyDir()
        dir.radec_deg(258.3, -39.7)
        spatial = gammalib.GModelSpatialRadialGauss(dir, 0.5)
        spatial['RA'].free()
        spatial['DEC'].free()
    elif 'shell' in srcmodel:
        dir     = gammalib.GSkyDir()
        dir.radec_deg(258.3, -39.7)
        spatial = gammalib.GModelSpatialRadialShell(dir, 0.4, 0.2)
        spatial['RA'].free()
        spatial['DEC'].free()

    # Set source model
    source = gammalib.GModelSky(spatial, spectral)
    source.name('RX J1713.7-3946')
    source.tscalc(True)

    # Return source model
    return source


# ========= #
# Set model #
# ========= #
def set_model(srcmodel, spec, bkgname, alpha=1.0):
    """
    Set model

    Parameters
    ----------
    srcmodel : str
        Source model name
    spec : str
        Spectral model
    bkgname : str
        Background filename
    alpha : float, optional
        Map scaling factor

    Returns
    -------
    models : `~gammalib.GModels()`
        Models
    """
    # Initialise model container
    models = gammalib.GModels()

    # Set source model
    if srcmodel != '':
        source = set_source_model(srcmodel, spec=spec, alpha=alpha)
        if not models.contains(source.name()):
            models.append(source)

    # Append background models
    if bkgname != '':
        bkg_models = gammalib.GModels(bkgname)
        for bkg_model in bkg_models:
            if not models.contains(bkg_model.name()):
                models.append(bkg_model)

    # Return models
    return models


# =================== #
# Set analysis string #
# =================== #
def set_analysis(srcmodel, bkgname, spec, pars, alpha=1.0, select=None):
    """
    Set analysis string

    Parameters
    ----------
    srcmodel : str
        Source model
    bkgname : str
        Background name
    spec : str
        Spectral model
    pars : dict
        Analysis parameters
    alpha : float, optional
        Map scaling factor
    select : list of int, optional
        Indices for observation selection
    """
    # Set source model
    if srcmodel == 'map':
        if alpha == 1.0:
            _srcmodel = srcmodel
        else:
            _srcmodel = 'map%.2f' % alpha
    else:
        _srcmodel = srcmodel

    # Set analysis
    analysis = '%s_%s_lookup_grad_hess' % (_srcmodel, spec)
    if select != None:
        analysis += '_select'
        for i in select:
            analysis += str(i)
    if pars['onoff']:
        if pars['stacked']:
            analysis += '_onoff-stacked%2.2d' % pars['ebins']
        else:
            analysis += '_onoff-joint%2.2d' % pars['ebins']
        if bkgname == '':
            if pars['use_cstat_off_weight']:
                analysis += '_wstat-cstat-off-weight'
            else:
                analysis += '_wstat'
        else:
            analysis += '_cstat'
    elif pars['stacked']:
        analysis += '_stacked%2.2d' % pars['ebins']
    elif pars['binned']:
        analysis += '_binned%2.2d' % pars['ebins']
    if pars['edisp']:
        analysis += '_edisp'

    # Return
    return analysis


# ==================== #
# Set observation name #
# ==================== #
def set_observation(obsname, srcmodel, bkgname, pars, alpha=1.0):
    """
    Set observation name

    Parameters
    ----------
    obsname : str
        Observation definition file name
    srcmodel : str
        Source model (only used for On/Off)
    bkgname : str
        Background name
    pars : dict
        Analysis parameters
    alpha : float, optional
        Map scaling factor
    """
    # Set source model
    if srcmodel == 'map':
        if alpha == 1.0:
            _srcmodel = srcmodel
        else:
            _srcmodel = 'map%.2f' % alpha
    else:
        _srcmodel = srcmodel

    # Set observation name
    if pars['onoff']:
        if bkgname == '':
            if pars['use_cstat_off_weight']:
                statistic = 'cstat' # Use cstat file
            else:
                statistic = 'wstat'
        else:
            statistic = 'cstat'
        if pars['stacked']:
            _obsname = add_attribute(obsname, '_onoff-stacked%2.2d_%s_%s' %
                                     (pars['ebins'], statistic, _srcmodel))
        else:
            _obsname = add_attribute(obsname, '_onoff-joint%2.2d_%s_%s' %
                                     (pars['ebins'], statistic, _srcmodel))
    elif pars['stacked']:
        _obsname = add_attribute(obsname, '_stacked%2.2d' % pars['ebins'])
        if pars['edisp']:
            _obsname = add_attribute(_obsname, '_edisp')
    elif pars['binned']:
        _obsname = add_attribute(obsname, '_binned%2.2d' % pars['ebins'])
    else:
        _obsname = obsname

    # Return
    return _obsname


# ======================= #
# Set analysis parameters #
# ======================= #
def set_pars(stacked=False, binned=False, onoff=False, ebins=20, edisp=False,
             use_cstat_off_weight=False):
    """
    Set analysis parameters

    Parameters
    ----------
    stacked : boolean, optional
        Perform stacked analysis
    binned : boolean, optional
        Perform joint binned analysis
    onoff : boolean, optional
        Perform On/Off analysis
    ebins : int, optional
        Number of energy bins for binned or stacked analysis
    edisp : boolean, optional
        Use energy dispersion
    use_cstat_off_weight : boolean, optional
        Use cstat off weights
    """
    # Set analysis parameters
    pars = {'stacked':              stacked,
            'binned' :              binned,
            'onoff':                onoff,
            'ebins':                ebins,
            'edisp':                edisp,
            'use_cstat_off_weight': use_cstat_off_weight}

    # Return
    return pars


# ================ #
# Fit observations #
# ================ #
def fit(obs, srcmodel, spec, bkgname, analysis, pars, alpha=1.0):
    """
    Analyse one observation and determine the background model

    Parameters
    ----------
    obs : `~gammalib.GObservations()`
        Observation container
    srcmodel : str
        Source model
    spec : str
        Spectral model
    bkgname : str
        Background model definition XML file
    analysis : str
        Analyse name
    pars : dict
        Dictionary of analysis parameters
    alpha : float, optional
        Map scaling factor
    """
    # Set file names
    outmodel = 'rx_results_%s.xml' % analysis
    logfile  = 'rx_results_%s.log' % analysis

    # Continue only if result file does not exist
    if not os.path.isfile(outmodel):

        # Set models
        models = set_model(srcmodel, spec, bkgname, alpha=alpha)

        # Attach models
        obs.models(models)

        # Perform maximum likelihood fitting with initial models
        like = ctools.ctlike(obs)
        like['edisp']    = pars['edisp']
        like['outmodel'] = outmodel
        like['max_iter'] = 500
        like['logfile']  = logfile
        like['debug']    = True
        like.logFileOpen()

        # Use wstat if no background model is provided
        if bkgname == '':
            like['statistic'] = 'WSTAT'

        # Execute ctlike
        like.execute()

    # Return
    return


# ========================= #
# Determine Npred of source #
# ========================= #
def npred(obs, analysis, pars):
    """
    Determine Npred of source

    Parameters
    ----------
    obs : `~gammalib.GObservations()`
        Observation container
    analysis : str
        Analyse name
    pars : dict
        Dictionary of analysis parameters
    """
    # Set file names
    inmodel  = 'rx_results_%s.xml' % analysis
    outmodel = 'rx_npred_%s.xml' % analysis
    logfile  = 'rx_npred_%s.log' % analysis

    # Continue only if result file does not exist
    if not os.path.isfile(outmodel) and os.path.isfile(inmodel):

        # Set models
        models = gammalib.GModels(inmodel)
        model  = models['RX J1713.7-3946']
        for par in model:
            par.fix()
        model.tscalc(False) # Otherwise leads to an error
        npred_models = gammalib.GModels()
        npred_models.append(model)

        # Attach models
        obs.models(npred_models)

        # If statistic is wstat then switch to cstat (since wstat does not
        # correctly compute Npred)
        for o in obs:
            if o.statistic() == 'wstat':
                o.statistic('cstat')

        # Perform maximum likelihood fitting
        like = ctools.ctlike(obs)
        like['edisp']    = pars['edisp']
        like['outmodel'] = outmodel
        like['logfile']  = logfile
        like['debug']    = True
        like.logFileOpen()
        like.execute()

    # Return
    return


# ==================== #
# Analyse observations #
# ==================== #
def analyse(obsname, bkgname, srcmodel, spec, pars, alphas=[1.0]):
    """
    Analyse observations

    Parameters
    ----------
    obsname : str
        Observation definition XML file
    bkgname : str
        Background model definition XML file
    srcmodels : list of str
        List of source models
    spec : str
        Spectral model
    pars : dict
        Dictionary of analysis parameters
    alphas : list of float, optional
        List of template map scaling factors
    """
    # Handle map
    if srcmodel == 'map':

        # Loop over alphas
        for alpha in alphas:

            # Load observations
            _obsname = set_observation(obsname, srcmodel, bkgname, pars, alpha=alpha)
            obs = gammalib.GObservations(_obsname)

            # Set analysis
            analysis = set_analysis(srcmodel, bkgname, spec, pars, alpha=alpha)

            # Fit model
            fit(obs, srcmodel, spec, bkgname, analysis, pars, alpha=alpha)

            # Determine Npred
            npred(obs, analysis, pars)

    # Handle other models
    else:

        # Load observations
        _obsname = set_observation(obsname, srcmodel, bkgname, pars, alpha=1.0)
        obs = gammalib.GObservations(_obsname)

        # Set analysis
        analysis = set_analysis(srcmodel, bkgname, spec, pars, alpha=1.0)

        # Fit model
        fit(obs, srcmodel, spec, bkgname, analysis, pars, alpha=1.0)

        # Determine Npred
        npred(obs, analysis, pars)

    # Return
    return


# ================== #
# Generate butterfly #
# ================== #
def butterfly(obsname, bkgname, srcmodel, spec, pars, alpha=1.0):
    """
    Generate butterfly

    Parameters
    ----------
    obsname : str
        Observation definition XML file
    bkgname : str
        Background model definition XML file
    srcmodel : str
        Source model
    spec : str
        Spectral model
    pars : dict
        Dictionary of analysis parameters
    alpha : float, optional
        Template map scaling factor
    """
    # Set observation name
    _obsname = set_observation(obsname, srcmodel, bkgname, pars)

    # Set analysis
    analysis = set_analysis(srcmodel, bkgname, spec, pars, alpha=alpha)

    # Set file names
    inmodel = 'rx_results_%s.xml' % analysis
    outfile = 'rx_butterfly_%s.txt' % analysis
    logfile = 'rx_butterfly_%s.log' % analysis

    # Continue only if result file does not exist
    if not os.path.isfile(outfile) and os.path.isfile(inmodel):

        # Generate butterfly
        butterfly = ctools.ctbutterfly()
        butterfly['inobs']    = _obsname
        butterfly['inmodel']  = inmodel
        butterfly['srcname']  = 'RX J1713.7-3946'
        butterfly['outfile']  = outfile
        butterfly['edisp']    = pars['edisp']
        butterfly['emin']     = 0.3
        butterfly['emax']     = 50.0
        butterfly['logfile']  = logfile
        butterfly['debug']    = True
        butterfly.logFileOpen()
        butterfly.execute()

    # Return
    return


# ================= #
# Generate spectrum #
# ================= #
def spectrum(obsname, bkgname, srcmodel, pars, emin=0.3, emax=20.0, alpha=1.0):
    """
    Generate spectrum

    Parameters
    ----------
    obsname : str
        Observation definition XML file
    bkgname : str
        Background model definition XML file
    srcmodel : str
        Source model
    pars : dict
        Dictionary of analysis parameters
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    alpha : float, optional
        Template map scaling factor
    """
    # Set observation name
    _obsname = set_observation(obsname, srcmodel, bkgname, pars)

    # Set analysis
    analysis = set_analysis(srcmodel, bkgname, 'plaw', pars, alpha=alpha)

    # Set file names
    outfile = 'rx_spectrum_%s.fits' % (analysis)
    logfile = 'rx_spectrum_%s.log'  % (analysis)

    # Continue only if result file does not exist
    if not os.path.isfile(outfile):

        # Load observations
        obs = gammalib.GObservations(_obsname)

        # Set background model for stacked analysis
        if pars['stacked']:
            bname  = add_attribute(bkgname, '_stacked')
            models = set_model(srcmodel, 'plaw', bname, alpha=alpha)
            for model in models:
                if model.type() == 'CTACubeBackground':
                    model.spectral()['Index'].fix()

        # ... otherwise create H.E.S.S. background model with power law
        # spectral component
        else:

            # Setup plaw background model
            bkg = cscripts.csbkgmodel()
            bkg['inobs']      = '$HESSDATA/obs/obs_rx.xml'
            bkg['instrument'] = 'HESS'
            bkg['spatial']    = 'LOOKUP'
            bkg['slufile']    = 'off_lookup.fits'
            bkg['gradient']   = True
            bkg['spectral']   = 'PLAW'
            bkg['runwise']    = True
            bkg['emin']       = emin
            bkg['emax']       = emax
            bkg['rad']        = 2.0
            bkg.run()

            # Initialise model container
            models = gammalib.GModels()

            # Set source model
            source = set_source_model(srcmodel, spec='plaw', alpha=alpha)
            models.append(source)

            # Set background model
            bkg_models = bkg.models()
            for model in bkg_models:
                if model.type() == 'CTABackground':
                    model.spectral()['Index'].fix()
                    model.spatial()[0].fix()  # Normalization
                    model.spatial()[1].free() # DETX
                    model.spatial()[2].free() # DETY

            # Append background models
            models.extend(bkg_models)

        # Append models
        obs.models(models)

        # Generate spectrum
        spec = cscripts.csspec(obs)
        spec['srcname']  = 'RX J1713.7-3946'
        spec['edisp']    = pars['edisp']
        spec['outfile']  = outfile
        spec['method']   = 'SLICE'
        spec['ebinalg']  = 'LOG'
        spec['emin']     = emin
        spec['emax']     = emax
        spec['enumbins'] = 15
        spec['logfile']  = logfile
        spec['debug']    = True
        spec['chatter']  = 4
        spec.logFileOpen()
        spec.execute()

    # Return
    return


# ===================== #
# Generate residual map #
# ===================== #
def generate_resmap(emin=0.3, emax=50.0, ebins=40, edisp=True):
    """
    Generate residual map

    Parameters
    ----------
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    ebins : int, optional
        Number of energy bins
    edisp : bool, optional
        Enable energy dispersion?
    """
    # Set filenames
    inobs    = 'obs_rx_selected.xml'
    inmodel  = 'rx_results_map0.48_eplaw_lookup_grad_hess_edisp.xml'
    outmodel = 'rx_background_map0.48_eplaw_lookup_grad_hess_edisp.xml'
    outmap   = 'rx_resmap.fits'
    logfile  = 'rx_resmap.log'

    # Continue only if inmodel exists
    if os.path.isfile(inobs) and os.path.isfile(inmodel):

        # Create model definition XML file without source
        models = gammalib.GModels(inmodel)
        models.remove('RX J1713.7-3946')
        models.save(outmodel)

        # Setup task parameters
        resmap = cscripts.csresmap()
        resmap['inobs']     = inobs
        resmap['inmodel']   = outmodel
        resmap['edisp']     = edisp
        resmap['algorithm'] = 'SUB'
        resmap['ebinalg']   = 'LOG'
        resmap['emin']      = emin
        resmap['emax']      = emax
        resmap['enumbins']  = ebins
        resmap['coordsys']  = 'CEL'
        resmap['proj']      = 'CAR'
        resmap['xref']      = 258.1125
        resmap['yref']      = -39.6867
        resmap['nxpix']     = 300
        resmap['nypix']     = 300
        resmap['binsz']     = 0.02
        resmap['outmap']    = outmap
        resmap['logfile']   = logfile
        resmap.logFileOpen()

        # Generate residual map
        resmap.execute()

    # Return
    return


# ================ #
# Generate sky map #
# ================ #
def generate_skymap(emin=0.3, emax=50.0):
    """
    Generate sky map

    Parameters
    ----------
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    """
    # Set filenames
    inobs    = 'obs_rx_selected.xml'
    outmap   = 'rx_skymap.fits'
    logfile  = 'rx_skymap.log'

    # Continue only if inmodel exists
    if os.path.isfile(inobs) and not os.path.isfile(outmap):

        # Setup task parameters
        skymap = ctools.ctskymap()
        skymap['inobs']       = inobs
        skymap['emin']        = emin
        skymap['emax']        = emax
        skymap['nxpix']       = 300
        skymap['nypix']       = 300
        skymap['binsz']       = 0.02
        skymap['coordsys']    = 'CEL'
        skymap['proj']        = 'CAR'
        skymap['xref']        = 258.1125
        skymap['yref']        = -39.6867
        skymap['bkgsubtract'] = 'RING'
        skymap['roiradius']   = 0.01
        skymap['inradius']    = 0.8
        skymap['outradius']   = 1.0
        skymap['iterations']  = 0
        skymap['inexclusion'] = 'rx_srcreg_map.fits'
        skymap['outmap']      = outmap
        skymap['logfile']     = logfile
        skymap.logFileOpen()

        # Generate sky map
        skymap.execute()

    # Return
    return


# ========================== #
# Generate background lookup #
# ========================== #
def generate_background_lookup():
    """
    Generate background lookup from empty field observations
    """
    # Set filenames
    obsname  = '$HESSDATA/obs/obs_off.xml'
    filename = 'off_lookup.fits'

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


# ==================== #
# Prepare observations #
# ==================== #
def prepare(obsname, bkgname, rad=2.0, emin=0.3, emax=50.0, ebins=20):
    """
    Prepare events for analysis

    Parameters
    ----------
    obsname : str
        Observation definition XML file
    bkgname : str
        Background model definition XML file
    rad : float, optional
        Selection radius (degrees)
    emin : float, optional
        Minimum energy for analysis (TeV)
    emax : float, optional
        Maximum energy for analysis (TeV)
    ebins : int, optional
        Number of energy bins
    """
    # Set filenames
    cntcube               = 'rx_stacked%2.2d_cntcube.fits'   % ebins
    expcube               = 'rx_stacked%2.2d_expcube.fits'   % ebins
    psfcube               = 'rx_stacked%2.2d_psfcube.fits'   % ebins
    edispcube             = 'rx_stacked%2.2d_edispcube.fits' % ebins
    bkgcube               = 'rx_stacked%2.2d_bkgcube.fits'   % ebins
    obsname_binned        = add_attribute(obsname, '_binned%2.2d' % ebins)
    bkgname_stacked       = add_attribute(bkgname, '_stacked')
    obsname_stacked       = add_attribute(obsname, '_stacked%2.2d' % ebins)
    obsname_stacked_edisp = add_attribute(obsname_stacked, '_edisp')

    # Generate background lookup
    generate_background_lookup()

    # Continue only if selected events do not exist
    if not os.path.isfile(obsname):

        # Setup task parameters
        select = ctools.ctselect()
        select['inobs']    = '$HESSDATA/obs/obs_rx.xml'
        select['outobs']   = obsname
        select['ra']       = 'UNDEF'
        select['dec']      = 'UNDEF'
        select['rad']      = rad
        select['tmin']     = 'UNDEF'
        select['tmax']     = 'UNDEF'
        select['emin']     = emin
        select['emax']     = emax
        select['usethres'] = 'DEFAULT'
        select['logfile']  = 'rx_hess_select_events.log'
        select.logFileOpen()

        # Select events
        select.execute()

    # Continue only if background model does not exist
    if not os.path.isfile(bkgname):

        # Setup task parameters
        bkg = cscripts.csbkgmodel()
        bkg['inobs']      = '$HESSDATA/obs/obs_rx.xml'
        bkg['outmodel']   = bkgname
        bkg['instrument'] = 'HESS'
        bkg['spatial']    = 'LOOKUP'
        bkg['slufile']    = 'off_lookup.fits'
        bkg['gradient']   = True
        bkg['spectral']   = 'NODES'
        bkg['ebinalg']    = 'LOG'
        bkg['emin']       = emin
        bkg['emax']       = 30.0
        bkg['enumbins']   = 8
        bkg['runwise']    = True
        bkg['rad']        = rad
        bkg['logfile']    = 'rx_hess_create_background.log'
        bkg.logFileOpen()

        # Generate background model
        bkg.execute()

    # Continue only if counts cube does not exist
    if not os.path.isfile(cntcube):

        # Setup task parameters
        ctbin = ctools.ctbin()
        ctbin['inobs']    = obsname
        ctbin['outobs']   = cntcube
        ctbin['ebinalg']  = 'LOG'
        ctbin['emin']     = emin
        ctbin['emax']     = emax
        ctbin['enumbins'] = ebins
        ctbin['coordsys'] = 'CEL'
        ctbin['proj']     = 'TAN'
        ctbin['xref']     = 258.1125
        ctbin['yref']     = -39.6867
        ctbin['nxpix']    = 300
        ctbin['nypix']    = 300
        ctbin['binsz']    = 0.02
        ctbin['logfile']  = 'rx_hess_create_cntcube.log'
        ctbin.logFileOpen()

        # Generate counts cube
        ctbin.execute()

    # Continue only if counts cubes for binned analysis do not exist
    if not os.path.isfile(obsname_binned):

        # Setup task parameters
        ctbin = ctools.ctbin()
        ctbin['inobs']    = obsname
        ctbin['outobs']   = obsname_binned
        ctbin['stack']    = False
        ctbin['usepnt']   = True
        ctbin['ebinalg']  = 'LOG'
        ctbin['emin']     = emin
        ctbin['emax']     = emax
        ctbin['enumbins'] = ebins
        ctbin['coordsys'] = 'CEL'
        ctbin['proj']     = 'TAN'
        ctbin['nxpix']    = 200
        ctbin['nypix']    = 200
        ctbin['binsz']    = 0.02
        ctbin['logfile']  = 'rx_hess_create_cntcube_binned.log'
        ctbin.logFileOpen()

        # Generate counts cubes
        ctbin.execute()

    # Continue only if exposure cube does not exist
    if not os.path.isfile(expcube):

        # Setup task parameters
        ctexpcube = ctools.ctexpcube()
        ctexpcube['inobs']   = obsname
        ctexpcube['incube']  = 'NONE'
        ctexpcube['ebinalg']  = 'LOG'
        ctexpcube['emin']     =   0.1 # Full energy range
        ctexpcube['emax']     = 100.0 # Full energy range
        ctexpcube['enumbins'] =   300 # Factor ~3 oversampling of IRF
        ctexpcube['coordsys'] = 'CEL'
        ctexpcube['proj']     = 'TAN'
        ctexpcube['xref']     = 258.1125
        ctexpcube['yref']     = -39.6867
        ctexpcube['nxpix']    = 300
        ctexpcube['nypix']    = 300
        ctexpcube['binsz']    = 0.02
        ctexpcube['outcube'] = expcube
        ctexpcube['logfile'] = 'rx_hess_create_expcube.log'
        ctexpcube.logFileOpen()

        # Generate exposure cube
        ctexpcube.execute()

    # Continue only if PSF cube does not exist
    if not os.path.isfile(psfcube):

        # Setup task parameters
        ctpsfcube = ctools.ctpsfcube()
        ctpsfcube['inobs']     = obsname
        ctpsfcube['incube']    = 'NONE'
        ctpsfcube['ebinalg']   = 'LOG'
        ctpsfcube['emin']      =   0.1 # Full energy range
        ctpsfcube['emax']      = 100.0 # Full energy range
        ctpsfcube['enumbins']  =   300 # Factor ~3 oversampling of IRF
        ctpsfcube['coordsys']  = 'CEL'
        ctpsfcube['proj']      = 'TAN'
        ctpsfcube['xref']      = 258.1125
        ctpsfcube['yref']      = -39.6867
        ctpsfcube['nxpix']     = 30
        ctpsfcube['nypix']     = 30
        ctpsfcube['binsz']     = 0.2
        ctpsfcube['amax']      = 0.7 # Full H.E.S.S. PSF range
        ctpsfcube['anumbins']  = 300 # Factor ~2 oversampling of IRF
        ctpsfcube['outcube']   = psfcube
        ctpsfcube['logfile']   = 'rx_hess_create_psfcube.log'
        ctpsfcube.logFileOpen()

        # Generate PSF cube
        ctpsfcube.execute()

    # Continue only if energy dispersion cube does not exist
    if not os.path.isfile(edispcube):

        # Setup task parameters
        ctedispcube = ctools.ctedispcube()
        ctedispcube['inobs']     = obsname
        ctedispcube['incube']    = 'NONE'
        ctedispcube['ebinalg']   = 'LOG'
        ctedispcube['emin']      =   0.1 # Full energy range
        ctedispcube['emax']      = 100.0 # Full energy range
        ctedispcube['enumbins']  =   300 # Factor ~3 oversampling of IRF
        ctedispcube['coordsys']  = 'CEL'
        ctedispcube['proj']      = 'TAN'
        ctedispcube['xref']      = 258.1125
        ctedispcube['yref']      = -39.6867
        ctedispcube['nxpix']     = 30
        ctedispcube['nypix']     = 30
        ctedispcube['binsz']     = 0.2
        ctedispcube['migramax']  = 5.0
        ctedispcube['migrabins'] = 300
        ctedispcube['outcube']   = edispcube
        ctedispcube['logfile']   = 'rx_hess_create_edispcube.log'
        ctedispcube.logFileOpen()

        # Generate energy dispersion cube
        ctedispcube.execute()

    # Continue only if background cube does not exist
    if not os.path.isfile(bkgcube):

        # Setup task parameters
        ctbkgcube = ctools.ctbkgcube()
        ctbkgcube['inobs']    = obsname
        ctbkgcube['incube']   = cntcube
        ctbkgcube['inmodel']  = bkgname
        ctbkgcube['outcube']  = bkgcube
        ctbkgcube['outmodel'] = bkgname_stacked
        ctbkgcube['logfile']  = 'rx_hess_create_bkgcube.log'
        ctbkgcube.logFileOpen()

        # Generate background cube
        ctbkgcube.execute()

    # Continue only if stacked observation definition XML file does not
    # exist
    if not os.path.isfile(obsname_stacked):

        # Build stacked observation
        run = gammalib.GCTAObservation(cntcube, expcube, psfcube, bkgcube)
        run.name('RX J1713.7-3946')
        run.instrument('HESS')

        # Append to observation container
        obs = gammalib.GObservations()
        obs.append(run)

        # Save observation container
        obs.save(obsname_stacked)

    # Continue only if stacked observation definition XML file with energy
    # energy dispersion enabled does not exist
    if not os.path.isfile(obsname_stacked_edisp):

        # Build stacked observation
        run = gammalib.GCTAObservation(cntcube, expcube, psfcube, edispcube, bkgcube)
        run.name('RX J1713.7-3946')
        run.instrument('HESS')

        # Append to observation container
        obs = gammalib.GObservations()
        obs.append(run)

        # Save observation container
        obs.save(obsname_stacked_edisp)

    # Continue only if stacked model definition XML file does exist
    if os.path.isfile(bkgname_stacked):

        # Load model definition XML files
        joint   = gammalib.GModels(bkgname)
        stacked = gammalib.GModels(bkgname_stacked)

        # Get spectral component of joint file and remplace it as spectral
        # component of stacked file
        spectrum = joint[0].spectral()
        for i in range(spectrum.nodes()):
            spectrum.intensity(i, 1.0)
        spectrum.autoscale()
        stacked[0].spectral(spectrum)

        # Save stacked model
        stacked.save(bkgname_stacked)

    # Return
    return


# =========================== #
# Prepare On/Off observations #
# =========================== #
def prepare_onoff(obsname, srcmodel, spec, bkgname, emin=0.3, emax=50.0,
                 ebins=20, alpha=1.0):
    """
    Prepare On/Off observations

    Parameters
    ----------
    obsname : str
        Observation definition XML file
    srcmodel : str
        Source model
    spec : str
        Spectral model
    bkgname : str
        Background model definition XML file
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    ebins : int, optional
        Number of energy bins
    alpha : float, optional
        Template map scaling factor
    """
    # Load observations
    obs = gammalib.GObservations(obsname)

    # Set source models
    models = set_model(srcmodel, spec, bkgname, alpha=alpha)

    # Attach models
    obs.models(models)

    # Set source model
    if srcmodel == 'map':
        if alpha == 1.0:
            _srcmodel = srcmodel
        else:
            _srcmodel = 'map%.2f' % alpha
    else:
        _srcmodel = srcmodel

    # Loop over statistics
    for statistic in ['wstat','cstat']:

        # Set background model usage
        if statistic == 'wstat':
            use_model_bkg = False
        else:
            use_model_bkg = True

        # Loop over joint/stacking
        for method in ['joint','stacked']:

            # Set stacking flag
            if method == 'joint':
                stack = False
            else:
                stack = True

            # Set filenames
            outobs   = 'obs_rx_selected_onoff-%s%2.2d_%s_%s.xml'  % \
                       (method, ebins, statistic, _srcmodel)
            outmodel = 'model_rx_onoff-%s%2.2d_%s_%s.xml' % \
                       (method, ebins, statistic, _srcmodel)
            prefix   = 'onoff-%s%2.2d_%s_%s' % \
                       (method, ebins, statistic, _srcmodel)
            logfile  = 'obs_rx_selected_onoff-%s%2.2d_%s_%s.log'  % \
                       (method, ebins, statistic, _srcmodel)

            # Continue only if observation definition file does not exist
            if not os.path.isfile(outobs):

                # Generate source and background region files
                create_region_maps()

                # Setup task parameters
                csphagen = cscripts.csphagen(obs)
                csphagen['srcname']       = 'RX J1713.7-3946'
                csphagen['ebinalg']       = 'LOG'
                csphagen['emin']          = emin
                csphagen['emax']          = emax
                csphagen['enumbins']      = ebins
                csphagen['coordsys']      = 'CEL'
                csphagen['stack']         = stack
                csphagen['use_model_bkg'] = use_model_bkg
                csphagen['bkgmethod']     = 'CUSTOM'
                csphagen['bkgregmin']     = 1
                csphagen['srcregfile']    = 'rx_srcreg_map.fits'
                csphagen['bkgregfile']    = 'rx_bkgreg_map.fits'
                csphagen['etruemin']      = 0.1
                csphagen['etruemax']      = 100.0
                csphagen['etruebins']     = 300
                csphagen['outobs']        = outobs
                csphagen['outmodel']      = outmodel
                csphagen['prefix']        = prefix
                csphagen['logfile']       = logfile
                csphagen['debug']         = True
                csphagen['chatter']       = 4
                csphagen.logFileOpen()

                # Generate On/Off data
                csphagen.execute()

            # Continue only if model definition file exists
            if os.path.isfile(outmodel):

                # If we have a cstat model then set the node values and scales
                # for all nodes to unity and remove all spectral nodes above
                # 10 TeV from background model to cope with limited statistics
                # in case we have a joint fit
                if statistic == 'cstat':
                    models = gammalib.GModels(outmodel)
                    elim   = gammalib.GEnergy(10.0, 'TeV')
                    for model in models:
                        if model.type() == 'CTABackground':
                            nodes = model.spectral()
                            for i in range(nodes.nodes()-1,-1,-1):
                                if method == 'joint' and nodes.energy(i) > elim:
                                    nodes.remove(i)
                                else:
                                    nodes[i*2+1].scale(1.0)
                                    nodes[i*2+1].value(1.0)
                    models.save(outmodel)

    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Check if $HESSDATA has been set
    if 'HESSDATA' not in os.environ:
        print('HESSDATA environment variable not set. Please set HESSDATA to '
              'the root directory of the H.E.S.S. data.')
        exit(1)

    # Set parameters
    obsname          = 'obs_rx_selected.xml'
    bkgname          = 'model_rx_lookup_grad_hess.xml'
    bkgname_stacked  = 'model_rx_lookup_grad_hess_stacked.xml'
    ebins            = 40  # About 20 energy bins per decade
    spec             = 'eplaw'
    alpha            = 0.48

    # Prepare observation
    prepare(obsname, bkgname, ebins=ebins)
    prepare_onoff(obsname, 'map', spec, bkgname, ebins=ebins, alpha=alpha)

    # Unbinned analysis
    par = set_pars(edisp=True)
    analyse(obsname, bkgname, 'disk',  spec, par)
    analyse(obsname, bkgname, 'gauss', spec, par)
    analyse(obsname, bkgname, 'shell', spec, par)
    analyse(obsname, bkgname, 'map',   spec, par)
    analyse(obsname, bkgname, 'map',   spec, par, alphas=[alpha])
    analyse(obsname, bkgname, 'map',   spec, par, alphas=[0.35,0.4,0.45,0.5,0.55,0.60])

    # Create spectrum
    butterfly(obsname, bkgname, 'map', spec, par, alpha=alpha)
    spectrum(obsname,  bkgname, 'map', par, alpha=alpha)

    # Joint binned analysis
    par = set_pars(binned=True, ebins=ebins, edisp=True)
    analyse(obsname, bkgname, 'map', spec, par, alphas=[alpha])

    # Stacked analysis
    par = set_pars(stacked=True, ebins=ebins, edisp=True)
    analyse(obsname, bkgname_stacked, 'map', spec, par, alphas=[alpha])

    # On/off joint analysis
    onoff_pars_joint   = set_pars(onoff=True, ebins=ebins)
    onoff_pars_stacked = set_pars(onoff=True, stacked=True, ebins=ebins)
    analyse(obsname, '', 'map', spec, onoff_pars_joint, alphas=[alpha])
    analyse(obsname, '', 'map', spec, onoff_pars_stacked, alphas=[alpha])
    analyse(obsname, 'model_rx_onoff-joint40_cstat_map0.48.xml',   'map', spec, onoff_pars_joint, alphas=[alpha])
    analyse(obsname, 'model_rx_onoff-stacked40_cstat_map0.48.xml', 'map', spec, onoff_pars_stacked, alphas=[alpha])

    # On/off joint analysis using cstat off region weights
    par = set_pars(onoff=True, ebins=ebins, use_cstat_off_weight=True)
    analyse(obsname, '', 'map', spec, par, alphas=[alpha])
    par = set_pars(onoff=True, stacked=True, ebins=ebins, use_cstat_off_weight=True)
    analyse(obsname, '', 'map', spec, par, alphas=[alpha])

    # Generate residual map
    generate_resmap()

    # Generate sky map
    generate_skymap()
