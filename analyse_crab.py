#! /usr/bin/env python
# ==========================================================================
# Analyse Crab data
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


# ================ #
# Set source model #
# ================ #
def set_source_model(srcmodel, spec='plaw'):
    """
    Set source model

    Parameters
    ----------
    srcmodel : str
        Source model name
    spec : str, optional
        Spectral model

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
    elif spec == 'plaw145':
        spectral = gammalib.GModelSpectralPlaw(1.0e-17, -2.0,
                                               gammalib.GEnergy(1.45,'TeV'))
        spectral['Prefactor'].min(1.0e-25)
    elif spec == 'eplaw':
        spectral = gammalib.GModelSpectralExpPlaw(1.0e-17, -2.0,
                                                  gammalib.GEnergy(1.0,'TeV'),
                                                  gammalib.GEnergy(10.0,'TeV'))
        spectral['Prefactor'].min(1.0e-25)
        spectral['CutoffEnergy'].min(3.0e6)
        spectral['CutoffEnergy'].max(1.0e12)
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
    elif spec == 'fplaw':
        spectral = gammalib.GModelSpectralPlawPhotonFlux(1.0e-11, -2.0,
                                                         gammalib.GEnergy(1.0,'TeV'),
                                                         gammalib.GEnergy(100.0,'TeV'))
        spectral['PhotonFlux'].min(1.0e-25)
    elif spec == 'aharonian2006':
        spectral = gammalib.GModelSpectralExpPlaw(3.84e-17, -2.41,
                                                  gammalib.GEnergy(1.0,'TeV'),
                                                  gammalib.GEnergy(15.1,'TeV'))
        spectral['Prefactor'].fix()
        spectral['Index'].fix()
        spectral['CutoffEnergy'].fix()
    elif spec == 'nigro2019':
        spectral = gammalib.GModelSpectralLogParabola(4.47e-17, -2.39,
                                                      gammalib.GEnergy(1.0,'TeV'),
                                                      -0.16)
        spectral['Prefactor'].fix()
        spectral['Index'].fix()
        spectral['Curvature'].fix()

    # Set spatial component
    dir = gammalib.GSkyDir()
    dir.radec_deg(83.6331, 22.0145)
    if srcmodel == 'ptsrc':
        spatial = gammalib.GModelSpatialPointSource(dir)
    elif srcmodel == 'gauss':
        spatial = gammalib.GModelSpatialRadialGauss(dir, 0.1)
        spatial['Sigma'].min(0.0001)
        spatial['Sigma'].free()
    elif srcmodel == 'egauss':
        spatial = gammalib.GModelSpatialEllipticalGauss(dir, 0.02, 0.01, 0.0)
        spatial['MinorRadius'].min(0.0001)
        spatial['MinorRadius'].free()
        spatial['MajorRadius'].min(0.0001)
        spatial['MajorRadius'].free()
        spatial['PA'].free()
    spatial['RA'].free()
    spatial['DEC'].free()

    # Set source model
    source = gammalib.GModelSky(spatial, spectral)
    source.name('Crab')
    source.tscalc(True)

    # Return source model
    return source


# ========= #
# Set model #
# ========= #
def set_model(srcmodel, spec, bkgname):
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

    Returns
    -------
    models : `~gammalib.GModels()`
        Models
    """
    # Initialise model container
    models = gammalib.GModels()

    # Set source model
    if srcmodel != '':
        source = set_source_model(srcmodel, spec=spec)
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
def set_analysis(srcmodel, bkgname, spec, pars, select=None):
    """
    Set analysis string

    Parameters
    ----------
    srcmodels : str
        Source model
    bkgname : str
        Background name
    spec : str
        Spectral model
    pars : dict
        Analysis parameters
    select : list of int, optional
        Indices for observation selection
    """
    # Set analysis
    analysis = '%s_%s_lookup_grad_hess' % (srcmodel, spec)
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
def set_observation(obsname, srcmodel, bkgname, pars):
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
    """
    # Set observation name
    if pars['onoff']:
        if bkgname == '':
            statistic = 'wstat'
        else:
            statistic = 'cstat'
        if pars['stacked']:
            _obsname = add_attribute(obsname, '_onoff-stacked%2.2d_%s_%s' %
                                     (pars['ebins'], statistic, srcmodel))
        else:
            _obsname = add_attribute(obsname, '_onoff-joint%2.2d_%s_%s' %
                                     (pars['ebins'], statistic, srcmodel))
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
def set_pars(stacked=False, binned=False, onoff=False, ebins=20, edisp=False):
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
    """
    # Set analysis parameters
    pars = {'stacked': stacked,
            'binned' : binned,
            'onoff':   onoff,
            'ebins':   ebins,
            'edisp':   edisp}

    # Return
    return pars


# ================ #
# Fit observations #
# ================ #
def fit(obs, srcmodel, spec, bkgname, analysis, pars):
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
    """
    # Set file names
    outmodel = 'crab_results_%s.xml' % analysis
    logfile  = 'crab_results_%s.log' % analysis

    # Continue only if result file does not exist
    if not os.path.isfile(outmodel):

        # Set models
        models = set_model(srcmodel, spec, bkgname)

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
    inmodel  = 'crab_results_%s.xml' % analysis
    outmodel = 'crab_npred_%s.xml' % analysis
    logfile  = 'crab_npred_%s.log' % analysis

    # Continue only if result file does not exist
    if not os.path.isfile(outmodel) and os.path.isfile(inmodel):

        # Set models
        models = gammalib.GModels(inmodel)
        model  = models['Crab']
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
def analyse(obsname, bkgname, srcmodel, spec, pars, select=None):
    """
    Analyse observations

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
    select : list of int, optional
        Indices for observation selection
    """
    # Load observations
    _obsname = set_observation(obsname, srcmodel, bkgname, pars)
    obs = gammalib.GObservations(_obsname)

    # Optionally select observations
    if select != None:
        obs_select = gammalib.GObservations()
        for i in select:
            obs_select.append(obs[i])
        obs = obs_select

    # Set analysis
    analysis = set_analysis(srcmodel, bkgname, spec, pars, select=select)

    # Fit model
    fit(obs, srcmodel, spec, bkgname, analysis, pars)

    # Determine Npred
    npred(obs, analysis, pars)

    # Return
    return


# ================== #
# Generate butterfly #
# ================== #
def butterfly(obsname, bkgname, srcmodel, spec, pars):
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
    """
    # Set observation name
    _obsname = set_observation(obsname, srcmodel, bkgname, pars)

    # Set analysis
    analysis = set_analysis(srcmodel, bkgname, spec, pars)

    # Set file names
    inmodel = 'crab_results_%s.xml' % analysis
    outfile = 'crab_butterfly_%s.txt' % analysis
    logfile = 'crab_butterfly_%s.log' % analysis

    # Continue only if result file does not exist
    if not os.path.isfile(outfile) and os.path.isfile(inmodel):

        # Generate butterfly
        butterfly = ctools.ctbutterfly()
        butterfly['inobs']    = _obsname
        butterfly['inmodel']  = inmodel
        butterfly['srcname']  = 'Crab'
        butterfly['outfile']  = outfile
        butterfly['edisp']    = pars['edisp']
        butterfly['emin']     = 0.3
        butterfly['emax']     = 40.0
        butterfly['logfile']  = logfile
        butterfly['debug']    = True
        butterfly.logFileOpen()
        butterfly.execute()

    # Return
    return


# ================= #
# Generate spectrum #
# ================= #
def spectrum(obsname, bkgname, srcmodel, pars, emin=0.7, emax=25.0):
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
    """
    # Set observation name
    _obsname = set_observation(obsname, srcmodel, bkgname, pars)

    # Set analysis
    analysis = set_analysis(srcmodel, bkgname, 'plaw', pars)

    # Set file names
    outfile = 'crab_spectrum_%s.fits' % (analysis)
    logfile = 'crab_spectrum_%s.log'  % (analysis)

    # Continue only if result file does not exist
    if not os.path.isfile(outfile):

        # Load observations
        obs = gammalib.GObservations(_obsname)

        # Set background model for stacked analysis
        if pars['stacked']:
            bname  = add_attribute(bkgname, '_stacked')
            models = set_model(srcmodel, 'plaw', bname)
            for model in models:
                if model.type() == 'CTACubeBackground':
                    model.spectral()['Index'].fix()

        # ... otherwise create H.E.S.S. background model with power law
        # spectral component
        else:

            # Setup plaw background model
            bkg = cscripts.csbkgmodel()
            bkg['inobs']      = '$HESSDATA/obs/obs_crab.xml'
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
            source = set_source_model(srcmodel, spec='plaw')
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
        spec['srcname']  = 'Crab'
        spec['edisp']    = pars['edisp']
        spec['outfile']  = outfile
        spec['method']   = 'SLICE'
        spec['ebinalg']  = 'LOG'
        spec['emin']     = emin
        spec['emax']     = emax
        spec['enumbins'] = 10
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
def generate_resmap(emin=0.670, emax=30.0, ebins=40, edisp=True):
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
    inobs    = 'obs_crab_selected.xml'
    inmodel  = 'crab_results_gauss_eplaw_lookup_grad_hess_edisp.xml'
    outmodel = 'crab_background_gauss_eplaw_lookup_grad_hess.xml'
    outmap   = 'crab_resmap.fits'
    logfile  = 'crab_resmap.log'

    # Continue only if inmodel exists
    if os.path.isfile(inobs) and os.path.isfile(inmodel):

        # Create model definition XML file without source
        models = gammalib.GModels(inmodel)
        models.remove('Crab')
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
        resmap['xref']      = 83.63
        resmap['yref']      = 22.01
        resmap['nxpix']     = 350
        resmap['nypix']     = 250
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
def generate_skymap(emin=0.670, emax=30.0):
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
    inobs    = 'obs_crab_selected.xml'
    outmap   = 'crab_skymap.fits'
    logfile  = 'crab_skymap.log'

    # Continue only if inmodel exists
    if os.path.isfile(inobs) and not os.path.isfile(outmap):

        # Setup task parameters
        skymap = ctools.ctskymap()
        skymap['inobs']       = inobs
        skymap['emin']        = emin
        skymap['emax']        = emax
        skymap['nxpix']       = 200
        skymap['nypix']       = 200
        skymap['binsz']       = 0.02
        skymap['coordsys']    = 'CEL'
        skymap['proj']        = 'CAR'
        skymap['xref']        = 83.63
        skymap['yref']        = 22.01
        skymap['bkgsubtract'] = 'RING'
        skymap['roiradius']   = 0.1
        skymap['inradius']    = 0.6
        skymap['outradius']   = 0.8
        skymap['iterations']  = 3
        skymap['threshold']   = 5.0
        skymap['inexclusion'] = 'NONE'
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
def prepare(obsname, bkgname, rad=2.0, emin=0.67, emax=30.0, ebins=20):
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
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    ebins : int, optional
        Number of energy bins
    """
    # Set filenames
    cntcube               = 'crab_stacked%2.2d_cntcube.fits'   % ebins
    expcube               = 'crab_stacked%2.2d_expcube.fits'   % ebins
    psfcube               = 'crab_stacked%2.2d_psfcube.fits'   % ebins
    edispcube             = 'crab_stacked%2.2d_edispcube.fits' % ebins
    bkgcube               = 'crab_stacked%2.2d_bkgcube.fits'   % ebins
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
        select['inobs']    = '$HESSDATA/obs/obs_crab.xml'
        select['outobs']   = obsname
        select['ra']       = 'UNDEF'
        select['dec']      = 'UNDEF'
        select['rad']      = rad
        select['tmin']     = 'UNDEF'
        select['tmax']     = 'UNDEF'
        select['emin']     = emin
        select['emax']     = emax
        select['usethres'] = 'DEFAULT'
        select['logfile']  = 'crab_hess_select_events.log'
        select.logFileOpen()

        # Select events
        select.execute()

    # Continue only if background model does not exist
    if not os.path.isfile(bkgname):

        # Setup task parameters
        bkg = cscripts.csbkgmodel()
        bkg['inobs']      = '$HESSDATA/obs/obs_crab.xml'
        bkg['outmodel']   = bkgname
        bkg['instrument'] = 'HESS'
        bkg['spatial']    = 'LOOKUP'
        bkg['slufile']    = 'off_lookup.fits'
        bkg['gradient']   = True
        bkg['spectral']   = 'NODES'
        bkg['ebinalg']    = 'LOG'
        bkg['emin']       = emin
        bkg['emax']       = emax
        bkg['enumbins']   = 8
        bkg['runwise']    = True
        bkg['rad']        = rad
        bkg['logfile']    = 'crab_hess_create_background.log'
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
        ctbin['xref']     = 83.63
        ctbin['yref']     = 22.01
        ctbin['nxpix']    = 350
        ctbin['nypix']    = 250
        ctbin['binsz']    = 0.02
        ctbin['logfile']  = 'crab_hess_create_cntcube.log'
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
        ctbin['logfile']  = 'crab_hess_create_cntcube_binned.log'
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
        ctexpcube['xref']     = 83.63
        ctexpcube['yref']     = 22.01
        ctexpcube['nxpix']    = 350
        ctexpcube['nypix']    = 250
        ctexpcube['binsz']    = 0.02
        ctexpcube['outcube'] = expcube
        ctexpcube['logfile'] = 'crab_hess_create_expcube.log'
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
        ctpsfcube['xref']      = 83.63
        ctpsfcube['yref']      = 22.01
        ctpsfcube['nxpix']     = 35
        ctpsfcube['nypix']     = 25
        ctpsfcube['binsz']     = 0.2
        ctpsfcube['amax']      = 0.7 # Full H.E.S.S. PSF range
        ctpsfcube['anumbins']  = 300 # Factor ~2 oversampling of IRF
        ctpsfcube['outcube']   = psfcube
        ctpsfcube['logfile']   = 'crab_hess_create_psfcube.log'
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
        ctedispcube['xref']      = 83.63
        ctedispcube['yref']      = 22.01
        ctedispcube['nxpix']     = 35
        ctedispcube['nypix']     = 25
        ctedispcube['binsz']     = 0.2
        ctedispcube['migramax']  = 5.0
        ctedispcube['migrabins'] = 300
        ctedispcube['outcube']   = edispcube
        ctedispcube['logfile']   = 'crab_hess_create_edispcube.log'
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
        ctbkgcube['logfile']  = 'crab_hess_create_bkgcube.log'
        ctbkgcube.logFileOpen()

        # Generate background cube
        ctbkgcube.execute()

    # Continue only if stacked observation definition XML file does not
    # exist
    if not os.path.isfile(obsname_stacked):

        # Build stacked observation
        run = gammalib.GCTAObservation(cntcube, expcube, psfcube, bkgcube)
        run.name('Crab')
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
        run.name('Crab')
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
def prepare_onoff(obsname, srcmodel, spec, bkgname, emin=0.67, emax=30.0, ebins=20):
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
    """
    # Load observations
    obs = gammalib.GObservations(obsname)

    # Set source models
    models = set_model(srcmodel, spec, bkgname)

    # Attach models
    obs.models(models)

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
            outobs   = 'obs_crab_selected_onoff-%s%2.2d_%s_%s.xml'  % \
                       (method, ebins, statistic, srcmodel)
            outmodel = 'model_crab_onoff-%s%2.2d_%s_%s.xml' % \
                       (method, ebins, statistic, srcmodel)
            prefix   = 'onoff-%s%2.2d_%s_%s' % \
                       (method, ebins, statistic, srcmodel)
            logfile  = 'obs_crab_selected_onoff-%s%2.2d_%s_%s.log'  % \
                       (method, ebins, statistic, srcmodel)

            # Continue only if observation definition file does not exist
            if not os.path.isfile(outobs):

                # Setup task parameters
                csphagen = cscripts.csphagen(obs)
                csphagen['srcname']       = 'Crab'
                csphagen['ebinalg']       = 'LOG'
                csphagen['emin']          = emin
                csphagen['emax']          = emax
                csphagen['enumbins']      = ebins
                csphagen['coordsys']      = 'CEL'
                csphagen['ra']            = 83.6330
                csphagen['dec']           = 22.0145
                csphagen['rad']           = 0.2
                csphagen['stack']         = stack
                csphagen['use_model_bkg'] = use_model_bkg
                csphagen['bkgmethod']     = 'REFLECTED'
                csphagen['etruemin']      = 0.1
                csphagen['etruemax']      = 100.0
                csphagen['etruebins']     = 300
                csphagen['outobs']        = outobs
                csphagen['outmodel']      = outmodel
                csphagen['prefix']        = prefix
                csphagen['logfile']       = logfile
                csphagen.logFileOpen()

                # Generate On/Off data
                csphagen.execute()

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
    obsname          = 'obs_crab_selected.xml'
    bkgname          = 'model_crab_lookup_grad_hess.xml'
    bkgname_stacked  = 'model_crab_lookup_grad_hess_stacked.xml'
    ebins            = 40  # About 20 energy bins per decade

    # Prepare observation
    prepare(obsname, bkgname, ebins=ebins)
    prepare_onoff(obsname, 'ptsrc', 'eplaw', bkgname, ebins=ebins)

    # Unbinned analysis
    par = set_pars(edisp=True)
    analyse(obsname,   bkgname, 'ptsrc', 'plaw', par)
    analyse(obsname,   bkgname, 'ptsrc', 'eplaw', par)
    analyse(obsname,   bkgname, 'ptsrc', 'logparabola', par)
    analyse(obsname,   bkgname, 'ptsrc', 'aharonian2006', par)
    analyse(obsname,   bkgname, 'ptsrc', 'nigro2019', par)
    butterfly(obsname, bkgname, 'ptsrc', 'eplaw', par)
    spectrum(obsname,  bkgname, 'ptsrc', par)

    # Spatial fitting of a Gaussian (unbinned)
    par = set_pars(edisp=True)
    analyse(obsname, bkgname, 'gauss', 'eplaw', par)

    # Joint binned analysis
    par = set_pars(binned=True, ebins=ebins, edisp=True)
    analyse(obsname, bkgname, 'ptsrc', 'eplaw', par)

    # Stacked analysis
    par = set_pars(stacked=True, ebins=ebins, edisp=True)
    analyse(obsname, bkgname_stacked, 'ptsrc', 'eplaw', par)

    # On/off analysis
    pars_joint   = set_pars(onoff=True, ebins=ebins)
    pars_stacked = set_pars(onoff=True, stacked=True, ebins=ebins)
    analyse(obsname, '',                                           'ptsrc', 'eplaw', pars_joint)
    analyse(obsname, 'model_crab_onoff-joint40_cstat_ptsrc.xml',   'ptsrc', 'eplaw', pars_joint)
    analyse(obsname, '',                                           'ptsrc', 'eplaw', pars_stacked)
    analyse(obsname, 'model_crab_onoff-stacked40_cstat_ptsrc.xml', 'ptsrc', 'eplaw', pars_stacked)

    # Generate residual map
    generate_resmap()

    # Generate sky map
    generate_skymap()

