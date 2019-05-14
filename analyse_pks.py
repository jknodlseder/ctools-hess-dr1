#! /usr/bin/env python
# ==========================================================================
# Analyse PKS 2155-304 data
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
        spectral = gammalib.GModelSpectralPlaw(7.0e-17, -3.0,
                                               gammalib.GEnergy(1.0,'TeV'))
        spectral['Prefactor'].min(1.0e-25)
    elif spec == 'eplaw':
        spectral = gammalib.GModelSpectralExpPlaw(7.0e-17, -3.0,
                                                  gammalib.GEnergy(1.0,'TeV'),
                                                  gammalib.GEnergy(3.0,'TeV'))
        spectral['Prefactor'].min(1.0e-25)
        spectral['CutoffEnergy'].min(1.0e5)
        spectral['CutoffEnergy'].max(1.0e8)
    elif spec == 'logparabola':
        spectral = gammalib.GModelSpectralLogParabola(7.0e-17, -3.0,
                                                      gammalib.GEnergy(1.0,'TeV'),
                                                      -0.8)
        spectral['Prefactor'].min(1.0e-25)
    elif spec == 'fplaw200':
        spectral = gammalib.GModelSpectralPlawPhotonFlux(1.0e-10, -2.0,
                                                         gammalib.GEnergy(0.2,'TeV'),
                                                         gammalib.GEnergy(100.0,'TeV'))
        spectral['PhotonFlux'].min(1.0e-25)
    elif spec == 'fplaw300':
        spectral = gammalib.GModelSpectralPlawPhotonFlux(1.0e-10, -2.0,
                                                         gammalib.GEnergy(0.3,'TeV'),
                                                         gammalib.GEnergy(100.0,'TeV'))
        spectral['PhotonFlux'].min(1.0e-25)
    elif spec == 'fplaw700':
        spectral = gammalib.GModelSpectralPlawPhotonFlux(1.0e-10, -3.3,
                                                         gammalib.GEnergy(0.7,'TeV'),
                                                         gammalib.GEnergy(100.0,'TeV'))
        spectral['PhotonFlux'].min(1.0e-25)
        spectral['Index'].min(-5.0)
        spectral['Index'].max(-1.5)

    # Set spatial component
    dir = gammalib.GSkyDir()
    dir.radec_deg(329.71694, -30.22559)
    if srcmodel == 'ptsrc':
        spatial = gammalib.GModelSpatialPointSource(dir)
    elif srcmodel == 'gauss':
        spatial = gammalib.GModelSpatialRadialGauss(dir, 0.1)
        spatial['Sigma'].min(0.0001)
        spatial['Sigma'].free()
    spatial['RA'].free()
    spatial['DEC'].free()

    # Set source model
    source = gammalib.GModelSky(spatial, spectral)
    source.name('PKS 2155-304')
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
def set_pars(stacked=False, binned=False, onoff=False, ebins=0, edisp=False):
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
def fit(obs, srcmodel, spec, bkgname, analysis, pars, fitname='pks_results'):
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
    fitname : str, optional
        Fit result prefix
    """
    # Set file names
    outmodel = '%s_%s.xml' % (fitname, analysis)
    logfile  = '%s_%s.log' % (fitname, analysis)

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
        like['logfile']  = logfile
        like['max_iter'] = 200
        like['debug']    = True
        like.logFileOpen()

        # Use wstat if no background model is provided
        if bkgname == '':
            like['statistic'] = 'WSTAT'

        # Execute ctlike
        like.execute()

    # Compute energy flux and add it to log file
    if os.path.isfile(outmodel):

        # Load result file and compute energy flux
        emin   = gammalib.GEnergy(0.3, 'TeV')
        emax   = gammalib.GEnergy(3.0, 'TeV')
        models = gammalib.GModels(outmodel)
        eflux  = models['PKS 2155-304'].spectral().eflux(emin, emax)

        # Append result to output file
        with open(logfile, "a") as myfile:
            myfile.write('\n')
            myfile.write('Energy flux (0.3-3 TeV): %e' % eflux)

    # Return
    return


# ========================= #
# Determine Npred of source #
# ========================= #
def npred(obs, analysis, pars, fitname='pks_results', npredname='pks_npred'):
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
    fitname : str, optional
        Fit result prefix
    npredname : str, optional
        Npred result prefix
    """
    # Set file names
    inmodel  = '%s_%s.xml' % (fitname,  analysis)
    outmodel = '%s_%s.xml' % (npredname, analysis)
    logfile  = '%s_%s.log' % (npredname, analysis)

    # Continue only if result file does not exist
    if not os.path.isfile(outmodel) and os.path.isfile(inmodel):

        # Set models
        models = gammalib.GModels(inmodel)
        model  = models['PKS 2155-304']
        for par in model:
            par.fix()
        model.tscalc(False) # Otherwise leads to an error
        npred_models = gammalib.GModels()
        npred_models.append(model)

        # Attach models
        obs.models(npred_models)

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
def analyse(obsname, bkgname, srcmodel, spec, pars, select=None,
            fitname='pks_results', npredname='pks_npred'):
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
    fitname : str, optional
        Fit result prefix
    npredname : str, optional
        Npred result prefix
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
    fit(obs, srcmodel, spec, bkgname, analysis, pars, fitname=fitname)

    # Determine Npred
    npred(obs, analysis, pars, fitname=fitname, npredname=npredname)

    # Return
    return


# =================== #
# Generate lightcurve #
# =================== #
def lightcurve(obsname, bkgname, srcmodel, spec, method, pars,
               emin=0.7, emax=10.0):
    """
    Generate lightcurve

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
    method : str
        Light curve method
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
    analysis = set_analysis(srcmodel, bkgname, spec, pars)

    # Set file names
    if method == 'ONOFF':
        outfile = 'pks_lightcrv_%s_%s_%dbins.fits' % (analysis, method, pars['ebins'])
        logfile = 'pks_lightcrv_%s_%s_%dbins.log'  % (analysis, method, pars['ebins'])
    else:
        outfile = 'pks_lightcrv_%s_%s.fits' % (analysis, method)
        logfile = 'pks_lightcrv_%s_%s.log'  % (analysis, method)

    # Continue only if result file does not exist
    if not os.path.isfile(outfile):

        # Create GTIs for light curve
        create_lightcurve_gti(obsname_t700, 'lightcurve_gti.fits')

        # Load observations
        obs = gammalib.GObservations(_obsname)

        # Set models
        models = set_model(srcmodel, spec, bkgname)

        # Set and fix source position
        models['PKS 2155-304']['RA'].value(329.71694)   # true
        models['PKS 2155-304']['DEC'].value(-30.225588) # true
        models['PKS 2155-304']['RA'].fix()
        models['PKS 2155-304']['DEC'].fix()

        # Set and fix spectral index
        models['PKS 2155-304']['Index'].value(-3.4)
        models['PKS 2155-304']['Index'].fix()

        # Attach models
        obs.models(models)

        # Generate light curve
        lightcurve = cscripts.cslightcrv(obs)
        lightcurve['srcname']   = 'PKS 2155-304'
        lightcurve['edisp']     = pars['edisp']
        lightcurve['outfile']   = outfile
        lightcurve['tbinalg']   = 'FILE'
        lightcurve['tbinfile']  = 'lightcurve_gti.fits'
        lightcurve['method']    = method
        lightcurve['emin']      = emin
        lightcurve['emax']      = emax
        lightcurve['coordsys']  = 'CEL'
        lightcurve['xref']      = 329.71694
        lightcurve['yref']      = -30.22559
        lightcurve['rad']       =   0.2
        lightcurve['enumbins']  = pars['ebins']
        lightcurve['etruebins'] = 100
        lightcurve['logfile']   = logfile
        lightcurve['clobber']   = 4
        lightcurve['debug']     = True
        if method == 'ONOFF':
            lightcurve['statistic']     = 'WSTAT'
            lightcurve['use_model_bkg'] = False
        lightcurve.logFileOpen()
        lightcurve.execute()

    # Return
    return


# =================== #
# Generate lightcurve #
# =================== #
def create_lightcurve_gti(obsname, gtiname):
    """
    Generate lightcurve for GTIs

    Parameters
    ----------
    obsname : str
        Observation definition XML file
    gtiname : str
        GTI output filename
    """
    # Load observations
    obs = gammalib.GObservations(obsname)

    # Initialize GTIs
    gti = gammalib.GGti()

    # Loop over observations
    for run in obs:

        # Get start time and ontime of GTIs
        tstart = run.gti().tstart().copy()
        ontime = run.gti().ontime()

        # Set time bins
        tbin = ontime / 14.0

        # Loop over time bins
        for i in range(14):

            # Compute stop time
            tstop = tstart + tbin

            # Append GTI
            gti.append(tstart, tstop)

            # Update start time
            tstart = tstop.copy()

    # Save GTIs
    gti.save(gtiname, True)

    # Return
    return


# ================== #
# Generate butterfly #
# ================== #
def butterfly(obsname, bkgname, srcmodel, spec, pars, emin=0.3, emax=10.0,
              fitname='pks_results', buttername='pks_butterfly'):
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
    emin : float, optional
        Minimum energy (TeV)
    emax : float, optional
        Maximum energy (TeV)
    fitname : str, optional
        Fit result prefix
    buttername : str, optional
        Butterfly diagram result prefix
    """
    # Set observation name
    _obsname = set_observation(obsname, srcmodel, bkgname, pars)

    # Set analysis
    analysis = set_analysis(srcmodel, bkgname, spec, pars)

    # Set file names
    inmodel = '%s_%s.xml'   % (fitname, analysis)
    outfile = '%s_%s.txt' % (buttername, analysis)
    logfile = '%s_%s.log' % (buttername, analysis)

    # Continue only if result file does not exist
    if not os.path.isfile(outfile) and os.path.isfile(inmodel):

        # Generate butterfly
        butterfly = ctools.ctbutterfly()
        butterfly['inobs']    = _obsname
        butterfly['inmodel']  = inmodel
        butterfly['srcname']  = 'PKS 2155-304'
        butterfly['outfile']  = outfile
        butterfly['edisp']    = pars['edisp']
        butterfly['emin']     = emin
        butterfly['emax']     = emax
        butterfly['logfile']  = logfile
        butterfly['debug']    = True
        butterfly.logFileOpen()
        butterfly.execute()

    # Return
    return


# ================= #
# Generate spectrum #
# ================= #
def spectrum(obsname, bkgname, srcmodel, pars, emin=0.7, emax=10.0, ebins=10,
             fitname='pks_results', specname='pks_spectrum'):
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
    ebins : int, optional
        Number of energy bins
    fitname : str, optional
        Fit result prefix
    specname : str, optional
        Spectrum result prefix
    """
    # Set observation name
    _obsname = set_observation(obsname, srcmodel, bkgname, pars)

    # Set analysis
    analysis = set_analysis(srcmodel, bkgname, 'plaw', pars)

    # Set file names
    outfile = '%s_%s_%dbins.fits' % (specname, analysis, ebins)
    logfile = '%s_%s_%dbins.log'  % (specname, analysis, ebins)

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
            bkg['inobs']      = obsname
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
        spec['srcname']  = 'PKS 2155-304'
        spec['edisp']    = pars['edisp']
        spec['outfile']  = outfile
        spec['method']   = 'SLICE'
        spec['ebinalg']  = 'LOG'
        spec['emin']     = emin
        spec['emax']     = emax
        spec['enumbins'] = ebins
        spec['logfile']  = logfile
        spec['debug']    = True
        spec['chatter']  = 4
        spec.logFileOpen()
        spec.execute()

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
def prepare(obsname, bkgname, rad=2.0, emin=0.70, emax=10.0, ebins=20,
            tmin=None, tmax=None):
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
    tmin : str, optional
        Start time
    tmax : str, optional
        Stop time
    """
    # Set filenames
    thresmin              = emin*1000.0
    cntcube               = 'pks_t%3.3d_stacked%2.2d_cntcube.fits'   % (thresmin, ebins)
    expcube               = 'pks_t%3.3d_stacked%2.2d_expcube.fits'   % (thresmin, ebins)
    psfcube               = 'pks_t%3.3d_stacked%2.2d_psfcube.fits'   % (thresmin, ebins)
    edispcube             = 'pks_t%3.3d_stacked%2.2d_edispcube.fits' % (thresmin, ebins)
    bkgcube               = 'pks_t%3.3d_stacked%2.2d_bkgcube.fits'   % (thresmin, ebins)
    obsname_binned        = add_attribute(obsname, '_binned%2.2d' % ebins)
    bkgname_stacked       = add_attribute(bkgname, '_stacked')
    obsname_stacked       = add_attribute(obsname, '_stacked%2.2d' % ebins)
    obsname_stacked_edisp = add_attribute(obsname_stacked, '_edisp')

    # Generate background lookup
    generate_background_lookup()

    # Continue only if selected events do not exist
    if not os.path.isfile(obsname):

        # Select observations if tmin/tmax is specified
        if tmin != None and tmax != None:

            # Select observations
            select = cscripts.csobsselect()
            select['inobs']     = '$HESSDATA/obs/obs_pks_flare.xml'
            select['pntselect'] = 'CIRCLE'
            select['coordsys']  = 'GAL'
            select['glon']      =   0.0
            select['glat']      =   0.0
            select['rad']       = 180.0
            select['tmin']      = tmin
            select['tmax']      = tmax
            select['debug']     = True
            select.run()

            # Setup task parameters
            select = ctools.ctselect(select.obs())
            select['outobs']   = obsname
            select['prefix']   = 'selected_t%3.3d_' % (thresmin)
            select['ra']       = 'UNDEF'
            select['dec']      = 'UNDEF'
            select['rad']      = rad
            select['tmin']     = tmin
            select['tmax']     = tmax
            select['emin']     = emin
            select['emax']     = emax
            select['usethres'] = 'DEFAULT'
            select['logfile']  = 'pks_t%3.3d_hess_select_events.log' % (thresmin)
            select.logFileOpen()

            # Select events
            select.execute()

        else:

            # Setup task parameters
            select = ctools.ctselect()
            select['inobs']    = '$HESSDATA/obs/obs_pks_flare.xml'
            select['outobs']   = obsname
            select['prefix']   = 'selected_t%3.3d_' % (thresmin)
            select['ra']       = 'UNDEF'
            select['dec']      = 'UNDEF'
            select['rad']      = rad
            select['tmin']     = 'UNDEF'
            select['tmax']     = 'UNDEF'
            select['emin']     = emin
            select['emax']     = emax
            select['usethres'] = 'DEFAULT'
            select['logfile']  = 'pks_t%3.3d_hess_select_events.log' % (thresmin)
            select.logFileOpen()

            # Select events
            select.execute()

    # Continue only if background model does not exist
    if os.path.isfile(obsname) and not os.path.isfile(bkgname):

        # Setup task parameters
        bkg = cscripts.csbkgmodel()
        bkg['inobs']      = obsname
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
        bkg['logfile']    = 'pks_t%3.3d_hess_create_background.log' % (thresmin)
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
        ctbin['xref']     = 329.71694
        ctbin['yref']     = -30.22559
        ctbin['nxpix']    = 250
        ctbin['nypix']    = 250
        ctbin['binsz']    = 0.02
        ctbin['logfile']  = 'pks_t%3.3d_hess_create_cntcube.log' % (thresmin)
        ctbin.logFileOpen()

        # Generate counts cube
        ctbin.execute()

    # Continue only if counts cubes for binned analysis do not exist
    if not os.path.isfile(obsname_binned):

        # Setup task parameters
        ctbin = ctools.ctbin()
        ctbin['inobs']    = obsname
        ctbin['prefix']   = 'cntcube_t%3.3d_' % (thresmin)
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
        ctbin['logfile']  = 'pks_t%3.3d_hess_create_cntcube_binned.log' % (thresmin)
        ctbin.logFileOpen()

        # Generate counts cubes
        ctbin.execute()

    # Continue only if exposure cube does not exist
    if not os.path.isfile(expcube):

        # Setup task parameters
        ctexpcube = ctools.ctexpcube()
        ctexpcube['inobs']    = obsname
        ctexpcube['incube']   = 'NONE'
        ctexpcube['ebinalg']  = 'LOG'
        ctexpcube['emin']     =   0.1 # Full energy range
        ctexpcube['emax']     = 100.0 # Full energy range
        ctexpcube['enumbins'] =   300 # Factor ~3 oversampling of IRF
        ctexpcube['coordsys'] = 'CEL'
        ctexpcube['proj']     = 'TAN'
        ctexpcube['xref']     = 329.71694
        ctexpcube['yref']     = -30.22559
        ctexpcube['nxpix']    = 250
        ctexpcube['nypix']    = 250
        ctexpcube['binsz']    = 0.02
        ctexpcube['outcube']  = expcube
        ctexpcube['logfile']  = 'pks_hess_create_expcube.log'
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
        ctpsfcube['xref']      = 329.71694
        ctpsfcube['yref']      = -30.22559
        ctpsfcube['nxpix']     = 25
        ctpsfcube['nypix']     = 25
        ctpsfcube['binsz']     = 0.2
        ctpsfcube['amax']      = 0.7 # Full H.E.S.S. PSF range
        ctpsfcube['anumbins']  = 300 # Factor ~2 oversampling of IRF
        ctpsfcube['outcube']   = psfcube
        ctpsfcube['logfile']   = 'pks_hess_create_psfcube.log'
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
        ctedispcube['xref']      = 329.71694
        ctedispcube['yref']      = -30.22559
        ctedispcube['nxpix']     = 25
        ctedispcube['nypix']     = 25
        ctedispcube['binsz']     = 0.2
        ctedispcube['migramax']  = 5.0
        ctedispcube['migrabins'] = 300
        ctedispcube['outcube']   = edispcube
        ctedispcube['logfile']   = 'pks_hess_create_edispcube.log'
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
        ctbkgcube['logfile']  = 'pks_t%3.3d_hess_create_bkgcube.log' % (thresmin)
        ctbkgcube.logFileOpen()

        # Generate background cube
        ctbkgcube.execute()

    # Continue only if stacked observation definition XML file does not
    # exist
    if not os.path.isfile(obsname_stacked):

        # Build stacked observation
        run = gammalib.GCTAObservation(cntcube, expcube, psfcube, bkgcube)
        run.name('PKS 2155-304')
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
        run.name('PKS 2155-304')
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
def prepare_onoff(obsname, srcmodel, spec, bkgname, emin=0.70, emax=10.0, ebins=20):
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
    # Set Txxx value
    thresmin = emin*1000.0

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
            outobs   = 'obs_pks_t%3.3d_selected_onoff-%s%2.2d_%s_%s.xml'  % \
                       (thresmin, method, ebins, statistic, srcmodel)
            outmodel = 'model_pks_t%3.3d_onoff-%s%2.2d_%s_%s.xml' % \
                       (thresmin, method, ebins, statistic, srcmodel)
            prefix   = 't%3.3d_onoff-%s%2.2d_%s_%s' % \
                       (thresmin, method, ebins, statistic, srcmodel)
            logfile  = 'obs_pks_t%3.3d_selected_onoff-%s%2.2d_%s_%s.log'  % \
                       (thresmin, method, ebins, statistic, srcmodel)

            # Continue only if observation definition file does not exist
            if not os.path.isfile(outobs):

                # Setup task parameters
                csphagen = cscripts.csphagen(obs)
                csphagen['srcname']       = 'PKS 2155-304'
                csphagen['ebinalg']       = 'LOG'
                csphagen['emin']          = emin
                csphagen['emax']          = emax
                csphagen['enumbins']      = ebins
                csphagen['coordsys']      = 'CEL'
                csphagen['ra']            = 329.71694
                csphagen['dec']           = -30.22559
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
    obsname_t200         = 'obs_pks_t200_selected.xml'
    obsname_t300         = 'obs_pks_t300_selected.xml'
    obsname_t700         = 'obs_pks_t700_selected.xml'
    bkgname_t200         = 'model_pks_t200_lookup_grad_hess.xml'
    bkgname_t300         = 'model_pks_t300_lookup_grad_hess.xml'
    bkgname_t700         = 'model_pks_t700_lookup_grad_hess.xml'
    bkgname_t200_stacked = 'model_pks_t200_lookup_grad_hess_stacked.xml'
    bkgname_t300_stacked = 'model_pks_t300_lookup_grad_hess_stacked.xml'
    bkgname_t700_stacked = 'model_pks_t700_lookup_grad_hess_stacked.xml'
    ebins                = 40  # About 20 energy bins per decade

    # Prepare observations
    prepare(obsname_t200, bkgname_t200, ebins=ebins, emin=0.2,
            tmin='MJD 53945.934618', tmax='MJD 53946.107844')
    prepare(obsname_t300, bkgname_t300, ebins=ebins, emin=0.3,
            tmin='MJD 53945.913536', tmax='MJD 53946.129121')
    prepare(obsname_t700, bkgname_t700, ebins=ebins, emin=0.7)
    prepare_onoff(obsname_t200, 'ptsrc', 'plaw', bkgname_t200, ebins=ebins, emin=0.2)
    prepare_onoff(obsname_t300, 'ptsrc', 'plaw', bkgname_t300, ebins=ebins, emin=0.3)
    prepare_onoff(obsname_t700, 'ptsrc', 'plaw', bkgname_t700, ebins=ebins, emin=0.7)

    # Determine light curve for T700 selection (unbinned)
    par = set_pars(edisp=True)
    lightcurve(obsname_t700, bkgname_t700, 'ptsrc', 'fplaw700', '3D', par)

    # Determine light curve for T700 selection (On/Off)
    par = set_pars(edisp=True, ebins=10)
    lightcurve(obsname_t700, bkgname_t700, 'ptsrc', 'fplaw700', 'ONOFF', par)

    # Fit spectra for T200 selection (unbinned)
    par = set_pars(edisp=True)
    analyse(obsname_t200, bkgname_t200, 'ptsrc', 'plaw', par,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, bkgname_t200, 'ptsrc', 'eplaw', par,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, bkgname_t200, 'ptsrc', 'logparabola', par,
            fitname='pks_t200_results', npredname='pks_t200_npred')

    # Determine spectra for T200 selection (unbinned)
    par = set_pars(edisp=True)
    butterfly(obsname_t200, bkgname_t200, 'ptsrc', 'plaw', par, emin=0.2,
              fitname='pks_t200_results', buttername='pks_t200_butterfly')
    butterfly(obsname_t200, bkgname_t200, 'ptsrc', 'eplaw', par, emin=0.2,
              fitname='pks_t200_results', buttername='pks_t200_butterfly')
    butterfly(obsname_t200, bkgname_t200, 'ptsrc', 'logparabola', par, emin=0.2,
              fitname='pks_t200_results', buttername='pks_t200_butterfly')
    spectrum(obsname_t200,  bkgname_t200, 'ptsrc', par,
             emin=0.2, emax=7.6695, ebins=19,
             fitname='pks_t200_results', specname='pks_t200_spectrum')

    # Fit position for T300 & T700 selections (unbinned)
    par = set_pars(edisp=True)
    analyse(obsname_t300, bkgname_t300, 'ptsrc', 'plaw', par,
            fitname='pks_t300_results', npredname='pks_t300_npred')
    analyse(obsname_t300, bkgname_t300, 'ptsrc', 'eplaw', par,
            fitname='pks_t300_results', npredname='pks_t300_npred')
    analyse(obsname_t300, bkgname_t300, 'ptsrc', 'logparabola', par,
            fitname='pks_t300_results', npredname='pks_t300_npred')
    analyse(obsname_t700, bkgname_t700, 'ptsrc', 'plaw', par,
            fitname='pks_t700_results', npredname='pks_t700_npred')
    analyse(obsname_t700, bkgname_t700, 'ptsrc', 'eplaw', par,
            fitname='pks_t700_results', npredname='pks_t700_npred')
    analyse(obsname_t700, bkgname_t700, 'ptsrc', 'logparabola', par,
            fitname='pks_t700_results', npredname='pks_t700_npred')

    # Fit spectra for T200 selection (joint binned)
    par = set_pars(binned=True, ebins=ebins, edisp=True)
    analyse(obsname_t200, bkgname_t200, 'ptsrc', 'plaw', par,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, bkgname_t200, 'ptsrc', 'eplaw', par,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, bkgname_t200, 'ptsrc', 'logparabola', par,
            fitname='pks_t200_results', npredname='pks_t200_npred')

    # Fit spectra for T200 selection (stacked binned)
    par = set_pars(stacked=True, ebins=ebins, edisp=True)
    analyse(obsname_t200, bkgname_t200_stacked, 'ptsrc', 'plaw', par,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, bkgname_t200_stacked, 'ptsrc', 'eplaw', par,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, bkgname_t200_stacked, 'ptsrc', 'logparabola', par,
            fitname='pks_t200_results', npredname='pks_t200_npred')

    # Fit spectra for T200 selection (On/Off)
    pars_joint   = set_pars(onoff=True, ebins=ebins)
    pars_stacked = set_pars(onoff=True, stacked=True, ebins=ebins)
    analyse(obsname_t200, '',
            'ptsrc', 'plaw', pars_joint,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, 'model_pks_t200_onoff-joint40_cstat_ptsrc.xml',
            'ptsrc', 'plaw', pars_joint,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, '',
            'ptsrc', 'plaw', pars_stacked,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, 'model_pks_t200_onoff-stacked40_cstat_ptsrc.xml',
            'ptsrc', 'plaw', pars_stacked,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    #
    analyse(obsname_t200, '',
            'ptsrc', 'eplaw', pars_joint,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, 'model_pks_t200_onoff-joint40_cstat_ptsrc.xml',
            'ptsrc', 'eplaw', pars_joint,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, '',
            'ptsrc', 'eplaw', pars_stacked,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, 'model_pks_t200_onoff-stacked40_cstat_ptsrc.xml',
            'ptsrc', 'eplaw', pars_stacked,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    #
    analyse(obsname_t200, '',
            'ptsrc', 'logparabola', pars_joint,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, 'model_pks_t200_onoff-joint40_cstat_ptsrc.xml',
            'ptsrc', 'logparabola', pars_joint,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, '',
            'ptsrc', 'logparabola', pars_stacked,
            fitname='pks_t200_results', npredname='pks_t200_npred')
    analyse(obsname_t200, 'model_pks_t200_onoff-stacked40_cstat_ptsrc.xml',
            'ptsrc', 'logparabola', pars_stacked,
            fitname='pks_t200_results', npredname='pks_t200_npred')
