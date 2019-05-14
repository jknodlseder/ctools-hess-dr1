# script to prepare Femi LAT data
# requires Fermi LAT Science Tools
from GtApp import GtApp
import os

# define the env variable FERMI_DATA to point to the directory where you have Fermi data
# and FERMI_DIFFUSE_DIR to the directory where you have the Fermi diffuse models

# input files
ft1 = '@' + os.environ['FERMI_DATA'] + 'ft1_files.txt'  # prepend @ to txt file list of FT1 files
ft2 = os.environ['FERMI_DATA'] + 'L181201100111FD084EAD23_SC00.fits' # FT2 (spacecraft file)
diffmodel = 'diffmodel.xml'

# load Science Tools
selEv = GtApp('gtselect')
selGti = GtApp('gtmktime')
calcLt = GtApp('gtltcube')
binEv = GtApp('gtbin')
expMap = GtApp('gtexpcube2')
srcMap = GtApp('gtsrcmaps')

# analysis configuration
# event selection
evclass = 128
evtype = 3
irfs = 'P8R3_SOURCE_V2'  # must match evclass and evtype
emin = 5.e4  # MeV
emax = 1.e6
zmax = 105
# binning
nebins = 16
npix = 60
binsize = 0.05
coordsys = 'CEL'
xref = 83.633
yref = 22.015
proj = 'CAR'

# event selection
filter_evts = 'filter_evts.fits'

selEv.run(evclass=evclass, evtype=evtype, infile=ft1, outfile=filter_evts, ra='INDEF',
          dec='INDEF', rad='INDEF', tmin='INDEF', tmax='INDEF', emin=emin, emax=emax,
          zmax=zmax)

# GTI selection
filter_gti_evts = 'filter_gti_evts.fits'

selGti.run(scfile=ft2, filter='(DATA_QUAL>0)&&(LAT_CONFIG==1)', roicut='no',
           evfile=filter_evts, outfile=filter_gti_evts)

# livetime cube
ltcube = 'ltcube.fits'

calcLt.run(evfile=filter_gti_evts, scfile=ft2, outfile=ltcube, dcostheta=0.025, binsz=1,
           zmax=zmax)

# count cube
cntmap = 'cntmap.fits'

binEv.run(algorithm='CCUBE', evfile=filter_gti_evts, scfile=ft2,
          outfile=cntmap, ebinalg='LOG', emin=emin,
          emax=emax, enumbins=nebins, nxpix=npix, nypix=npix,
          binsz=binsize, coordsys=coordsys, xref=xref, yref=yref,
          axisrot=0, proj=proj)

# exposure map
expmap = 'expmap.fits'

expMap.run(infile=ltcube, cmap='none', outfile=expmap,
           irfs=irfs, nxpix=3 * npix, nypix=3 * npix, binsz=binsize, coordsys=coordsys,
           xref=xref, yref=yref, axisrot=0., proj=proj, ebinalg='LOG', emin=emin,
           emax=emax, enumbins=nebins)

# diffuse source maps
srcmap = 'srcmaps.fits'

srcMap.run(scfile=ft2, expcube=ltcube, cmap=cntmap, srcmdl=diffmodel, bexpmap=expmap,
           outfile=srcmap, irfs='CALDB', binsz=" ")