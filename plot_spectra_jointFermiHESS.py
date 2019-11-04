import matplotlib.pyplot as plt
import gammalib
from astropy.table import Table
from astropy.io import fits
import astropy.units as u
import numpy as np

TeV2erg = 1.60218

twofhlfilename = "/Users/ltibaldo/Fermi/Catalogs/LAT/2FHL/gll_psch_v09.fit"
hessfilename = '/Users/ltibaldo/Software/GitHub/gamma-cat/input/data/2006/2006A%26A...457..899A/tev-000025-sed.ecsv'
magicfilename = '/Users/ltibaldo/HESS/DR1/Fermi_extended/MAGIC_2015_Crab_Nebula.fits'

tfhl = fits.getdata(twofhlfilename, 1)
tfhl_ebounds = np.array([0.05, 0.171, 0.585, 2.])
tfhl_emean = np.sqrt(tfhl_ebounds[:-1] * tfhl_ebounds[1:])

def sed_from2fhl(srcid):
    flux = np.zeros(3)
    flow = np.zeros(3)
    fhigh = np.zeros(3)
    flux[0] = tfhl[srcid]['Flux50_171GeV']
    flux[1] = tfhl[srcid]['Flux171_585GeV']
    flux[2] = tfhl[srcid]['Flux585_2000GeV']
    flow[0] = - tfhl[srcid]['Unc_Flux50_171GeV'][0]
    flow[1] = - tfhl[srcid]['Unc_Flux171_585GeV'][0]
    flow[2] = - tfhl[srcid]['Unc_Flux585_2000GeV'][0]
    fhigh[0] = tfhl[srcid]['Unc_Flux50_171GeV'][1]
    fhigh[1] = tfhl[srcid]['Unc_Flux171_585GeV'][1]
    fhigh[2] = tfhl[srcid]['Unc_Flux585_2000GeV'][1]
    flux *= TeV2erg * tfhl_emean ** 2 / (tfhl_ebounds[1:] - tfhl_ebounds[:-1])
    flow *= TeV2erg * tfhl_emean ** 2 / (tfhl_ebounds[1:] - tfhl_ebounds[:-1])
    fhigh *= TeV2erg * tfhl_emean ** 2 / (tfhl_ebounds[1:] - tfhl_ebounds[:-1])
    return flux, flow, fhigh

def plot_butterfly(filename, ax, color,label):

    csv = gammalib.GCsv(filename)

    # Initialise arrays to be filled
    butterfly_x = []
    butterfly_y = []
    line_x = []
    line_y = []

    # Loop over rows of the file
    nrows = csv.nrows()
    for row in range(nrows):
        # Get conversion coefficient
        conv = csv.real(row, 0) * csv.real(row, 0) * gammalib.MeV2erg

        # Compute upper edge of confidence band
        butterfly_x.append(csv.real(row, 0) / 1.0e6)  # TeV
        butterfly_y.append(csv.real(row, 2) * conv)

        # Set line values
        line_x.append(csv.real(row, 0) / 1.0e6)  # TeV
        line_y.append(csv.real(row, 1) * conv)

    # Loop over the rows backwards to compute the lower edge of the
    # confidence band
    for row in range(nrows):
        index = nrows - 1 - row
        conv = csv.real(index, 0) * csv.real(index, 0) * gammalib.MeV2erg
        butterfly_x.append(csv.real(index, 0) / 1.0e6)
        low_error = max(csv.real(index, 3) * conv, 1e-26)
        butterfly_y.append(low_error)

    ax.plot(line_x, line_y, color=color, ls='-',label=label)
    ax.fill(butterfly_x, butterfly_y, color=color, alpha=0.4)

def plot_spectrum(filename, ax, color):

    # Read spectrum file
    fits     = gammalib.GFits(filename)
    table    = fits.table(1)
    c_energy = table['Energy']
    c_ed     = table['ed_Energy']
    c_eu     = table['eu_Energy']
    c_flux   = table['Flux']
    c_eflux  = table['e_Flux']
    c_ts     = table['TS']
    c_upper  = table['UpperLimit']

    # Initialise arrays to be filled
    energies    = []
    flux        = []
    ed_engs     = []
    eu_engs     = []
    e_flux      = []
    ul_energies = []
    ul_ed_engs  = []
    ul_eu_engs  = []
    ul_flux     = []

    # Loop over rows of the file
    nrows = table.nrows()
    for row in range(nrows):

        # Get Test Statistic, flux and flux error
        ts    = c_ts.real(row)
        flx   = c_flux.real(row)
        e_flx = c_eflux.real(row)

        # If Test Statistic is larger than 9 and flux error is smaller than
        # flux then append flux plots ...
        if ts > 9.0 and e_flx < flx:
            energies.append(c_energy.real(row))
            flux.append(c_flux.real(row))
            ed_engs.append(c_ed.real(row))
            eu_engs.append(c_eu.real(row))
            e_flux.append(c_eflux.real(row))

        # ... otherwise append upper limit
        else:
            ul_energies.append(c_energy.real(row))
            ul_flux.append(c_upper.real(row))
            ul_ed_engs.append(c_ed.real(row))
            ul_eu_engs.append(c_eu.real(row))

    # Set upper limit errors
    yerr = [0.6 * x for x in ul_flux]

    # Plot the spectrum
    ax.errorbar(energies, flux, yerr=e_flux, xerr=[ed_engs, eu_engs],
                 fmt=color+'o',alpha=0.8)
    ax.errorbar(ul_energies, ul_flux, xerr=[ul_ed_engs, ul_eu_engs],
                 yerr=yerr, uplims=True, fmt=color+'o',alpha=0.8)


fig = plt.figure('Spectra', (7, 4))
fig.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.15)
ax = plt.subplot()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("Energy (TeV)", fontsize=14)
ax.set_ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)', fontsize=14)
ax.set_ylim(5.e-13, 2.e-10)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.grid()

# my analysis
plot_butterfly('butterfly_unbinned.txt',ax,'r','ctools: joint Fermi LAT + H.E.S.S. (unbinned)')
plot_spectrum('spectrum_unbinned.fits',ax,'r')

# HESS
table = Table.read(hessfilename, format='ascii.ecsv')

energies = table['e_ref'].quantity
flux = table['dnde'].quantity
flux_errp = table['dnde_errp'].quantity
flux_errn = table['dnde_errn'].quantity
flux = flux * energies**2
flux_errp = flux_errp * energies**2
flux_errn = flux_errn * energies**2
flux = flux.to(u.Unit("erg cm-2 s-1"))
flux_errp = flux_errp.to(u.Unit("erg cm-2 s-1"))
flux_errn = flux_errn.to(u.Unit("erg cm-2 s-1"))

ax.errorbar(energies.base,flux.base,xerr=0,yerr=[flux_errn.base,flux_errp.base],fmt='bo',alpha=0.8)

energies = np.logspace(np.log10(0.41),np.log10(40),50)
flux = 3.76e-11 * np.power(energies,-2.39) * np.exp(-energies/14.3)
flux = flux * energies**2
flux *= 1.6021765 # TeV 2 erg

ax.plot(energies,flux,color='b',alpha=0.8, linestyle='--',
        label='H.E.S.S. 2006 A&A 457 899A')

# LAT 2FHL
energies = np.logspace(np.log10(0.05),np.log10(2),50)
flux = 1.31e-9 * -1.13 * np.power(energies,-2.13) / (2**-1.13-0.05**-1.13)
flux = flux * energies**2
flux *= TeV2erg # TeV 2 erg


ax.plot(energies,flux,color='g',alpha=0.8, linestyle='--',
        label='Fermi LAT 2016 ApJS 222 1')

srcid = 85
flux, flow, fhigh = sed_from2fhl(srcid)
ax.errorbar(tfhl_emean, flux, yerr=[flow, fhigh],
            xerr=0,
            fmt='go', alpha=0.8)

#MAGIC
energies = np.logspace(np.log10(0.05),np.log10(30),100)
flux = 3.23e-11 * np.power(energies,-2.47 - 0.24 * np.log10(energies))
flux = flux * energies**2
flux *= TeV2erg # TeV 2 erg

ax.plot(energies,flux,color='m',alpha=0.7, linestyle='--',
        label='MAGIC 2015 JHEA 5-6 30-38')

spectrum = fits.getdata(magicfilename,'SPECTRUM')
energies = 1.e-3 * spectrum['energy']
flux = TeV2erg * spectrum['flux']
flux_err = TeV2erg * spectrum['Dflux']

ax.errorbar(energies,flux,xerr=0,yerr=flux_err,fmt='mo',alpha=0.7)

ax.legend()

fig.savefig("joint_fermi-hess_comparison.png",dpi=300)

fig1 = plt.figure('Spectra ctols', (7, 4))
fig1.subplots_adjust(left=0.16, right=0.95, top=0.95, bottom=0.15)
ax1 = plt.subplot()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel("Energy (TeV)", fontsize=14)
ax1.set_ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)', fontsize=14)
ax1.set_ylim(5.e-13, 2.e-10)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)
ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')
ax1.grid()

plot_butterfly('butterfly_unbinned.txt',ax1,'C1','ctools: joint Fermi LAT + H.E.S.S. (unbinned)')
plot_spectrum('spectrum_unbinned.fits',ax1,'C1')

plot_butterfly('butterfly_wstat.txt',ax1,'C0','ctools: joint Fermi LAT + H.E.S.S. (On-Off wstat)')
plot_spectrum('spectrum_wstat.fits',ax1,'C0')

ax1.legend()

fig1.savefig("joint_fermi-hess_unbinned_wstat.png",dpi=300)

#plt.show()


