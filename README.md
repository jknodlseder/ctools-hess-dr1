ctools analysis of the H.E.S.S. public data release
===================================================

This repository contains all scripts that were used to generate the results
of the paper "ctools analysis of the H.E.S.S. public data release" written
by J. Knödlseder, L. Tibaldo, D. Tiziani, A. Specovious, J. Cardenzana,
M. Mayer, N. Kelley-Hoskins, L. De Venere, S. Bonnefoy, A. Ziegler, S. Eschbach,
P. Martin, T. Louge, F. Brun, M. Haupt, and R. Bühler and published in A&A.


Prerequisites
-------------

In order to run the scripts in this repository, please set the following
environment variables

```
$ export HESSDATA=/path/where/the/hess/data/reside
$ export HESSGIT=/path/where/the/scripts/reside
```

You may of course add this setup to your ``.profile`` or ``.bashrc`` file to
avoid doing this manually all the time.


Empty field observations
------------------------

To generate a plot of the lookup table (Figure 1) execute

```
$ show_bkg_lookup.py
```

This script will generate the lookup table file ``bkg_lookup.fits`` and
display a smoothed version of that file.

To generate the plot of spatial gradients (Figure 2) execute

```
$ show_off_fits.py
```

This script will generate the lookup table file ``bkg_lookup.fits`` and use
this file for the maximum likelihood analysis of each empty field observation.
The fit results for all empty field observations are then gathered in a single
plot. The fit results are also stored in the file ``show_off_fits.dmp`` which
is simply read in when running the script again. If the fits should be repeated,
remove the ``show_off_fits.dmp`` file.

To generate the residual plots (Figures 3, 4 and 5) for all empty field
observations stacked together, as well as the residual plots for the 5 random
samples of 10 stacked empty field observations (Figures E.1 to E.15) execute 

```
$ show_stacked_offs_residuals.py
```

If in addition you want to see the residual plots for each individual
observation then execute

```
$ show_off_residuals.py
```


Crab nebula observations
------------------------

To analyse the Crab observations, run

```
$ analyse_crab.py
```

in a dedicated directory. The directory will be populated by lots of files
which contain the results of the analysis. Please note that the execution
of the script may take a long time, since the analysis comprises all methods
available in ctools and takes the energy dispersion into account.

To generate the Crab map (Figure 6), run

```
$ show_crab_map.py
```

while will display the file ``crab_resmap.fits``.


The source SED (Figure 7) is generated using

```
$ show_crab_spectrum.py crab_spectrum_ptsrc_plaw_lookup_grad_hess_edisp.fits crab_butterfly_ptsrc_eplaw_lookup_grad_hess_edisp.txt
```

The source extension map (Figure 8) is generated using

```
$ show_crab_extension.py
```

The residual plots (Figures 9, E.16 and E.17) are generated using

```
$ show_analysis_residuals.py Crab
```


MSH 15-52 observations
----------------------

To analyse the MSH 15-52 observations, run

```
$ analyse_msh.py
```

in a dedicated directory. The directory will be populated by lots of files
which contain the results of the analysis. Please note that the execution
of the script may take a long time, since the analysis comprises all methods
available in ctools and takes the energy dispersion into account.

The MSH 15-52 sky map with fitted morphology (Figure 10) is generated using

```
$ show_msh_extension.py
```

The residual plots for the power law spectrum (Figures 11, E.18 and E.19) are
generated using

```
$ show_analysis_residuals.py MSH-PL
```

and for the exponentially cut-off power law spectrum (Figures 12, E.20 and E.21)
are generated using

```
$ show_analysis_residuals.py MSH
```

The source SED (Figure 13) is generated using

```
$ show_msh_spectrum.py msh_spectrum_egauss_plaw_lookup_grad_hess_edisp.fits msh_butterfly_egauss_plaw_lookup_grad_hess_edisp.txt msh_butterfly_egauss_eplaw_lookup_grad_hess_edisp.txt
```


RX J1713.7-3946 observations
----------------------------

To analyse the RX J1713.7-3946 observations, run

```
$ analyse_rx.py
```

in a dedicated directory. The directory will be populated by lots of files
which contain the results of the analysis. Please note that the execution
of the script may take a long time, since the analysis comprises all methods
available in ctools and takes the energy dispersion into account.

The source SED (Figure 14) is generated using

```
$ show_rx_spectrum.py rx_spectrum_map0.48_plaw_lookup_grad_hess_edisp.fits rx_butterfly_map0.48_eplaw_lookup_grad_hess_edisp.txt
```

The residual plots (Figures 15, E.22 and E.23) are generated using

```
$ show_analysis_residuals.py RX
```

The background subtracted sky map (left panel of Figure 16) is generated using

```
$ show_map.py -s 0.06 rx_resmap.fits
```

The ring-background sky map (right panel of Figure 16) is generated using

```
$ show_map.py -s 0.06 rx_skymap.fits
```


PKS 2155-304 observations
-------------------------

To analyse the PKS 2155-304 observations, run

```
$ analyse_pks.py
```

in a dedicated directory. The directory will be populated by lots of files
which contain the results of the analysis. Please note that the execution
of the script may take a long time, since the analysis comprises all methods
available in ctools and takes the energy dispersion into account.

The map of PKS 2155-304 with the position fitting results (Figure 17) is
generated using

```
$ show_pks_position.py
```

The PKS 2155-304 light curves (Figure 18) are generated using

```
$ show_pks_lightcurve.py
```

The PKS 2155-304 light curve correlation plots (Figure 19) are generated using

```
$ show_pks_correlations.py
```

The source SED (Figure 20) is generated using

```
$ show_pks_spectrum.py pks_t200_spectrum_ptsrc_plaw_lookup_grad_hess_edisp_19bins.fits pks_t200_butterfly_ptsrc_logparabola_lookup_grad_hess_edisp.txt
```

The residual plots (Figures 21, E.24 and E.25) are generated using

```
$ show_analysis_residuals.py PKS
```


Joint Fermi-LAT and H.E.S.S. analysis
-------------------------------------

To prepare the joint Fermi-LAT and H.E.S.S. analysis run

```
$ prepare_fermi_data.py
```

Note that you need to define the ``FERMI_DIFFUSE_DIR`` environment variable
to point to the directory where the Fermi diffuse models are.

Now you can run

```
$ joint_hess_fermi.py
```

The joint SED plots (Figure 22) are generated using

```
$ plot_spectra_jointFermiHESS.py
```

(you need to adapt some file paths in this script to your local environment).
