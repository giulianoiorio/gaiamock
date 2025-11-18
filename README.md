# This fork 

In this fork version we pull the recent changes from the upstream repository  + 
we added the new module gaiamock_mod. It contains all the function of gaiamock with 
the same name but handle the ruwe estimate for single stars in a more detailed way, 
see the main documentation below. Notice that to use it, we need additional files (see below).


# gaiamock

This is a package for simulating Gaia astrometry at the epoch level. To install it, do the following: 

(0) Install the required packages: numpy, matplotlib, os, ctypes, healpy, and joblib. In my experience these are usually relatively painless to install via pip. If you want to simulate sources distributed throughout the Galaxy with dust, you also need to install [mwdust](https://github.com/jobovy/mwdust), which is also available via pip. 

(1) clone this repository (click the green "code" button in the upper right of this page).

(2) download the file healpix_scans.zip from [this](https://caltech.box.com/s/4f7q6qdh0bku881bzvzxc4cm5u0902cf) link.
Unzip it inside gaiamock/ so that you have a directory gaiamock/healpix_scans/. That directory should have 49152 fits files inside it. The total size after unzipping will be 984 MB. 

(3) If you don't already have it installed, install [GSL](https://www.gnu.org/software/gsl/). On a Mac, this is likely most easily accomplished via Homebrew. 

(4) Inside the gaiamock/ directory, compile the file kepler_solve_astrometry.c. This will require linking GSL. The exact command will depend on where you installed GSL and on your compiler. On my Macbook, the command to compile was: 

gcc -shared -o kepler_solve_astrometry.so kepler_solve_astrometry.c -I/opt/homebrew/Cellar/gsl/2.7.1/include  -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas -lm 

On my local cluster (where GSL is already linked by default), the command to compile was:

gcc -shared -o kepler_solve_astrometry.so kepler_solve_astrometry.c -lgsl -lgslcblas -lm -fPIC 

If everything works, this will create a compiled file gaiamock/kepler_solve_astrometry.so. If it didn't work, there is probably a problem with the linking of your GSL installation.

(5) You are ready to go! Some basic functionality is demonstrated in the demo.ipynb notebook. 

A simple bash submission script to run on a cluster is provided in example_bash_submission.py. That reproduces Fig 13 of the paper. 

# A modified version to predict RUWE more reliably for small orbits and single stars

The default version of gaiamock does a pretty good job of predicting RUWE for binaries with "large" photocenter orbits -- see e.g. Figure 1 of [this](https://arxiv.org/abs/2504.11528) paper. However, it doesn't accurately predict the RUWE distribution of single stars or binaries with barely-detectable orbital motion: it predicts a RUWE distribution that is narrower than observed. There are at least two reason for this: 

(a) The default gaiamock bins (i.e. averages) the 8-9 measurements from individual CCDs during a single FOV transit for computational speed. This reduces the variance in the predicted RUWE due to shot noise.

(b) The observed Gaia data (at least in DR3) displays systematic trends in the median RUWE with sky position, likely in part due to crowding (e.g. the median RUWE for bright stars is lowest near the Galactic center). The origin of these trends is not well understood.

A modified version of gaiamock is available to improve the reliability of RUWE predictions. It does away with binning (and therefore is slower by a factor of a few, but still fast enough for most applications), and implements an empirical position-dependent rescaling of the epoch-level astrometric uncertainties. This is probably not the optimal way to model RUWE, but it is significantly more reliable than the default version of the code. 

To use the modified version of gaiamock, do the following: 

(1) From [this](https://caltech.box.com/s/lnszhrytqjt4f28l6eoyghsxrs7vg25n) link, download the individual_ccds.zip file and unzip the contents into the healpix_scans/ directory of gaiamock. You can leave the other files that are already there.

(2) Download the healpix_16_med_ruwe.npz file and put it in your gaiamock/ directory.

(3) Download the gaiamock_mod.py file and put it in your gaiamock/ directory.

Now you can 
```
import gaiamock_mod as gaiamock
```
and use the same functions you would use in gaiamock, e.g. for predicting epoch astrometry and computing RUWE. 

# A modified version to include the chromaticity effect due to star Variability

IN Gaia, the measured along-scan (AL) centroid of a star depends not only on its true astrometric motion but also on its spectral
energy distribution (SED).
Because Gaia PSF/LSF and optical response are wavelength-dependent, stars with different colours (BP-RP) produce slightly different
AL image centroids. 
This systematic shift is known as chromaticity. 

For non-variable stars, Gaia's calibraion pipeline remove the average chroamtic offest using colour information BP/RP or directly 
the specta. 
However, for stars with time-varying colours (e.g., RR Lyrae or Cepheids), the chromaticity will not change randomly but radially periodically, and the a-posteriori correction will not remove this effect (but just centre the periodic shift around the mean). 
This will likely increase the residual of the astrometic fit inflating the RUWE. 
Indeed, in Gaia for RR Lyrae there is the tendency to have inflated RUWE at high amplitudes ([Belokurov et al., 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.1922B/abstract))

As discussed in the Gaia thecnical document [Chromaticity in Gaia](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://dms.cosmos.esa.int/COSMOS/doc_fetch.php%3Fid%3D2694426&ved=2ahUKEwjslpWj9fqQAxXlVaQEHUg4AdUQFnoECB0QAQ&usg=AOvVaw19wUyqqsQqgIK0X1OY4fHR) the chromaticity effect depends mostly on the effective wavenumber $\nu_{eff}$ reported in Gaia in the columns pseudocolor. 
With the pseudocolor, the shight along the AL can be modeled as $\Delta \eta = C_1 (\nu_{eff}-\nu_{eff,0}) + C_2 (\nu_{eff}-\nu_{eff,0})^2 + ....$ with $\nu_{eff,0}$ that is a reference frequency.

The pseudocolor depends on the colour and can be approximated by 
$$
\nu_{eff} \approx 1.76 - \frac{1.61}{\pi} \mathrm{atan} \left( 0.531 (G_{BP} - G_{RP}) \right) \mu m^{-1}
$$
as in Eq. 3 of [Lindegren+21](https://ui.adsabs.harvard.edu/abs/2021A%26A...649A...2L/abstract).

Anyway, for a simple first order prescription of the chromaticity effect, we can assume that the offset is a linear function of the color shift:

$$
\eta = k_G ((BP-RP)- (BP-RP)_0) = k_G ((BP-RP)_0 + \delta_C) = k_G(BP-RP)_0 + k_G\delta_C
$$
where $(BP-RP)_0$ is the reference mean color and the first term can be excluded by the model because it is the term it is corrected for in the calibarion, so
$$
\Delta \eta = k_G \delta_c
$$, 
and since for variable stars the color change is time dependent, we have
$$
\Delta \eta (t) = k_G \delta_c(t,t_0)
$$
where $t_0$ is the reference epogh for the light curve. 

The value of the coefficient $k_G$ this is essentially a free paramter, however we can set it on the same order of the typical  random uncertantines on the calibrated $\eta$. This depends on G and can be find in Fig. 1 of [El-Badry+24](https://ui.adsabs.harvard.edu/abs/2024OJAp....7E.100E/abstract) (Fig. 1), based on the analysis in [Holl+23](https://ui.adsabs.harvard.edu/abs/2023A%26A...674A..10H/abstract) (their Fig. 3), that is relatated to the astrometric paper by [Lindegren+21](https://www.aanda.org/articles/aa/pdf/2021/05/aa39709-20.pdf). 

So the final model could be 

$$
\Delta \eta (t) = f_\sigma \sigma_\eta(G) \delta_c(t,t_0)
$$

Varying f, we can explore the effect in term of signal over the noise, typical color amplitude for RR Lyrae are  < 1 mag. >
