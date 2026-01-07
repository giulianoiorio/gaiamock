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

# A modified version to include a comprehensive treatment of variable stars in the astrometric analysis.

The main purpose of this fork is to include variable stars and their effects or biases on astrometric solutions.

Here, we consider two possible configurations: a variable star in isolation, and a variable star in a binary system in which the companion is not variable. We remain agnostic about the source of variability, so any type of variability can be included in the analysis, including irregular behavior such as flares

The effect of a variable star on astrometric analysis involves three main contributions (ordered from least to most likely relevant):

- **Heteroskedastic astrometric errors**:  
  Since Gaia astrometric and photometric uncertainties depend on the \(G\) magnitude, a variable source will have different errors at different transits. This effect is negligible for sources with $12 \lesssim G \lesssim 15$, where the error curve is nearly flat, but it may become relevant for very bright or very faint sources, where the errors vary rapidly with magnitude.

- **Variability-Induced Mover (VIM)**:  
  In binary systems, Gaia traces the motion of the photocentre. If one component is variable, changes in its flux will shift the photocentre position. This produces an astrometric signal that can be detectable even when the binary orbital period is too long to be directly detected through orbital motion alone. More details can be found in [Halbwachs+23](http://arxiv.org/abs/2206.05726).

- **Chromatic shift**:  
  Gaia’s optics are not perfectly achromatic, meaning that the position of a star on the detector depends on its spectral energy distribution (i.e. its color). The Gaia calibration pipeline accounts for this effect assuming a constant source luminosity (see [Chromaticity in Gaia](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://dms.cosmos.esa.int/COSMOS/doc_fetch.php%3Fid%3D2694426&ved=2ahUKEwjslpWj9fqQAxXlVaQEHUg4AdUQFnoECB0QAQ&usg=AOvVaw19wUyqqsQqgIK0X1OY4fHR)). As a result, random photometric variability is largely absorbed into the astrometric noise. However, variable sources that also exhibit color variations will introduce residual, time-dependent patterns, leading to additional time-correlated astrometric noise.

We have created a new module called *gaiamock_var* in which we have updated all the function to take into account al the three effects.
To use this version:
```
import gaiamock_var as gaiamock
```

This is almost full retrocompatible (except for the function `predict_astrometry_single_source` and `predict_astrometry_luminous_binary`) now returning two additional arrays (see [here](#heteroskedastic-astrometric-errors)). Most of the relevant function now accept an extra parameter, that is a class including the variability model (see  [VariabilityTool class description](#the-variabilitytool-class)). However in all the functions there is a default with a non-variable model, therefore all the codes using standard *gaiamock* should continue to work. 
To properly include variability effects, the user must define a variability model (see below).

## Heteroskedastic astrometric errors

We have added the epoch dependent estimate of errors in the function 

- `predict_astrometry_single_source`
- `predict_astrometry_luminous_binary`
- `predict_astrometry_and_rvs_simultaneously`

In addition to the module *gaiamock_var*  this has been added also to the same functions in the module *gaiamock_mod*.

### Implementation 

Very simply, at each observation epoch t, we have

$$G(t) =G_\mathrm{mean} + dG(t)$$

where $G_\mathrm{mean}$ is the mean magnitude G observed in Gaia (the one stored in the variable *phot_g_mean_mag*) and $dG(t)$ is the scaled light curve used to model the variability (see  [VariabilityTool class description](#the-variabilitytool-class)).

So the epoch erros is $\sigma_\eta=f_{\sigma_{eta}}(G(t))$. The functions now also estimate the photometric errors in the same way $\sigma_G=f_{\sigma_{G}}(G(t))$ (based on tool described [here](https://www.cosmos.esa.int/web/gaia/fitted-dr3-photometric-uncertainties-tool)).
The two functions also estimate the observed G, $G_\mathrm{obs}$ convolving the true $G(t)$ with the error $\sigma_G(t)$, and with respect to the same functions in the original *gaiamock* it returs two additional arrays containg $G_\mathrm{obs}$ and $G_\mathrm{err}$.


## VIM: Astrometric shift induced by variability in binary 

### VIM astrometric bias

As described in [Halbwachs+23](http://arxiv.org/abs/2206.05726), 
in the case of a binary system hosting a variable source, the source variability causes a change in the flux ratio of the two sources changing the along-scan position of the photocentre, the so called "variability induced movers" (VIM). 
This effect if actually searched for in the astrometric binary analysis in Gaia DR3 with a dedicated model.  
Concerning the gaiamock, this effect can be easily introduced 
by simply estimating the flux luminosity at each scanning time, rather than using a simple value. 

The flux luminosity is defined such as (see [El-Badry+24](https://ui.adsabs.harvard.edu/abs/2024OJAp....7E.100E/abstract))

$$
f =  \frac{L_2}{L_1} = 10^{-\frac{G_2 -G_1}{2.5}} =  10^{\frac{G_1 -G_2}{2.5}}, 
$$
where the subscript 1 refers to the most luminous star, so that we have always $0 \leq f \leq1$. 
If one of the two stars (let's assume the most luminous one) is a variable star, the current  version of Gaiamock uses the average G-band magnitude, so
$$
f  = 10^{\frac{\langle G_1 \rangle -G_2}{2.5}}, 
$$. 

The true flux ratio will be instead 
$$
f_\mathrm{var}(t) = 10^{\frac{G_1(t) -G_2}{2.5}} = 10^{\frac{\langle G_1 \rangle -G_2}{2.5}} 10^{\frac{G_1(t) - \langle G_1 \rangle}{2.5}} =  
f 10^{\frac{d G_1(t)}{2.5}}
$$

So the final time dependent flux ratio is equal to the standard averaged one times a correction that depends of the magnitude variation rescaled to have mean magnitude 0. 
In case of large photometric variations and/or stars with similar magnitude, the role of the primary and secondary can switch. 

#### Implementation 

We further modify  the module **gaiamock_var**, updating the following functions:

- `predict_astrometry_luminous_binary`
- `predict_astrometry_and_rvs_simultaneously`

The most important change is that while before the flux ratio was constant now depends on the epoch (see above). In addition the phase dependent photometric errors are also consideres same as described in the [previous section](#heteroskedastic-astrometric-errors). Also in this case the functions returns two additional arrays containg $G_\mathrm{obs}$ and $G_\mathrm{err}$.

In case of large photometric variations and/or stars with similar magnitude, the role of the primary and secondary can switch, in this case $f>1$ and we switch primary and secondary by estimating the corrected flux ratio as $1/f$ and the corrected mass ratio as $1/q=m_1/m_2$.

We have modified the same functions also in  `gaiamock_mod`. 

### VIMF model fit

A model to account for the VIM effect is among the one used to check for astrometric binary in GaiaDR3 ([Halbwachs+23](http://arxiv.org/abs/2206.05726)). It assumes that the binary motion is negligible and only the VIM is present, for this reason it is defined as fixed VIM (VIMF).
The model is quite simple, the astrometric along-scan shift is:

$$
\eta(t) = \left[ \Delta \alpha + \mu_\alpha t + D_\alpha  \left( \frac{\bar{F_G}}{F_G(t)}-1 \right) \right] \sin (\psi) + \left[ \Delta \delta + \mu_\delta t + D_\delta  \left( \frac{\bar{F_G}}{F_G(t)}-1 \right) \right] \cos (\psi) + \Pi \omega
$$

with ($\Delta \alpha, \mu_\alpha, \Delta \delta, \mu_\delta, \omega$) the usual 5 paramters of the single star solution (2 positions, 2 proper motions, parallax), and $\vec{D}=(D_\alpha, D_\delta)$ is the parameter to fit to account for the the magnitude and direction of the VIM.
The ratio $\bar{F_G}{F_G(t)}$ is the inverse of the flux ratio between the epoch flux in the G-band and a reference one, that we assume is the mean related to the measured *phot_g_mean_mag*. To rewrite it as function of the magnitude and its variation we note that 
$$\frac{\bar{F_G}}{F_G(t)} = 10^\frac{G_\mathrm{obs}(t)-G_\mathrm{mean}}{2.5} $$

In this case however, the data errors are not only dependent on the astrometric error but also on the photometric error. Hence the final error is

$$
\sigma(t) = \sqrt{\sigma^2_\eta(t) + \sigma^2_\mathrm{mod}(t)}, 
$$
where from the propagation of errors 
$$
\sigma_\mathrm{mod} = \sigma_\mathrm{F} \frac{\bar{F_G}}{F^2_G(t)} | D_\alpha \sin (\psi)  + D_\delta \cos (\psi) |, 
$$

same for the propagation of errors we can write

$$
\sigma_F = \frac{\ln 10}{2.5} * \sigma_G * 10**(-G/2.5)= \frac{\ln 10}{2.5}  \sigma_G  10^{-G/2.5}=  \frac{\ln 10}{2.5} \sigma_G  F
$$, so

$$
\sigma_\mathrm{mod} = \frac{\ln 10}{2.5} \sigma_\mathrm{G} \frac{\bar{F_G}}{F_G(t)} | D_\alpha \sin (\psi)  + D_\delta \cos (\psi) | = \frac{\ln 10}{2.5} \sigma_\mathrm{G} 10^\frac{G_\mathrm{obs}(t)-G_\mathrm{mean}}{2.5} | D_\alpha \sin (\psi)  + D_\delta \cos (\psi) | , 
$$. 

Since the final errors depends on the fitting parameters ($D_\alpha ,D_\delta$), the an iterative procedure must be considered to find the best parameters. 


#### Implementation 

We have added the function 
 - `check_VIMF` 
 
 following the implementations of the other check function (in particular of `check_7par` since they share the same number of paramters). The difference with respect to the other check functions is that it requests two additional paramters containing the epoch photometry: $G_\mathrm{obs}$ and $G_\mathrm{err}$. These are present in the output of both `predict_astrometry_luminous_binary` and 
and `predict_astrometry_single_source`.

Using the epoch photometry, the function gets a first estimate of 
$D$ and $\sigma_\mathrm{mod}$, then it checks if at all the epochs $\sigma_\mathrm{mod}<0.01 \sigma_\eta$, if this is the case this first solution is accepted as the main one, otherwise it stars an iterative fitting procedure updating each time $D$ and $\sigma_\mathrm{mod}$, until both components of $D$ converge within 1% or the number of iterations reach 5.

The final output of the method includes:

- $F_2$ statistic (goodness of fit, Eq. 1 in [Halbwachs+23](http://arxiv.org/abs/2206.05726))
- $s$ statistic (signficance, Eq. 3 in [Halbwachs+23](http://arxiv.org/abs/2206.05726))
- 7 best-fit parameters with:
    - Position offset along ra
    - Proper motion along ra
    - $D$ component along ra 
    - Position offset along dec 
    - Proper motion along dec
    - $D$ component along dec
    - parallax 
- 7 errors pf the best-fit parameters (order as above).

Given the implementation of the new astrometric binary model, we have integrated it in the functions:

- `fit_full_astrometric_cascade`
- `run_full_astrometric_cascade`

In the first function there are now two additional parameters $G_\mathrm{obs}$ and $G_\mathrm{err}$ that by default are None. If both of them are not None and include the epoch photometry (output of the `predict_astrometry_luminous_binary function`), then the last step of the Gaia astrometric cascade is activated. 
This last model is the VIMF and was not included in the standard *gaiamock*. 

The VIMF model is checked only if the full binary solution is not accepted. Then to accept or not the VIFM, we implement the criteria discussed in [Halbwachs+23](http://arxiv.org/abs/2206.05726):

- $F_2<25$
- $s>12$
- $\frac{\omega}{\sigma_\omega}>30$ (parallax over parallax error)

**Note** In  [Halbwachs+23](http://arxiv.org/abs/2206.05726), the condition $s>12$ is used just as a pre-condition to accept the solution, but then no solutions with $s<20$ are included in the final catalogue
We still use 12 to be consistent with the rest of the astrometric cascade (also the acceleration and variable acceleration are treated in the same way in Gaia, but in *gaiamock* a threshould of 12 is used).

## Chromaticity effect

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

So the final model is

$$
\Delta \eta (t) = f_\sigma \sigma_\eta(G) \delta_c(t,t_0)
$$

Varying f, we can explore the effect in term of signal over the noise, typical color amplitude for RR Lyrae are  0.1-0.5

As an alternative, we can assume a constant  normalisation constant that does not depend on G, so

$$
\Delta \eta (t) = f_\sigma \delta_c(t,t_0)
$$

In thi case the effect of the chromatic shfit will be larger for bright source with small astrometric  uncertanties and becomes negligible for faint source.

### Implementation

We added the module **gaiamock_var**, which includes the same functions as the standard `gaiamock` module,  
but with variations in the following functions:

- `predict_astrometry_single_source`
- `predict_astrometry_luminous_binary`
- `predict_astrometry_binary_in_terms_of_a0`
- `run_only_5par_solution`
- `run_full_astrometric_cascade`

These functions now include an additional parameter, `variability_tool`, which must be an instance of the class  
`VariabilityTool` defined in the module `variability_tool.py`.  
By default, this parameter is set to the dummy instance `VariabilityTool(0.)`, which introduces **no variability-induced bias**.  
Therefore, running these functions without specifying a `variability_tool` instance yields **exactly the same results** as the standard `gaiamock` functions, ensuring full **backwards compatibility**.

In addition to the functions in `gaiamock_var`, the function  
`predict_astrometry_single_source` in the module `gaiamock_mod` also accepts this extra parameter.

To study the bias introduced by intrinsic variability, one must define a specialised class that inherits from  
`VariabilityTool`, or use one of the classes already implemented in `variability_tool.py`.  
Currently available classes include:

- `SimpleSinusoidal` — models variability following a simple sinusoidal pattern  
- `RRVariable` — reproduces RR Lyrae colour variations using the Gaia SOS harmonic decomposition

The section below describes how to use the available classes and how to implement new ones.


## The VariabilityTool class 

### VariabilityTool: framework for modelling variability-induced chromatic shifts

The `VariabilityTool` framework provides a unified and extensible
interface to model **intrinsic photometric variability** and its impact
on **chromaticity-driven astrometric bias** in Gaia-like simulations.\
It is designed to be integrated into the `gaiamock_var` module and any
astrometric prediction routine that accepts a `variability_tool`
parameter.

All variability models implemented within this framework must:

-   Inherit from the base class `VariabilityTool`
-   Define the attributes
    -   `fchrom`: scaling factor converting AL uncertainty into chromatic shift.
    -   `relative_norm`: boolean flag controlling the normalisation scheme
-   Implement the `__call__(self,time)` method returning a colour
    variability signal evaluated at BJD(TCB) times (rescaled for the colour average).
-   Implement the `g_lcurve_normalised(self,time)` method returning the 
    flux-variability signal evaluated at BJD(TCB) times (rescaled for the magnitude average).

The goal is to provide a standardised way to inject colour variability
and study its effect on Gaia astrometric solutions---while ensuring
backward compatibility when no variability is applied.

------------------------------------------------------------------------

#### Core concept and behaviour

The base class:

``` python
class VariabilityTool:
    def __init__(self, fchrom, relative_norm=True)
    def __call__(self, time)
    def g_lcurve_normalised(self,time)
```

implements the following logic:

-   **`fchrom`**\
    A multiplicative factor converting a color variation of the source
    to the chromatic shfit.
     Depending on the `relative_norm` value (see below) it could be a constant of units of mas/mag (`relative_norm=true`) or be interepreted as a strengh of the chromaticity sfhit in terms of astrometric noise at given G (`relative_norm=false`).

-   **`relative_norm=True`**\
    If enabled, the chromatic bias is normalised *relative to the AL
    uncertainty* at the star's G magnitude.\
    This ensures a *magnitude-independent* bias amplitude: fainter and
    brighter stars receive variability-driven offsets scaled to their
    intrinsic Gaia AL uncertainties.

-   **`relative_norm=False`**\
    The chromatic bias is expressed in absolute units (mas).\
    The induced astrometric shift becomes larger for bright stars (with
    smaller AL uncertainties) and smaller for faint stars. This value corresponds to the expected shift for a colour variability (with respect to the mean) of 1 mag

-   **`__call__(self,time)`**\
    This is the function that will be called in Gaiamock to get the colour at at given time of 
    Gaia observations, then it will be internally transformed to a shift based on the parameter 
    `fchrom` and `relative_norm`
    The base implementation returns zero variability, making
    `VariabilityTool` a no-op. Using the default instance
    `VariabilityTool(0.)` ensures perfect backward compatibility with
    the standard `gaiamock` behaviour.
    **Note**, when called within Gaiamock the time will be in BJD in TCB. 

- **`g_lcurve_normalised(self,time)`**\
    This is the function that will be called in Gaiamock to get the flux  at at given time of  Gaia observations, then it will be used to estimate the flux ratio and then along-scan position in combination with the shift due to the binary motion.
    The base implementation returns zero variability, making
    `VariabilityTool` a no-op. Using the default instance
    or not overloading the method ensures perfect backward compatibility with
    the standard `gaiamock` behaviour.
    **Note**, when called within Gaiamock the time will be in BJD in TCB. 

Typical usage:

``` python
from variability_tool import VariabilityTool
vt = VariabilityTool(fchrom=0.0)
signal = vt(time_array)   # Always zero
```

------------------------------------------------------------------------

### Implemented Variability Models

#### 1. SimpleSinusoidal

``` python
class SimpleSinusoidal(VariabilityTool)
```

A minimal model producing a **pure sinusoidal** colour or
flux variation.

**Parameters:**

-   `amp`: amplitude of the sinusoid
-   `period`: period in days
-   `fchrom`: chromatic-scaling factor
-   `relative_norm`: whether the normalisation is relative or absolute


**Example:**

``` python
from variability_tool import SimpleSinusoidal
vt = SimpleSinusoidal(amp=0.05, period=0.6, fchrom=1.0)
bias = vt(time_bjd_tcb)
```

#### 2. Burst-like signal

``` python
class ColorBurst(VariabilityTool)
```

A class to model a burst like increase (or decrease) of color
it is modelled as a split Gaussian with: 

**Parameters:**

-   `amp`: (BP-RP) value at the peak (positive if the burst is reddening the source, negative otherwise)
-   `tpeak`: time of the peak 
-   `dtrise`: time scale for the rising part (decreasing part if amp<0)
-   `dtdecay`: time scale for the decaying part (incresing part if amp<0)
-   `fchrom`: chromatic-scaling factor
-   `relative_norm`: whether the normalisation is relative or absolute

``` python
from variability_tool import ColorBurst
#JBDTCB ≈ 2456863.5 to 2459000.5 #Range of days in GaiaDR3
tpeak_tcbjd = 2457063.5
cb=vt.ColorBurst(amp=-1,tpeak=2457063,dtrise=5,dtdecay=30,fchrom=1)
bias = vt(time_bjd_tcb)
```


------------------------------------------------------------------------

#### 3. RRLVariable: RR Lyrae variability from Gaia SOS harmonic models

``` python
class RRLVariable(VariabilityTool)
```

A specialised model that reconstructs the **BP**, **RP**, and **G**
light curves of an RR Lyrae star using the **Gaia SOS harmonic
decomposition**.\
From this, it computes the **BP--RP colour variability** as:

$$
v(t) = (BP(t) - RP(t)) - \langle  BP(t) - RP(t) \rangle
$$


**How it works:**

-   Uses SOS parameters to build a Fourier-series template for BP, RP, and G.
-   Computes the instantaneous colour variation at any BJD(TCB) time.
-   Removes the mean colour so that `<v(t)> = 0`, as required for
    chromatic modelling.

**Expected `gaiadf` content (per SOS RRL entry):**

-   `zp_mag_<band>`
-   `fund_freq1`
-   `fund_freq1_harmonic_ampl_<band>`
-   `fund_freq1_harmonic_phase_<band>`
-   `num_harmonics_for_p1_<band>`
-   `reference_time_<band>`
-   `phot_<band>_mean_mag`

**Additional methods:**

-   `g_lcurve(time)` → G-band magnitude at time `t`
-   `bp_lcurve(time)` → BP-band magnitude at time `t`
-   `rp_lcurve(time)` → RP-band magnitude at time `t`

**Example usage:**

``` python
from variability_tool import RRLVariable
vt = RRLVariable(gaiadf=sos_row, fchrom=1.0, relative_norm=False)
color_variation = vt(time_bjd_tcb)
Gmag_variation = vt.g_lcurve_normalised(time_bjd_tcb)
```

This class is ideal forsimulations of RR Lyrae astrometric biases in Gaia-like pipelines.

------------------------------------------------------------------------

### How to Implement a New Variability Model

To define a new variability model, create a subclass of
`VariabilityTool`, call the `super().__init__(fchrom=fchrom, relative_norm=relative_norm)`
in the init to initiliase correctly fchrom and relative norm
and override the `__call__` method.

A minimal template:

``` python
from variability_tool import VariabilityTool
import numpy as np

class MyCustomVariability(VariabilityTool):
    def __init__(self, fchrom, relative_norm=True, <your parameters>):
        super().__init__(fchrom=fchrom, relative_norm=relative_norm)
        # store your parameters here

    def __call__(self, time):
        # Compute the colour variability signal at BJD(TCB) times
        # Must return an array-like with the same length as time
        return <your expression>

    def g_lcurve_normalised(self, time):
        # Compute the magnitude variability signal at BJD(TCB) times
        # Must return an array-like with the same length as time
        return <your expression>
```

### Guidelines for writing a new model

1.  **Input time standard must be BJD(TCB)**\
    This is the native Gaia astrometric timescale used throughout
    `gaiamock_var`.

2.  **Relative variation around the mean**, the model 
    assumes that the Gaia astrometric calibration already removes the shift due to the average colour, 
    therefore the implemented colour variation should subtract the mean or better should subtract directly the reported colour mean in Gaia. 

3.  **`fchrom` must always be passed to the base class**\
    This ensures  conversion from variability → chromatic
    astrometric shift.

4.  **Prefer vectorised NumPy expressions**\
    Your implementation should efficiently handle numpy arrays of transit
    times.


### Example: exponential flare model

``` python
class FlareVariability(VariabilityTool):
    def __init__(self, amp_g, amp_color, tau, t0, fchrom, relative_norm=True):
        super().__init__(fchrom=fchrom, relative_norm=relative_norm)
        self.amp_g = amp_g
        self.amp_color = amp_color
        self.tau = tau
        self.t0 = t0

    def __call__(self, time):
        dt = time - self.t0
        return self.amp_color * np.exp(-np.clip(dt, 0, None) / self.tau)

    def g_lcurve_normalised(self, time):
        # Compute the magnitude variability signal at BJD(TCB) times
        # Must return an array-like with the same length as time
        dt = time - self.t0
        return self.amp_g * np.exp(-np.clip(dt, 0, None) / self.tau)
```


