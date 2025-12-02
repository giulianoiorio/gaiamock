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

## Implementation

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

### The VariabilityTool class 

#### VariabilityTool: framework for modelling variability-induced chromatic shifts

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
-   Implement the `__call__(time)` method returning a colour- or
    flux-variability signal evaluated at BJD(TCB) times

The goal is to provide a standardised way to inject colour variability
and study its effect on Gaia astrometric solutions---while ensuring
backward compatibility when no variability is applied.

------------------------------------------------------------------------

#### \### Core concept and behaviour

The base class:

``` python
class VariabilityTool:
    def __init__(self, fchrom, relative_norm=True)
    def __call__(self, time)
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

-   **`__call__(time)`**\
    This is the function that will be called in Gaiamock to get the colour at at given time of 
    Gaia observations, then it will be internally transformed to a shift based on the parameter 
    `fchrom` and `relative_norm`
    The base implementation returns zero variability, making
    `VariabilityTool` a no-op. Using the default instance
    `VariabilityTool(0.)` ensures perfect backward compatibility with
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
------------------------------------------------------------------------

#### 2. RRLVariable: RR Lyrae variability from Gaia SOS harmonic models

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
        # Compute the variability signal at BJD(TCB) times
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
    def __init__(self, amp, tau, t0, fchrom, relative_norm=True):
        super().__init__(fchrom=fchrom, relative_norm=relative_norm)
        self.amp = amp
        self.tau = tau
        self.t0 = t0

    def __call__(self, time):
        dt = time - self.t0
        return self.amp * np.exp(-np.clip(dt, 0, None) / self.tau)
```


