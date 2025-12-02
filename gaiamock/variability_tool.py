import numpy as np


class VariabilityTool:
    """
    Base class for variability tools.
    It should be subclassed to implement specific variability models, but 
    all the classes must specify fchrom and relative_norm attributes.
    fchrom: scaling factor to convert AL uncertainty to chromatic shift uncertainty.
    relative_norm: if True, normalisaiton is relative to the AL uncertainty per CCD at given G mag; if False, absolute normalisation (mas).
    If relative_norm is True, the bias due to the chromatict shift is the same at all G mag, because it scales with the G-dependent AL uncertainty.
    If relative_norm is False, the bias due to the chromatict shift is larger where the G-uncertainty are smaller (G-mag ), because the AL uncertainty is smaller there.
    """

    def __init__(self ,fchrom, relative_norm=True):
        """
        Initialize the variability tool.
        Parameters:
        fchrom: scaling factor to convert AL uncertainty to chromatic shift uncertainty.
        relative_norm: if True, normalisaiton is relative to the AL uncertainty per CCD at given G mag; if False, absolute normalisation (mas).
        """
        self.fchrom=fchrom
        self.relative_norm=relative_norm

    def __call__(self, time):
        """
        Evaluate the variability model at given times.
        Parameters: time: array-like, times at which to evaluate the model.
                    The input time  is in Barycentric Julian Date (BJD) in Barycentric Coordinated Time (TCB) standard.
                    This is th stadard time used in Gaia astrometry.
        Returns: array-like, variability values at the given times.
        """
        return 0*time

class SimpleSinusoidal(VariabilityTool):
    
    def __init__(self, amp, period,fchrom,relative_norm=True):
        super().__init__(fchrom=fchrom,relative_norm=relative_norm)
        self.amp = amp
        self.period = period

    def __call__(self, time):
        return self.amp * np.sin(2 * np.pi * time / self.period)

### FOr RR Lyrae

class LCtemplate:
    """
    Class to model light curves using a Fourier series template.
    """
    def __init__(self, zp_mag, A, Phi, narm, reference_time,period,
                obs_average):
        """
        Initialize the light curve template.
        Parameters:
        zp_mag: zero-point magnitude
        A: array-like, amplitudes of the Fourier series components
        Phi: array-like, phases of the Fourier series components
        narm: number of harmonics to use
        reference_time: reference time for phase calculation (BJD in TCB in days)
        period: period of the variability (in days)
        obs_average: if True, the light curve is averaged over the observation duration.
        """
        
        self.REF_TIME_SCALE = 2455197.5 #The reference time in the Gaia archive are reported as BJD in TCB - 2455197.5, but here we assume that the time is always in tBJD in TCB
        self.zp_mag = zp_mag
        self.narm = narm
        self.A    = A[:narm]
        self.Phi  = Phi[:narm]
        self.reference_time = reference_time
        self.period = period
        self.iarm = np.arange(1,narm+1)
        self.obs_average   = obs_average

    def lc(self,x,folded=False):
        """
        Evaluate the light curve at given times.
        Parameters:
        x: array-like, times at which to evaluate the light curve.
           The input time  is in Barycentric Julian Date (BJD) in Barycentric Coordinated Time (TCB) standard.
        folded: if True, x is treated as phase (0 to 1), otherwise as time in days (see above).
        Returns: array-like, light curve values at the given times.
        """
        if folded:
            phase=x
        else:
            phase=(x-self.REF_TIME_SCALE-self.reference_time)/self.period 
        
        return self.zp_mag + np.sum([ self.A[i]*np.cos(2*np.pi*self.iarm[i]*phase+self.Phi[i]) for i in range(self.narm) ],axis=0)
        

    def lc_der(self,x,folded=False):
        """
        Evaluate the derivative of the light curve at given times.
        Parameters:
        x: array-like, times at which to evaluate the derivative.
           The input time  is in Barycentric Julian Date (BJD) in Barycentric Coordinated Time (TCB) standard.
        folded: if True, x is treated as phase (0 to 1), otherwise as time in days (see above).
        Returns: array-like, derivative values at the given times.
        """
        if folded:
            phase=x
            dmod_dx = 2*np.pi*self.iarm
        else:
            phase=(x-self.REF_TIME_SCALE-self.reference_time)/self.period 
            dmod_dx = 2*np.pi*self.iarm/self.period
            
        return np.sum([ -self.A[i]*dmod_dx[i]*np.sin(2*np.pi*self.iarm[i]*phase+self.Phi[i]) for i in range(self.narm) ],axis=0)
    
class RRLVariable(VariabilityTool):
    """
    Class to model RR Lyrae variability using Gaia SOS light curve harmonic decomposition.
    It reconstructs the BP and RP light curves from the Gaia SOS parameters, and computes the
    color variability as the difference between BP and RP light curves removing the average color.
    Parameters:
    gaiadf: pandas Series, containing the Gaia SOS parameters for the RR Lyrae
    fchrom: scaling factor to convert AL uncertainty to chromatic shift uncertainty.    
    relative_norm: if True, normalisaiton is relative to the AL uncertainty per CCD at given G mag; if False, absolute normalisation (mas).
    """

    def __init__(self,gaiadf,fchrom,relative_norm=True):
        """
        Initialize the RR Lyrae variability model.
        Parameters:
        gaiadf: pandas Series, containing the Gaia SOS parameters for the RR Lyrae
        fchrom: scaling factor to convert AL uncertainty to chromatic shift uncertainty.
        relative_norm: if True, normalisaiton is relative to the AL uncertainty per CCD at given G mag; if False, absolute normalisation (mas).
        """
        #Initialtize the base class with fchrom and relative_norm
        super().__init__(fchrom=fchrom,relative_norm=relative_norm)

        #Reconstruct BP and RP light curve templates
        self.LCtemplate_bp = self.reconstruct_lc(gaiadf,band="bp")
        self.LCtemplate_rp = self.reconstruct_lc(gaiadf,band="rp")
        self.LCtemplate_g = self.reconstruct_lc(gaiadf,band="g")

        self.average_g = self.LCtemplate_g.obs_average
        self.average_bp = self.LCtemplate_bp.obs_average
        self.average_rp = self.LCtemplate_rp.obs_average
        self.average_color = self.average_bp - self.average_rp
        

    def __call__(self, time):
        """
        Evaluate the color variability (BP-RP) at given times.
        Parameters:
        time: array-like, times at which to evaluate the color variability.
              The input time  is in Barycentric Julian Date (BJD) in Barycentric Coordinated Time (TCB) standard.
        Returns: array-like, color variability values at the given times which average to zero.
        """
        return (self.LCtemplate_bp.lc(time)  - self.LCtemplate_rp.lc(time))-self.average_color

    def g_lcurve(self,time):
        """
        Evaluate the G-band light curve at given times.
        Parameters:
        time: array-like, times at which to evaluate the G-band light curve.
              The input time  is in Barycentric Julian Date (BJD) in Bary
              Centric Coordinated Time (TCB) standard.
        Returns: array-like, G-band light curve values at the given times."""
        return self.LCtemplate_g.lc(time)
    
    def bp_lcurve(self,time):
        """
        Evaluate the BP-band light curve at given times.
        Parameters:
        time: array-like, times at which to evaluate the BP-band light curve.
              The input time  is in Barycentric Julian Date (BJD) in Barycentric Coordinated Time (TCB) standard.
        Returns: array-like, BP-band light curve values at the given times."""
        return self.LCtemplate_bp.lc(time)

    def rp_lcurve(self,time):
        """
        Evaluate the RP-band light curve at given times.
        Parameters:
        time: array-like, times at which to evaluate the RP-band light curve.
              The input time  is in Barycentric Julian Date (BJD) in Barycentric Coordinated Time (TCB) standard.
        Returns: array-like, RP-band light curve values at the given times."""
        return self.LCtemplate_rp.lc(time)
    
    
    @staticmethod
    def reconstruct_lc(df,band="g"):
        """
        In Gaia SOS the RRL light curve are fitted with a Fourier series as
        G(t) = zp + Sum A_i cos(2pi i nu_max (t-tref) + Phi_i)
        where zp is the zero poing magnitude, A_i and Phi_i are the amplitude 
        of the armonics, with nu_max the frequencey of the fundamental one with nu_max=1/Period
        """
        zp_mag = df[f"zp_mag_{band}"] #zero point for the Fourier curve decomposition
        numax=df["fund_freq1"] #numax (1/Period)
        narm=df[f"num_harmonics_for_p1_{band}"]
        reference_time = df[f"reference_time_{band}"] #Reference time for the fourier curve decomposition
        Ai=df[f"fund_freq1_harmonic_ampl_{band}"][:narm]
        Phii=df[f"fund_freq1_harmonic_phase_{band}"][:narm]
        iarm = np.arange(1,narm+1)
        obs_average = df[f"phot_{band}_mean_mag"] #Use this one instead of the one of the SOS because I think this is the one used in postprocessing for astometric calibration
    
        return LCtemplate(zp_mag=zp_mag, A=Ai, Phi=Phii, narm=narm, 
                          reference_time=reference_time, period=1/numax,
                         obs_average=obs_average)