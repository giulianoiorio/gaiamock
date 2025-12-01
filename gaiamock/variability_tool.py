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

