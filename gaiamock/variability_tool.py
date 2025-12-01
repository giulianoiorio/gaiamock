import numpy as np


class VariabilityTool:

    def __init__(self ,fchrom):
        self.fchrom=fchrom

    def __call__(self, time):
        return 0*time

class SimpleSinusoidal(VariabilityTool):
    
    def __init__(self, amp, period,fchrom):
        super().__init__(fchrom=fchrom)
        self.amp = amp
        self.period = period

    def __call__(self, time):
        return self.amp * np.sin(2 * np.pi * time / self.period)

