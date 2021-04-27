import numpy as np
from comms_filters import rrcosfilter
from fir_filter import FIRFilter

class TimingSyncMeyrOerder():
    '''
        A Meyr-Oerder algorithm based 
        timing recovery module
    '''

    def __init__(self, N=4):
        '''
            Create the object
        '''

        self.N = N

        # create filter coefficients
        self.L = 33
        self.M = 32
        self.alpha = 0.5
        self.h = rrcosfilter(N     = self.L,
                             alpha = self.alpha,
                             Ts    = 1,
                             Fs    = self.N)[1]
        self.filter = FIRFilter(self.h, complex=True)

        self.reset()


    def reset(self):
        '''
            Reset the internal state
        '''

        self.summation = 0
        self.k = 0


    def update(self, sample):
        '''
            Update a new sample in/out
        '''
        
        self.sample_filtered = self.filter.update(sample)/128.0

        squared = self.sample_filtered**4
        self.summation += squared * np.exp(-1j*2*np.pi*self.k/self.N)

        self.k += 1
        
        return np.angle(self.summation)/(2*np.pi)