import numpy as np

from timing_sync_polyphase import TimingSyncPolyphase
from timing_sync_meyroerder import TimingSyncMeyrOerder


class TimingRecovery():
    '''
        A Python model of timing recovery module.
        Contains multiple algorithms.
    '''

    def __init__(self,  k_p=0.5, k_i=1/2048):
        '''
            Create the object
        '''


        self.poly_obj = TimingSyncPolyphase(k_p=k_p, k_i=k_i)
        self.mo_obj   = TimingSyncMeyrOerder()

        self.reset()


    def reset(self):
        '''
            Reset the internal state
        '''

        self.poly_obj.reset()
        self.mo_obj.reset()

        self.state = 'acquistion'

        self.count = 0


    def update(self, sample, high_snr=True):
        '''
            Update the timing recovery object
        '''

        if self.state == 'acquistion':

            estimate = self.mo_obj.update(sample[0])

            if self.count == 4*8*4-1:
                self.state = 'tracking'
            else:
                self.count += 1
            
            self.poly_obj.controller_phase_index_ = (0.5-estimate)*self.poly_obj.M % self.poly_obj.M 


        return self.poly_obj.update(sample, high_snr=high_snr)
