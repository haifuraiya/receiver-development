from comms_filters import rrcosfilter


class TimingSyncPolyphase():
    '''
        A Python model of the polyphase based timing
        synchronisation
    '''

    def __init__(self, N=4, M=32, L=33, alpha=0.5):
        '''
            Create the object
        '''

        # save values internally
        self.N      = N
        self.M      = M
        self.L      = L
        self.alpha  = alpha

        # create filter coefficients
        self.h_int = rrcosfilter(N     = self.L * self.M,
                                 alpha = self.alpha,
                                 Ts    = 1,
                                 Fs    = self.N * self.M)[1]
        
        self.dh_int = []
        for n in range(len(self.h_int)-1):
            self.dh_int.append(self.h_int[n+1]-self.h_int[n])

        # create the multiple phases
        self.h = []
        self.dh = []
        for m in range(self.M):
            self.h.append( [ self.h_int[self.M*n + m] for n in range(self.L-1) ] )
            self.dh.append([self.dh_int[self.M*n + m] for n in range(self.L-1) ] )



        
