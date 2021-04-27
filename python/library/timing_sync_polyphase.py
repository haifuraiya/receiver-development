
import numpy as np
import matplotlib.pyplot as plt
from comms_filters import rrcosfilter
from fir_filter import FIRFilter
from generic_modem import GenericModem
from utilities import phase_shift
import control

class TimingSyncPolyphase():
    '''
        A Python model of the polyphase based timing
        synchronisation
    '''

    def __init__(self, N=4, M=32, L=33, alpha=0.5, k_p=1.00, k_i=0.125, D=1):
        '''
            Create the object
        '''

        # save values internally
        self.N      = N         # samples per symbol
        self.M      = M         # number of phases
        self.L      = L         # filter length
        self.alpha  = alpha     # RRC factor
        self.k_p    = k_p       # proportional gain
        self.k_i    = k_i       # integral gain
        self.D      = D         # excess delay

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

        # create the polyphase filters
        self.MF = []
        for i in range(self.M):
            self.MF.append(FIRFilter(self.h[i], complex=True))
        
        self.dMF = []
        for i in range(self.M):
            self.dMF.append(FIRFilter(self.dh[i], complex=True))

        # reset the controller state
        self.reset()


    def sensitivity(self, high_snr=True):
        '''
            Calculate the error detector sensitivty
        '''

        # create test data
        data_length = 256
        data = [np.random.randint(2) for _ in range(data_length)]
        modem_bpsk_obj = GenericModem(  modulation_type = ["129", "QPSK 2/9"],
                                        samples_per_symbol = 4,
                                        pulse_factor = 0.5,
                                        pulse_length = 33,
                                        filename = "../library/DVB-S2X_constellations.json")
        tx = modem_bpsk_obj.modulate(data)

        # calculate the energy in the signal for normalisation
        energy = 0.0
        for n in tx:
            energy += abs(n)**2

        # measure sensitivity
        sensitivity = []
        phase     = np.linspace(0,self.N,512)
        phase_rad = np.linspace(0,2*np.pi,512)
        for phase_i in phase:
            MF  = FIRFilter(self.h[int(self.M/2)], complex=True)
            dMF = FIRFilter(self.h[int(self.M/2)], complex=True)
            filtered = []
            for n, samp in enumerate(phase_shift(tx,phase_i)):
                MF_samp  = MF.update(samp)
                dMF_samp = dMF.update(samp)
                if n % 4 == 0:
                    if high_snr:
                        filtered.append( np.sign(np.real(MF_samp)) * np.real(dMF_samp) 
                                            + np.sign(np.imag(MF_samp)) * np.imag(dMF_samp) )
                    else:
                        filtered.append( np.real(MF_samp) * np.real(dMF_samp)
                                            + np.imag(MF_samp) * np.imag(dMF_samp) )
            sensitivity.append(abs(sum(filtered))/energy)

        return [sensitivity, phase_rad] 


    def reset(self):
        '''
            Reset the controller state
        '''

        self.controller_i_sum        = 0
        self.controller_phase_index  = int(self.M/2)
        self.controller_sample_count = 0
        self.error                   = 0
        self.control                 = 0
        self.under_over_flow         = 0
        self.delay                   = [0.0 for _ in range(self.D)]



    def update(self, input_samples, high_snr=True):
        '''
            Update the control loop of the 
        '''

        output_samples = []

        for input_sample in input_samples:
        
            output  = []
            doutput = []

            # skip this sample
            if self.under_over_flow == 1:
                samples_number_times = 0

            # stuff an extra sample into the stream
            elif self.under_over_flow == -1:
                samples_number_times = 2

            # update normally
            else:
                samples_number_times = 1


            for i in range(samples_number_times):

                for i in range(self.M):
                    output.append(  self.MF[i].update(input_sample)  )
                    doutput.append( self.dMF[i].update(input_sample) )

                # add the excess delay
                if self.D > 0:
                    self.delay.append(self.controller_phase_index)
                    index = int(self.delay[0])
                    del self.delay[0]
                else:
                    index = int(self.controller_phase_index)

                # select the correct filter
                output_sel_real  = np.real( output[index])
                doutput_sel_real = np.real(doutput[index])
                output_sel_imag  = np.imag( output[index])
                doutput_sel_imag = np.imag(doutput[index])

                # handle the error differently depending on the SNR
                if high_snr:
                    self.error = np.sign(output_sel_real) * doutput_sel_real + np.sign(output_sel_imag) * doutput_sel_imag
                else:
                    self.error = output_sel_real * doutput_sel_real + output_sel_imag * doutput_sel_imag

                # update the controller
                if self.controller_sample_count != self.N-1:
                    self.controller_sample_count += 1
                else: 
                    self.controller_sample_count = 0

                    # update the PI controller
                    self.controller_i_sum += self.k_i * self.error;
                    self.control = self.k_p * self.error + self.controller_i_sum

                    # select new phase
                    self.controller_phase_index = (self.controller_phase_index - self.control)

                    # wrap around
                    if self.controller_phase_index >= self.M:
                        self.controller_phase_index -= self.M
                        self.under_over_flow = 1
                    elif self.controller_phase_index < 0:
                        self.controller_phase_index += self.M
                        self.under_over_flow = -1
                    else:
                        self.under_over_flow = 0

                # create the output samples
                output_samples.append(output_sel_real + 1j*output_sel_imag)

                if self.under_over_flow == 1:
                    self.under_over_flow = 0

        return output_samples


    def transfer_function(self, high_snr=True):
        '''
            Determine the transfer functionm of the loop
        '''

        # create z term
        z = control.TransferFunction.z
        z_1 = 1/z

        # find the maximum phase detector gain
        [sensitivity, phase] = self.sensitivity(high_snr)

        # k_pd = max( [abs(sensitivity[_+1] - sensitivity[_]) for _ in range(len(sensitivity)-1)] )/(phase[1]-phase[0])
        # k_pd = np.mean( [abs(sensitivity[_+1] - sensitivity[_]) for _ in range(len(sensitivity)-1)] )/(phase[1]-phase[0])
        k_pd = np.mean( sensitivity )

        # gain of each component
        H_pd  = k_pd
        H_p   = self.k_p
        H_i   = self.k_i * z_1 / (1 - z_1)
        H_nco = (2*np.pi/self.M) * z_1 / (1 - z_1)
        H_d   = z_1**self.D

        H_ol  = H_pd * (H_p + H_i) * H_nco * H_d

        H_ol = H_ol.minreal()

        return H_ol


    def bandwidth(self, high_snr=True):
        '''
            Determine the bandwidth of the loop
        '''

        # find the maximum phase detector gain
        [sensitivity, phase] = self.sensitivity(high_snr)

        # k_pd = max( [abs(sensitivity[_+1] - sensitivity[_]) for _ in range(len(sensitivity)-1)] )/(phase[1]-phase[0])
        k_pd = np.mean( sensitivity )

        H_pd  = k_pd
        H_p   = self.k_p
        H_nco = 2*np.pi/self.M

        return H_pd * H_p * H_nco


    def stability(self, high_snr=True):
        '''
            Determine the stability of the loop
        '''

        H_ol = self.transfer_function(high_snr)

        gm, pm, wg, wp = control.margin(H_ol)
        print("Margins")
        print("  Gain\t:%f"  % gm)
        print("  Phase\t:%f" % pm)

        rlist, klist = control.root_locus(H_ol, plot=True)
        plt.show()

