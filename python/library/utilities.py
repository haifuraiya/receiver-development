from scipy.interpolate import interp1d

def phase_shift(input, phase_shift):
    '''
        Phase shift the input signal by a fractional number
    '''

    # create the original time steps
    x_orig = [_ for _ in range(len(input))]
    interpolate = interp1d(x_orig, input, kind='cubic')

    # split the phase shift into integer and fractional components
    phase_shift_int = int(phase_shift)
    phase_shift_frac = phase_shift % 1

    # create a new timestep list
    x_new = [_ + phase_shift_frac for _ in x_orig[phase_shift_int:-(phase_shift_int+1)]]

    # calculate new values
    return interpolate(x_new)
