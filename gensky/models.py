import numpy as np
import sys
import warnings
import os
import dummy_func

if not sys.warnoptions:
    warnings.simplefilter("ignore")

class photometricmodel:
    def __init__(self, time):
        self.time = time

    def v19(self, theta):
        """ Returns the parametric model from the Vilar et al. 2019 model with the continity descriptions from the Sanchez-Saez et al. (eqn A4)

        Input
        -----
        self.time (int): time of the transient in MJD (or phase?)
        t_rise (int): Rise time
        t_fall (int): Decline time
        t0 (int) [MJD]: "Start" Time
        A (int): Amplitude
        beta (int): Plateu slope
        t0-t1 (int): plateu duration

        Output
        ------
        model (float): flux modoel for the V19
        """
        # unpack my variables
        t_rise, t_fade, t0, A, beta, t1 = theta
        y = np.piecewise(self.time, [self.time<t1, self.time>=t1],
                        [lambda time: (A*(1-beta*((time-t0)/(t1-t0))))/(1+np.exp(-(time-t0)/t_rise)),
                        lambda time: (A*(1-beta)*np.exp(-(time-t1)/(t_fade)))/(1+np.exp(-(time-t0)/(t_rise)))]
                        )
        return (y)

    @staticmethod
    def cos_fouri(time, amp, index, phase, period, t0):
        chi = np.modf((time-t0)/period)[0]
        return amp*np.cos(2*np.pi*index*chi + phase)

    @staticmethod
    def sin_fouri(time, amp, index, phase, period, t0):
        chi = np.modf((time-t0)/period)[0]
        return amp*np.sin(2*np.pi*index*chi + phase)

    @staticmethod
    def sin_cos_combo(time, ai, bi, index, period, t0):
        c = (time-t0)/period
        chi = c - c.astype(int)
        return ai*np.cos(2*np.pi*index*chi) + bi*np.sin(2*np.pi*index*chi)

    def fourier_eb(self, theta, nfc=4):
        """Calculates the combined cosine and sine fourier series function.

           Input
           -----
           m0 (float): offset flux/magnitude of source
           ...
           ....

        """
        m0 = theta[0]
        Ai, Bi = theta[1], theta[2] # amplitudes
        P, T0 = theta[3], theta[4]

        # Evaluate model
        model = np.zeros(shape=np.shape(self.time)) # model
        for i in range(1, nfc+1):
            model += self.sin_cos_combo(self.time, Ai[i-1], Bi[i-1], i, P, T0)

        return m0 + model

    def fourier(self, theta, nfc=3, mode='sin'):
        """ Calculate the sine fourier series for N fourier terms (*modified to also include the period).

            Input
            -----
            m0 (float): average baseline value of function
            amps (float or list): list of amplitudes of each sine component
            phi (float or list): phases of each sine component
            T (float): Period of total sine curve

            Output
            ------
            fourier_cos (numpy.array): array of the fourier sine funcuntion (for N components)


            Example
            ------
            A typical Cepheid-like star would have the following parameters
            t = np.linspace(0, 10, 5200)

            Fourier_model = fourfam(t, 1,
                        np.array([1, 0.55, 1]), # amplitude list
                        np.array([1, 2.0, 1]), # k-index list
                        np.array([1, -4.5, 1]), # phase list
                         5) # Period of curve
        """

        # Unpack model parameters
        m0, amps, phases = theta[0], theta[1], theta[2]
        period, t0 = theta[3], theta[4]

        model = np.zeros(shape=np.shape(self.time)) # model
        if mode=='sin':
            for i in range(1, nfc+1):
                model += self.sin_fouri(self.time, amps[i-1], i, phases[i-1], period, t0)
        elif mode=='cos':
            for i in range(1, nfc+1):
                model += self.cos_fouri(self.time, amps[i-1], i, phases[i-1], period, t0)

        return (m0 + model)


    def sho(self, theta, noise=0.1):
        """ Return a simple harmonic oscillator function of the functional form:

            Y = A*sin(2piT/Period - phi ) + noise

            Input
            -----
            theta (list of floats): [amplitude, period, phase] of the simple harmonic osccilator
            noise (float): random gaussian noise

            Output
            -----
            sho(float): returns the values of a simple harmonic osccilator
        """

        A, period, phi = theta
        model = A*np.sin((2*np.pi*self.time)/period - phi)
        model += np.random.normal(0, np.sqrt(noise)) #add some random noise

        return model


# WIP:
# > Add some physical constraints to the model (i.e distance modulus)
