import numpy as np
import sys
import warnings
import os

if not sys.warnoptions:
    warnings.simplefilter("ignore")

class photometricmodel:
    def __init__(self, time):
        self.time = time

    def v19(self, theta):
        """ Returns the parametric model from the Vilar et al. 2019 model with the continity descriptions from the Sanchez-Saez et al. (eqn A4)

        Input
        -----
        time (int): time of the transient in MJD (or phase?)
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
