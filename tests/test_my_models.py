import matplotlib.pyplot as plt
import sys
sys.path.append('../gensky')
from models import photometricmodel
import numpy as np

def plot_v19():
    times = np.linspace(-50, 50, 10000)

    M = photometricmodel(times).v19([8, 23, 0, 500, 0, 0])

    plt.figure(figsize=(5,5))
    plt.plot(times, M)
    plt.show()


def plot_sho():
    times = np.linspace(-2*np.pi, +2*np.pi, 100000)
    M = photometricmodel(times).sho([20, 3, 0.3])

    plt.plot(times, M)
    plt.show()

def plot_sine():
    t = np.linspace(0, 10, 5200)
    F = photometricmodel(t).fourier_sine(1,
            np.array([1, 0.55, 1]), # amps
            np.array([1, 2.0, 1]), # k
            np.array([1, -4.5, 1]), 5)
    plt.plot(t, F)
    plt.show()



def main():
    """ display all our models"""
    plot_v19()
    plot_sine()
    plot_sho()

if __name__ == "__main__":
    main()
