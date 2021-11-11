import matplotlib.pyplot as plt
from .models import photometricmodel


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

plot_v19()
plot_sho()
