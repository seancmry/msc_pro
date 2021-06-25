import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


def ackley(x, a=20, b=0.2, c=2*np.pi):

    d = len(x) # dimension of input vector x
    sum_sq_term = -a * np.exp(-b * np.sqrt(sum(x*x) / d))
    cos_term = -np.exp(sum(np.cos(c*x) / d))
    return a + np.exp(1) + sum_sq_term + cos_term

def plot_ackley_3d():
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make data.
    X = np.arange(-32, 32, 0.25)
    Y = np.arange(-32, 32, 0.25)
    X, Y = np.meshgrid(X, Y)

    a = 20
    b = 0.2
    c = 2 * np.pi

    sum_sq_term = -a * np.exp(-b * np.sqrt(X*X + Y*Y) / 2)
    cos_term = -np.exp((np.cos(c*X) + np.cos(c*Y)) / 2)
    Z = a + np.exp(1) + sum_sq_term + cos_term

    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                                   linewidth=0, antialiased=False)

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

def plot_ackley_2d():
    x = np.arange(-32, 33)
    y = [ackley(np.array([x_i])) for x_i in x]
    plt.plot(x,y)
    plt.show()


plot_ackley_2d()