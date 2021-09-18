import math
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

#sphere
def objective(x1,x2):
    return x1**2 + x2**2

#ackley
def ackley(x, a=20, b=0.2, c=2*np.pi):
    # x = input vals
    d = len(x) # dimension of input vector x
    sum_sq_term = -a * np.exp(-b * np.sqrt(sum(x*x) / d))
    cos_term = -np.exp(sum(np.cos(c*x) / d))
    return a + np.exp(1) + sum_sq_term + cos_term


def plot_ackley_3d():
    # Plot figure
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    # Make data.
    X = np.arange(-32.8, 32.8, 0.1)
    Y = np.arange(-32.8, 32.8, 0.1)
    X, Y = np.meshgrid(X, Y)

    a = 20
    b = 0.2
    c = 2 * np.pi

    sum_sq_term = -a * np.exp(-b * np.sqrt(X*X + Y*Y) / 2)
    cos_term = -np.exp((np.cos(c*X) + np.cos(c*Y)) / 2)
    Z = a + np.exp(1) + sum_sq_term + cos_term

    surf = ax.plot_surface(X, Y, Z, cmap='viridis',
                           linewidth=0, antialiased=False)
    ax.set_title('Ackley')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_zlabel('Fitness');
    ax.view_init(20, 30)

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()


def plot_rosenbrock_3d():
    # Function define
    b = 10
    f = lambda x,y: (x-1)**2 + b*(y-x**2)**2

    # Plot figure
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    # Make data.
    X = np.arange(-2.048, 2.048, 0.15)
    Y = np.arange(-2.048, 2.048, 0.15)
    X, Y = np.meshgrid(X, Y)
    Z = f(X, Y)

    # Plot the surface
    surf = ax.plot_surface(X, Y, Z, cmap='viridis',
                          linewidth=0, antialiased=False)
    ax.set_title('Rosenbrock')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_zlabel('Fitness');
    ax.set_zlim(0, 500)
    ax.view_init(20, 30)
    fig.colorbar(surf, shrink=0.5, aspect=10)

    plt.show()

def plot_sphere_3d():

    # Make data.
    X = np.arange(-100, 100, 0.15)
    Y = np.arange(-100, 100, 0.15)
    x1, x2 = np.meshgrid(X, Y)

    #Compute function
    result = objective(x1,x2)

    # Plot figure
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    ax.view_init(20, 30)

    # Plot the surface
    surf = ax.plot_surface(x1, x2, result, cmap='viridis',
                           linewidth=0, antialiased=False)
    ax.set_title('Sphere')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_zlabel('Fitness');
    fig.colorbar(surf, shrink=0.5, aspect=10)

    plt.show()


def plot_rastrigin_3d():

    X = np.linspace(-2, 2, 50)
    Y = np.linspace(-2, 2, 50)
    X, Y = np.meshgrid(X, Y)

    #NOTE: the function used in this instance is a simplified
    #version of the Rastrigin function to illustrate the fine grained view
    #of how the Griewank function appears on the hypercube ONLY. X and Y are small dimensions
    #used for the evaluation of this function (usually between 5,-5), but the Griewank function
    #is typically evaluated between -600,600, so at a far higher dimensionality.

    Z = (X ** 2 - 10 * np.cos(2 * np.pi * X)) - \
    (Y ** 2 - 10 * np.cos(2 * np.pi * Y)) + 1

    fig = plt.figure()
    ax = plt.axes(projection='3d')

    ax.view_init(20, 30)

    # Plot the surface
    surf = ax.plot_surface(X, Y, Z, cmap='viridis',  rstride=1, cstride=1, linewidth=0.08,
                antialiased=True)
    ax.set_title('Rastrigin')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_zlabel('Fitness');
    fig.colorbar(surf, shrink=0.5, aspect=10)

    plt.show()


#Call plots
plot_ackley_3d()
plot_rosenbrock_3d()
plot_sphere_3d()
plot_rastrigin_3d()
