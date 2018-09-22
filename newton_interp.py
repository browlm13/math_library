#!/usr/bin/env python3

"""
    File name: newton_interp.py
    Python Version: 3.6

        Construct a newton interpolating polynomial for f(x) using linspaced data points for nodes.
        Where f(x*) = 0 is the root finding form of the fixed point problem function g(x*) = x(*).

                                    f(x*) = g(x*) - x* = 0.

                                        (Plot Results)

    TODO: use chebyshev nodes

    L.J. Brown
    Math5315 @ SMU
    Fall 2018
"""

__filename__ = "newton_interp.py"
__author__ = "L.J. Brown"

# external libraries
import numpy as np

def coeffients(x, y):
    """ 
        Computes and returns the coeffients of the interpolating polynomial of degree len(x).
        refrence sources: ['https://stackoverflow.com/questions/14823891/newton-s-interpolating-polynomial-python']

        :param x: 1d numpy array of x datapoints.
        :param y: 1d numpy array of f(x) datapoints.
        :returns: 1d numpy array of coeffiencts for newton interpolating polynomial.
    """

    # ensure floating point datatypes
    x.astype(float)
    y.astype(float)

    # degree of interpolating polynomial
    n = len(x)

    # intitilize list of coeffients for interpolating polynomial to y values
    c = y.tolist()

    # compute coeffients
    for j in range(1, n):
        for i in range(n-1, j-1, -1):
            c[i] = float(c[i]-c[i-1])/float(x[i]-x[i-j])

    # return an array of polynomial coefficient, note: reverse order for np.polyval function
    return np.array(c[::-1])

def testing_fixed_point_newton_interp(fixed_point_functions, n, m=400):
    """
        Plot/test accuracy interpolating polynomial construction using newton \
        interpolating polynomials.

        :param fixed_point_functions: dictonary of fixed point functions to test target \
                                     form: g(x*) = x*. converts to root finding problem \
                                     f(x*) = g(x*) = x* = 0.
        :param n: number of data points given for construction and degree of interpolating polynomial.
        :param m: number of datapoints to plot in display plot, default m=400.
    """

    # Function to convert to root finding problem given g(x). 'g(x*) = x*' -> 'f(x*) = 0'
    Ffun = lambda Gfun: lambda x: Gfun(x) -x

    import matplotlib.pylab as plt

    # setting up figure
    num_plots = len(fixed_point_functions)

    fig, axs = plt.subplots(1, num_plots, figsize=(15, 6), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace=.001)
    axs = axs.ravel()

    i = 0 # 'graph number'
    for Gfun_name, Gfun in fixed_point_functions.items():

        # <computation block>

        # convert to root finding problem
        f = Ffun(Gfun)

        # compute x and y data points
        x = np.linspace(-1,1,n)
        y = f(x)

        # compute coefficients of interpolating polynomial
        c = coeffients(x,y)

        # evaluate actual function points for graph
        ax = np.linspace(-1,1,m)
        ay = f(ax)

        # calculate y values using the interpolating polynomials coefficients
        y_hats = []
        for xi in ax:
            y_hati = np.polyval(c, xi)
            y_hats.append(y_hati)

        # <\computation block>

        # create plot for this function
        axs[i].plot( ax, ay, 'k' )     # function in black
        axs[i].plot( ax, y_hats, 'r' ) # interpolating polynomial in red
        axs[i].set_title(Gfun_name)

        # increment graph number
        i += 1

    plt.show()

# testing
if "__main__" in __name__:

    n = 10    # number of nodes to interpolate (uses linearly spaced x's for evaluations of f(x))
    m = 100   # number of datapoints to plot in display plot.

    # test functions:

    # use: 
    #    Gfun_a = test_functions['Gfun_a']
    #    f = Ffun(Gfun_a)

    fixed_point_functions = {

        'Gfun_a' : lambda x: (x**2)/4 -x/2 -1,
        'Gfun_b' : lambda x: np.cos(x),
        'Gfun_c' : lambda x: (x/3 +x)/2,
        'Gfun_d' : lambda x: np.cosh(x)/x -np.arctan(x)
    }

    testing_fixed_point_newton_interp(fixed_point_functions, n, m)

