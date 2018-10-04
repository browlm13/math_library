#!/usr/bin/env python3

"""
# Script to test compare 
# lagrange polynomial interpolation and hermite polynomial interpolation
# of uniformly-spaced nodes.
#
# LJ Brown
# Math5315 @ SMU
# Fall 2018
"""

# imports

# my libraries
from lagrange import *
from hermite import *

# external libraries
import numpy as np
import matplotlib.pyplot as plt

#
# testing
#

if __name__ == "__main__":

    #
    # general parameters
    #

    L = 3.0
    nvals = [ 5, 11, 21, 41]            # polynomial degrees
    zvals = np.linspace(-L, L, 201)     # evaluation points

    ##
    ## test problem 1 ##
    ##

    print("Running test problem 1, f(x) = arctan(2x^2), f'(x) = (4x) / [ 4(x^4) + 1] ")

    f = lambda x: np.arctan(2*x**2)
    df = lambda x: (4*x)/(4*(x**4) +1)

    #f = lambda x: np.sin(x)
    #df = lambda x: np.cos(x)

    # set data values
    fvals = f(zvals)

    # allocate storage for lagrange interpolant
    pl = np.zeros(zvals.shape, dtype=float)

    # allocate storage for hermite interpolant
    ph = np.zeros(zvals.shape, dtype=float)

    # run test for varying numbers of interpolation nodes
    for n in nvals:
        lagrange_nodes = np.linspace(-L, L, n+1)     # lagrange interpolation nodes 
        lagrange_f_yvals = f(lagrange_nodes)           # function values

        #hermite_nodes = np.linspace(-L, L, int(n+1//2)) 
        hermite_nodes = np.linspace(-L, L, (n+1)//2) 
        hermite_f_yvals = f(hermite_nodes)
        hermite_df_yvals = df(hermite_nodes)

        # evaluate interpolant lagrange
        #for i in range(zvals.size):
        #    pl[i] = lagrange(lagrange_nodes, lagrange_f_yvals, zvals[i])
        pl = lagrange(lagrange_nodes, lagrange_f_yvals, zvals)

        # evaluate interpolant hermite
        ph = hermite(hermite_nodes, hermite_f_yvals, hermite_df_yvals, zvals)

        # generate comparison plots
        fig, axarr = plt.subplots(1,2)

        #
        # overlayed interpolant (normal) plots
        #
        axarr[0].plot(zvals, fvals, 'r--', zvals, pl, 'b-', zvals, ph, 'g-')

        axarr[0].set_xlabel('x')
        axarr[0].set_ylabel('y')
        axarr[0].set_title('Problem 1: $f$ vs $p_{' + str(n) + '}$')

        # set y_lim temporary ***
        miny, maxy = np.min(fvals), np.max(fvals)
        padding = 1
        axarr[0].set_ylim([miny - padding, maxy + padding])

        #
        # overlayed semilogy plots
        #
        axarr[1].plot(zvals, abs(fvals-pl), 'b-', zvals, abs(fvals-ph), 'g-')

        axarr[1].set_xlabel('x')
        axarr[1].set_ylabel('y')
        axarr[1].set_title('Problem 1: $E_{' + str(n) + '}$')

    # display plots
    plt.show()
