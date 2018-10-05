#!/usr/bin/env python3

"""
File name: prob3.py

    Python Version: 3.6

        Script to compare lagrange polynomial interpolation 
    and hermite polynomial interpolation of uniformly-spaced nodes.

LJ Brown
Fall 2018
"""

__filename__ = "prob3.py"
__author__ = "L.J. Brown"

# imports

# internal library
import logging

# my library
from lagrange import *
from hermite import *

# external library
import numpy as np
import matplotlib.pyplot as plt

# initilize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

#
# testing
#

if __name__ == "__main__":

    # general parameters
    L = 3.0
    nvals = [ 41, 21, 11, 5]            # polynomial degrees
    zvals = np.linspace(-L, L, 201)     # evaluation points

    #
    # test problem 1 
    #

    logger.info("Running test problem 1, f(x) = arctan(2x^2), f'(x) = (4x) / [ 4(x^4) + 1] \n")

    f = lambda x: np.arctan(2*x**2)
    df = lambda x: (4*x)/(4*(x**4) +1)

    # set data values
    fvals = f(zvals)

    # allocate storage for lagrange interpolant
    pl = np.zeros(zvals.shape, dtype=float)

    # allocate storage for hermite interpolant (evaluated using modified divided diffrence method)
    phdd = np.zeros(zvals.shape, dtype=float)

    # allocate storage for hermite interpolant (evaluated using tmp method)
    phtmp = np.zeros(zvals.shape, dtype=float)

    # run test for varying numbers of interpolation nodes
    for n in nvals:
        lagrange_nodes = np.linspace(-L, L, n+1)     # lagrange interpolation nodes 
        lagrange_f_yvals = f(lagrange_nodes)           # function values

        hermite_nodes = np.linspace(-L, L, int((n+1)//2)) 
        hermite_f_yvals = f(hermite_nodes)
        hermite_df_yvals = df(hermite_nodes)

        # evaluate interpolant lagrange
        pl = lagrange(lagrange_nodes, lagrange_f_yvals, zvals)

        # evaluate interpolant hermite using modified divided diffrence method
        phdd = hermite(hermite_nodes, hermite_f_yvals, hermite_df_yvals, zvals, method='modified_divided_difference')

        # evaluate interpolant hermite using modified lagrange method
        phl = hermite(hermite_nodes, hermite_f_yvals, hermite_df_yvals, zvals, method='modified_lagrange')

        # generate comparison plots
        fig, axarr = plt.subplots(1,2)

        #
        # overlayed interpolant (normal) plots
        #

        ftitle = '$f(z)$'

        ltitle = '$p_L(z)$'
        htitle = '$p_H(z)$'
        hltitle = '$p_H(z) (l)$'
        
        fplot = axarr[0].plot(zvals, fvals, 'r--', label=ftitle)

        lplot = axarr[0].plot(zvals, pl, 'b-', label=ltitle)
        hplot = axarr[0].plot(zvals, phdd, 'g-', label=htitle)
        #hlplot = axarr[0].plot(zvals, phl, 'y-', label=hltitle)

        axarr[0].scatter(hermite_nodes,hermite_f_yvals,s=10,color='g', zorder=2)

        axarr[0].set_xlabel('x')
        axarr[0].set_ylabel('y')
        axarr[0].set_title('Problem 1: $f$ vs $p_{' + str(n) + '}$')

        leg = axarr[0].legend(loc='upper center', shadow=True)

        #
        # overlayed semilogy plots
        #

        l_error_title = '$ | f(z) - p_L(z) | $'
        h_error_title = '$ | f(z) - p_H(z) | $'
        hl_error_title = '$ | f(z) - p_H(z) (l)| $'

        lplot_error = axarr[1].plot(zvals, abs(fvals-pl), 'b-', label=l_error_title)
        hplot_error = axarr[1].plot(zvals, abs(fvals-phdd), 'g-', label=h_error_title)
        #hlplot_error = axarr[1].plot(zvals, abs(fvals-phl), 'y-', label=hl_error_title)

        axarr[1].set_xlabel('x')
        axarr[1].set_ylabel('y')
        axarr[1].set_title('Problem 1: $E_{' + str(n) + '}$')

        leg = axarr[1].legend(loc='upper center', shadow=True)


    # display plots
    plt.show()
