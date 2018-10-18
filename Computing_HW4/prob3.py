#!/usr/bin/env python3

"""
File name: prob3.py

    Python Version: 3.6

        Script to compare hermite polynomial interpolation 
        and cubic spline interpolationof uniformly-spaced nodes.

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
from cubic_spline import *
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
    zvals = np.linspace(-L, L, 800)     # evaluation points

    #
    # test problem 1 
    #

    logger.info("Running test problem 1, f(x) = arctan(2x^2), f'(x) = (4x) / [ 4(x^4) + 1] \n")

    # arctan(2*(x^2))
    f = lambda t: np.arctan(2*(t**2))

    # 4x/( 4(x^4) +1 )
    df = lambda t: (4*t)/( 4*(t**4) + 1)

    # [4(-12(x^4) +1)]/( 4(x^4) +1 )
    ddf = lambda t: ( 4*( -12*(t**4) +1 ) )/( 4*(t**4) +1 )


    # set data values
    fvals = f(zvals)

    # allocate storage for lagrange interpolant
    pl = np.zeros(zvals.shape, dtype=float)

    # allocate storage for hermite interpolant (evaluated using modified divided diffrence method)
    phdd = np.zeros(zvals.shape, dtype=float)

     # allocate storage for cubic spline interpolant
    cs = np.zeros(zvals.shape, dtype=float)

    # run test for varying numbers of interpolation nodes
    for n in nvals:

        # evaluation at the nodes
        lagrange_nodes = np.linspace(-L, L, n+1)     # lagrange interpolation nodes 
        lagrange_f_yvals = f(lagrange_nodes)         # function values

        # evaluation at the nodes
        hermite_nodes = np.linspace(-L, L, int((n+1)//2)) 
        hermite_f_yvals = f(hermite_nodes)
        hermite_df_yvals = df(hermite_nodes)

        # evaluation at the nodes
        #cubic_spline_nodes = np.linspace(-L, L, n+1)
        cubic_spline_nodes = np.linspace(-L, L, n)
        alpha = df(zvals[0])
        beta = ddf(zvals[-1])
        z_coefs = cubic_spline_coefficients(cubic_spline_nodes,f(cubic_spline_nodes),alpha,beta)
        cubic_spline_f_yvals = []
        for t in cubic_spline_nodes:
            yhat = cubic_spline_evaluate(cubic_spline_nodes,f(cubic_spline_nodes),z_coefs,t)
            cubic_spline_f_yvals.append(yhat)

        # evaluate interpolant lagrange
        pl = lagrange(lagrange_nodes, lagrange_f_yvals, zvals)

        # evaluate interpolant hermite using modified divided diffrence method
        phdd = hermite(hermite_nodes, hermite_f_yvals, hermite_df_yvals, zvals, method='modified_divided_difference')

        # evaluate cubic spline 
        cs = []
        for t in zvals:
            yhat = cubic_spline_evaluate(cubic_spline_nodes,f(cubic_spline_nodes),z_coefs,t)
            cs.append(yhat)


        # generate comparison plots
        fig, axarr = plt.subplots(1,2)

        #
        # overlayed interpolant (normal) plots
        #

        ftitle = '$f(z)$'

        ltitle = '$p_L(z)$'
        htitle = '$p_H(z)$'
        stitle = '$S(z)$'
        
        fplot = axarr[0].plot(zvals, fvals, 'r-', label=ftitle)

        lplot = axarr[0].plot(zvals, pl, 'g--', label=ltitle)
        hplot = axarr[0].plot(zvals, phdd, 'b--', label=htitle)
        splot = axarr[0].plot(zvals, cs, 'k--', label=stitle)

        axarr[0].scatter(lagrange_nodes,lagrange_f_yvals,s=10,color='g', zorder=2)
        axarr[0].scatter(hermite_nodes,hermite_f_yvals,s=10,color='b', zorder=2)
        axarr[0].scatter(cubic_spline_nodes,cubic_spline_f_yvals,s=10,color='k', zorder=2)

        axarr[0].set_xlabel('x')
        axarr[0].set_ylabel('y')
        axarr[0].set_title('Problem 1: $f$ vs $p_{' + str(n) + '}$')

        leg = axarr[0].legend(loc='upper center', shadow=True)

        #
        # overlayed semilogy plots
        #

        l_error_title = '$ | f(z) - p_L(z) | $'
        h_error_title = '$ | f(z) - p_H(z) | $'
        s_error_title = '$ | f(z) - S(z) | $'

        lplot_error = axarr[1].plot(zvals, abs(fvals-pl), 'g--', label=l_error_title)
        hplot_error = axarr[1].plot(zvals, abs(fvals-phdd), 'b--', label=h_error_title)
        splot_error = axarr[1].plot(zvals, abs(fvals-cs), 'k--', label=s_error_title)

        axarr[1].set_xlabel('x')
        axarr[1].set_ylabel('y')
        axarr[1].set_title('Problem 1: $E_{' + str(n) + '}$')

        leg = axarr[1].legend(loc='upper center', shadow=True)


    # display plots
    plt.show()
