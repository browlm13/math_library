#!/usr/bin/env python3

"""
File name: graph_compare.py

    Python Version: 3.6

        Comparing interpolation methods

LJ Brown
Fall 2018
"""

__filename__ = "graph_compare.py"
__author__ = "L.J. Brown"

# imports

# internal library
import logging

# my library
from LagrangeInterpolatingPolynomial import *
from HermiteInterpolatingPolynomial import *
from CubicSpline import *

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




    #
    # test problems
    #

    #
    # test interpolants
    #

    get_t = lambda domain, n: np.linspace(domain[0], domain[1], num=n)

    def create_lagrange_interpolating_polynomial(f, df, ddf, domain, degree):

        # generate nodes
        t = get_t(domain, degree)

        # construct hermite interpolating polynomial class passing nodes aand first derivatives
        lp = LagrangeInterpolatingPolynomial(t, f(t))

        return lp

    def create_hermite_interpolating_polynomial(f, df, ddf, domain, degree):

        # generate nodes
        t = get_t(domain, int((degree+1)//2))

        # construct hermite interpolating polynomial class passing nodes aand first derivatives
        hp = HermiteInterpolatingPolynomial(t, f(t), df(t))

        return hp

    def create_cubic_spline(f, df, ddf, domain, degree):

        # generate nodes
        t = get_t(domain, degree)

        # boundry conditions from HW problem 3
        alpha = df(t[0])
        beta = ddf(t[-1])

        # construct cubic spline class passing nodes and boundry conditions
        cs = CubicSpline(t, f(t), start_bc=(1,alpha), end_bc=(2,beta))

        return cs

    #
    # general parameters
    #

    domain = (-3,3)
    num_evaluations = 800
    xs = np.linspace(domain[0], domain[1], num_evaluations)     # evaluation points

    nvals = [ 41, 21, 11, 5]            # polynomial degrees


    #
    # test functions
    #

    class Function: pass

    # test function 1
    f1 = Function()
    f1.latex_title = '$ atan (2 x^{2}) $'
    f1.f = lambda t: np.arctan(2*(t**2)) # arctan(2*(x^2)),
    f1.df = lambda t: (4*t)/( 4*(t**4) + 1) # 4x/( 4(x^4) +1 )
    f1.ddf = lambda t: ( 4*( -12*(t**4) +1 ) )/( 4*(t**4) +1 ) # [4(-12(x^4) +1)]/( 4(x^4) +1 )

    f2 = Function()
    f2.latex_title = '$ sin(x) $'
    f2.f = lambda t: np.sin(t) 
    f2.df = lambda t: np.cos(t)
    f2.ddf = lambda t: -np.sin(t)

    test_functions = [f1, f2]
    #logger.info("Running test problem 1, f(x) = arctan(2x^2), f'(x) = (4x) / [ 4(x^4) + 1] \n")

    test_interpolants = [create_lagrange_interpolating_polynomial, create_hermite_interpolating_polynomial, create_cubic_spline]

    for n in reversed(nvals):

        # generate comparison plots
        num_test_functions = len(test_functions)
        fig, axarr = plt.subplots(num_test_functions,2)

        






    """

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
    """
