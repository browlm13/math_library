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
    for n in reversed(nvals):

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

        #leg = axarr[0].legend(loc='upper center', shadow=True)
        leg = axarr[0].legend(loc='lower left', shadow=True)

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


""""
problems:

Runge phenominon - especially for lagrange

Flat regions for lagrange and hermite, but not spline.
I think it becomes difficult when the degree of the polynomial is high, this is not true for cubic spline interpolation

   
If there are a large number of data points use piecewise polynomial interpolation
Use piecewise polynomials instead of using a single polynomial (of high degree) to interpolate these points

A spline of degree m is a piecewise polynomial (of degree m) with the maximum possible 
smoothness at each of the points where the polynomials join

optimality theorem?

Polynomials are not well suited to interpolation problems because of Runge's phenomenon -- 
the short version is that if you try to interpolate values using a polynomial, 
you get unexpected oscillations in the interpolation, and using higher degree polynomials can make things worse: 

Say you have the red function above, and you try to interpolate using a fifth order polynomial 
(the blue line) -- that looks pretty rough. So you say, hmmm, let me try a 9th order polynomial. 
And what you find is that it actually makes things worse -- there are weird oscillations between the 
interpolation points. 

A second reason when used in modelling is that cubic splines have local control. 
When you modify a control point for a spline, it only affects the area around the control points 
(for a cubic spline, it is limited by the convex hull of the two prior and two subsequent control points). 
Polynomials: you change one interpolation point and you can change the curve anywhere along the curve 
(e.g. in the above picture, you take the minimum at x=-0.95 and change it, the local minimum at x=0.95 
    will also shift.

for a lot of data sets it can be beneficial. The problem with having lots of data, 
especially if it’s roughly equally spaced apart, is that polynomial interpolation suffers from Runge’s Phenomena. 
This means that as you add more data, the derivatives at each of the data points tend to grow. 
This results in large oscillations between data points that typically don’t tend to be “right”.

Spline interpolation tends to greatly reduce oscillation between data points.
 Part of this is in the derivation of the splines. For typical 3rd order splines, 
 you can actually derive the spline equations using the 1D Beam Equations. Given that, 
 you can use some physical intuition to understand why the splines may have less error between data. 
 You can view the spline result as trying to bend a beam to touch each of the data points. 
 Since it’s a beam, it will deflect a minimal amount between the data points, doing just enough to 
 get to the next data point.

 Thus, spline interpolation often does not have as much error as polynomial interpolation. 
 However, it is possible to space data out in such a way that the polynomial fit will actually perform well.
  But splines do a pretty good job for a generic data set (assuming no outliers), so it’s a bit more reliable.

  One can use, Lagrange's Interpolation formula when there is an uneven spread in the domain set.
Example: When sets such as {0,1,4,8,11,12,17} or {1,2,4,6,9,10} are the domains of functions.
"""
