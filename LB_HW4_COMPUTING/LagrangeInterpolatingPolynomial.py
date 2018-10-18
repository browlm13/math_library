#!/usr/bin/env python3

"""
File name: LangrangeInterpolatingPolynomial.py

    Python Version: 3.6

            Lagrange Interpolation

        Construct a Lagrange interpolating polynomial, p(z), from n datapoints
     
                p(t) = y_0*l1(t) + y_1*l2(t) + ... + y_n-1*ln(z).

        Where li(t) is the lagrange basis polynomial.
    
 Source: Daniel R. Reynolds 

Daniel R. Reynolds 
LJ Brown
Fall 2018
"""

__filename__ = "LangrangeInterpolatingPolynomial.py"
__author__ = ["L.J. Brown", "Daniel R. Reynolds"]

# internal libraries
import logging
from functools import reduce

# external libraries
import numpy as np

# initilize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class LagrangeInterpolatingPolynomial:

    """
        use example:

            # test function and its derivatives:
            f = lambda t: np.arctan(2*(t**2))                       # arctan(2*(x^2))

            # general parameters
            L = 3.0 # range (-L, L) for nodes
            n = 41  # number of nodes

            # generate nodes
            t = np.linspace(-L, L, n) 
            y = f(t)

            # construct cubic spline class passing nodes and first derivatives
            lp = LagrangeInterpolatingPolynomial(t, y)

            # evaluation points to test cubic spline after construction
            xs = np.linspace(-L, L, 800)
        
            # evaluate
            yhats = lp(xs)

    """

    def __init__(self, t, y):

        self.t = t
        self.y = y 

        n = len(self.t)
        self.degree = n

        self.latex_title = self.get_latex_title()
        self.latex_error_title = self.get_latex_error_title()
  
    def __call__(self, xs):

        yhats = lagrange(self.t, self.y, xs)
        return yhats

    @staticmethod
    def num_nodes_given_degree(degree):
        return degree

    @staticmethod
    def get_latex_title():
        title = '$Lagrange_{p}(x)$'
        return title

    @staticmethod
    def get_latex_error_title():
        error_title = '$ | f(x) - Lagrange_{p}(x) | $'
        return error_title



def lagrange(x, y, z):
    """ Usage: p = lagrange(x,y,z)                               """
    """                                                          """
    """ This routine evaluates the Lagrange interpolating        """
    """ polynomial p(z) defined in the Lagrange basis by         """
    """   p(z) = y[0]*l1(z) + y[1]*l2(z) + ... + y[n-1]*ln(z).   """
    """                                                          """
    """ Inputs:   x  nodal locations [array of length n],        """
    """              (we assume that x(i) ~= x(j) for i != j)    """
    """           y  data values [array of length n]             """
    """           z  evaluation point                            """
    """ Outputs:  p  p(z)                                        """
    """                                                          """
    """ Daniel R. Reynolds                                       """
    """ Math5315 @ SMU                                           """
    """ Fall 2018                                                """

    # check inputs
    if (x.size != y.size):
        raise ValueError("lagrange error: (x,y) have different sizes")

    n = x.size

    # initialize output
    p = 0.0

    # iterate over Lagrange basis functions
    for k in range(n):
   
        # initialize l (the kth Lagrange basis function)
        l = 1.0
   
        # iterate over data to construct l(z)
        for j in range(n):
            # exclude the k-th data point
            if (j != k):
                l *= (z-x[j]) / (x[k]-x[j])
   
        # add contribution from this basis function (and data value) into p
        p += y[k]*l
   
    return p

#
#   Testing
#

import matplotlib.pyplot as plt

if __name__ == '__main__':

    #
    # testing
    #

    # test function and its derivatives:
    f = lambda t: np.arctan(2*(t**2))                       # arctan(2*(x^2))

    # general parameters
    L = 3.0 # range (-L, L) for nodes
    nvals = [ 41, 21, 11, 5] # number of nodes

    for n in nvals:

        # generate nodes
        t = np.linspace(-L, L, n+1)

        # construct hermite interpolating polynomial class passing nodes aand first derivatives
        lp = LagrangeInterpolatingPolynomial(t, f(t))

        # evaluation points to test cubic spline after construction
        xs = np.linspace(-L, L, 400)

        #
        # plot results
        #

        # generate comparison plots
        fig, axarr = plt.subplots(1,2)

        ftitle = '$f(z)$'
        ltitle = lp.latex_title() #'$L_{p}(z)$'

        # plot evaluations at nodes using scatter
        axarr[0].scatter(t,f(t),s=10,color='r', zorder=2)
        axarr[0].scatter(t,lp(t),s=5,color='y', zorder=2)

        # plot function graphs
        fplot = axarr[0].plot(xs, f(xs), 'r-', label=ftitle)
        lplot = axarr[0].plot(xs, lp(xs) , 'b--', label=ltitle)

        axarr[0].set_xlabel('x')
        axarr[0].set_ylabel('y')

        l_error_title = lp.latex_error_title()

        lplot_error = axarr[1].plot(xs, abs(f(xs)-lp(xs)), 'b--', label=l_error_title)

        axarr[1].set_xlabel('x')
        axarr[1].set_ylabel('y')
        axarr[1].set_title('Problem 1: $E_{' + str(n) + '}$')


        leg = axarr[0].legend(loc='upper center', shadow=True)
        leg = axarr[1].legend(loc='upper center', shadow=True)


    plt.show()


