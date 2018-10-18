#!/usr/bin/env python3

"""
File name: HermiteInterpolatingPolynomial.py

    Python Version: 3.6

            Hermite Interpolation

        Construct Hermite polynomial interpolant, ph, for the input dataset. 
        The parameters required are x values, and a function and its first derivative evaluated at those x values,

            x  = [x0,     x1,     ..., xn    ] 
            y  = [f(x0),  f(x1),  ..., f(xn) ] 
            dy = [f'(x0), f'(x1), ..., f'(xn)]

        methods construct ph where,

            ph(xi)  = f(xi),
            ph'(xi) = f'(xi)  for i = 0, 1, ..., n

        Methods implimented to construct the hermite polynomial are:

            1. A modified Divided Difference basis

            [TODO]: 2. A modified Lagrange basis
            [TODO]: 3. Method of undetermined coefficients

        Method details:

            Method 1. Modified Divided Difference:

                Given n+1 data points (for each xs, ys and dys) construct
                    a hermite polynomial interpolant of degree 2n + 1,

                Modified newton form: (degree 2n + 1)

                    ph(x) = b0 + b1(x - x0) + b2(x - x0)^2 + b3(x - x0)^2(x - x1) 
                        + ... + b2n+1(x - x0)^2(x - x1)^2...(x - xn-1)^2(x - xn)

                    where,

                        bi = f[x0,...,xi] ~ (2n+2 coefficients)

                    and where,

                        f[xi,xj,...,xk,xl] =  f[xi,...,xk] - f[xj,...,xl]
                                                -------------------
                                                      xl - xi

                    or in case of even indices in 3rd column (for 0/0 case),

                        f'(xi) = f[xi,xi+1]
        

LJ Brown
Fall 2018
"""

__filename__ = "HermiteInterpolatingPolynomial.py"
__author__ = "L.J. Brown"

# internal libraries
import logging
from functools import reduce

# external libraries
import numpy as np

# initilize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class HermiteInterpolatingPolynomial:

    """
        use example:

            # test function and its derivatives:
            f = lambda t: np.arctan(2*(t**2))                       # arctan(2*(x^2))
            df = lambda t: (4*t)/( 4*(t**4) + 1)                    # 4x/( 4(x^4) +1 )

            # general parameters
            L = 3.0 # range (-L, L) for nodes
            n = 41  # number of nodes

            # generate nodes
            t = np.linspace(-L, L, int((n+1)//2)) # for hermite nodes to have degree n
            y = f(t)

            # construct cubic spline class passing nodes and first derivatives
            hp = HermiteInterpolatingPolynomial(t, y, dy)

            # evaluation points to test cubic spline after construction
            xs = np.linspace(-L, L, 800)
        
            # evaluate
            yhats = hp(xs)

    """

    def __init__(self, t, y, dy):

        self.t = t
        self.y = y 
        self.dy = dy
  
        # construct hermite polynoial using modified divided diffrences method
        self.hp = construct_hermite_interpolating_polynomial_using_modified_divided_differences(self.t, self.y, self.dy)

    def __call__(self, xs):

        yhats = self.hp(xs)
        return yhats



def hermite_coeffs(xs, ys, dys, output=True):
    """
        Given n+1 data points (for each xs, ys and dys) construct 
        a hermite polynomial interpolant of degree 2n + 1,

            Newton form: (degree n)

                p(x) = b0 + b1(x - x0) + b2(x - x0)(x - x1) + ... + bn+1(x - x0)...(x - xn)

            Modified newton form: (degree 2n + 1)

                ph(x) = b0 + b1(x - x0) + b2(x - x0)^2 + b3(x - x0)^2(x - x1) 
                    + ... + b2n+1(x - x0)^2(x - x1)^2...(x - xn-1)^2(x - xn)

        Use Newtons Divded Diffrence Method, and return the coefficients, 

            bs = b0, b1, ..., b2n+1

        where,

            bi = f[x0,...,xi]

        -----------
        P3 example: (given 2 datapoints)
        -----------

            -------
            |xs|ys|
            -------     (b0)
                        /
            x0  y0 = f[x0]         (b1)
                                  /
                           dy_0 = f[x0,x1]        (b2)
                                                 /
            x0  y0 = f[x1]               f[x0,x1,x2]              (b3)
                                                                /
                            f[x1,x2]                 f[x0,x1,x2,x3]

            x1  y1 = f[x2]               f[x1,x2,x3]

                           dy_1 = f[x2,x3]

            x1  y1 = f[x3]


        Code Structure:

            M ~ nxn Upper triangular matrix to store 
                intermediate divided diffrence values

            where,

                M[i,j] = f[xi,...,xj]


                        | f[x0]     f[x0,y1]       f[x0,x1,x2]     f[x0,x1,x2,x3] | 0
                        |           f[x1]          f[x1,x2]        f[x1,x2,x3]    | 1
                M   =   |                          f[x2]           f[x2,x3]       | 2
                        |                                          f[x3]          | 3
                           0          1               2               3


            and where,

                bs = M[0,:]

            and where,

                f[xi,xj,...,xk,xl] =  f[xi,...,xk] - f[xj,...,xl]
                                        -------------------
                                              xl - xi

                f[xi,xj,...,xk,xl] = ( M[i,k] - M[j,l] )/( xs[l] - xs[i] )

            or in case of even indices in 3rd column,

                dy_i = f[xi,xi+1]

        :param xs: numpy array of length n+1 containg [x0,     x1,     ..., xn    ].
        :param ys: numpy array of length n+1 containg [f(x0),  f(x1),  ..., f(xn) ].
        :param dys: numpy array of length n+1 containg [f'(x0), f'(x1), ..., f'(xn)].
        :param output: display output boolean, default set to True.
        :returns: numpy array of length 2n+2 containing modified divided difference coefficients [b0, b1, ..., b2n+1].

    """

    # modify ys and xs so that they repeate values: [x0,x1,x2] -> [x0,x0,x1,x1,x2,x2]
    mod_xs = np.repeat(xs, 2)
    mod_ys = np.repeat(ys, 2)

    n = mod_ys.shape[0]

    # initilize Matrix to zeros
    M = np.zeros(shape=(n,n))

    # get row, col indices of diagnol
    row, col = np.diag_indices_from(M)

    # set main diagnol to y values
    M[row,col] = mod_ys

    # compute succsessive diagnoal values
    for diag in range(1,n):
        # r,c ~ row, col to store calculated divided diffrence in
        for r in range(0,n-diag):
            c = diag+r
            low, high = r,c

            # at even rows during the first diagnol use y' values
            if (r%2 ==0) and (diag == 1):
                M[r,c] = dys[r//2]
            else:
                # calculate normal divided diffrence
                M[r,c] = divided_diffrence(low,high,M,mod_xs)

    # grab coeffiecients from first row in difference table
    bs = M[0,:]

    # display matrix
    if output:
        logger.info("\n\n\tDivided Difference Table:\n\n%s\n" % M)
        logger.info("\n\n\tComputed Coeffcients:\n\n(\'bs = b0, b1, ..., bn-1\')\n%s\n" % bs)

    # return coeffcients
    return bs

def divided_diffrence(low, high, M, xs):
    """
        [low,...,high] = ( M[ low, high-1 ] - M[ low+1, high ] )/( xs[high] - xs[low] )

        return [low,...,high]
    """
    xs_flipped = np.flip(xs, axis=0)
    return (M[ low, high-1 ] - M[ low+1, high ])/( xs_flipped[high] - xs_flipped[low] )

def evaluate_newton_polynomial(x, xs,bs, modify=False):
    """
        Evaluates the polynomial at x given its nodes, xs, and its coefficients, bs.

        Unmodified form:
        
            Pn(x) = b0 + b1(x - x0) + b2(x - x0)(x - x1) + ... + bn-1(x - x0)...(x - xn)

        Modified form:
            Phn(x) = b0 + b1(x - x0) + b2(x - x0)^2 + b3(x - x0)^2(x - x1) 
                    + ... + b2n+1(x - x0)^2(x - x1)^2...(x - xn-1)^2(x - xn)

    """

    if modify:

        # modify xs: [x0, x0, x1, x1, ..., xn-1, xn-1, xn]
        xs = np.repeat(xs, 2)[:-1]

    #
    #          mimic: qk(x) = (1)(x - x0)(x - x1)...(x - xk)
    # or if modified: qk(x) = (1)(x - x0)^2(x - x1)^2...(x - xk)

    # [(1),(x - x0),(x - x1),...,(x - xn)]
    multipliers = [1] + [x - xi for xi in xs]   

    # q(k) = (1)(x - x0)(x - x1)...(x - xk)
    q = lambda k: reduce(lambda x, y: x*y, multipliers[:k+1])

    # compute and return pn(x)
    pnx = 0
    for i in range(len(bs)):
        pnx += bs[i]*q(i)

    return pnx

def construct_hermite_interpolating_polynomial_using_modified_divided_differences(xs, ys, dys):
    """
        Modified Divided Diffrences:

        Modified newton basis:

            Phn(x) = b0 + b1(x - x0) + b2(x - x0)^2 + b3(x - x0)^2(x - x1) 
                        + ... + b2n+1(x - x0)^2(x - x1)^2...(x - xn-1)^2(x - xn)
        use: 
            ph = construct_hermite_interpolating_polynomial_using_divided_differences(xs, ys, dys)
            approximation = ph(x)

        returns a function for the constructed hermite polynomial interpolant.
    """

    bs = hermite_coeffs(xs,ys,dys, output=False)
    ph = lambda x: evaluate_newton_polynomial(x, xs, bs, modify=True)

    return ph

def hermite(x, y, dy, z):
    """ 
        Usage: p = hermite(x,y,z)

        This routine evaluates the Hermite interpolating                               
        polynomial p(z) at all points in z and returns a numpy array of the results.

        :param x: numpy array of length n+1 containg [x0,     x1,     ..., xn    ].
        :param y: numpy array of length n+1 containg [f(x0),  f(x1),  ..., f(xn) ].
        :param dy: numpy array of length n+1 containg [f'(x0), f'(x1), ..., f'(xn)].
        :param z: numpy array of evaluation points (x values).
        :param method: method to use when constructing the hermite polynomial interpolant, deafult 'modified_divided_difference', other option 'modified_lagrange'. 
        :returns: numpy array of same length as z containing approximate values.
    """

    # check inputs
    if (x.size != y.size != dy.size):
        raise ValueError("lagrange error: (x,y) have different sizes")

    # construct hermite polynomial using modified divided differences
    ph = construct_hermite_interpolating_polynomial_using_modified_divided_differences(x, y, dy)


    return ph(z)


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
    df = lambda t: (4*t)/( 4*(t**4) + 1)                    # 4x/( 4(x^4) +1 )

    # general parameters
    L = 3.0 # range (-L, L) for nodes
    nvals = [ 41, 21, 11, 5] # number of nodes

    for n in nvals:

        # generate nodes
        t = np.linspace(-L, L, int((n+1)//2)) # for hermite nodes to have degree n

        # construct hermite interpolating polynomial class passing nodes aand first derivatives
        hp = HermiteInterpolatingPolynomial(t, f(t), df(t))

        # evaluation points to test cubic spline after construction
        xs = np.linspace(-L, L, 400)

        #
        # plot results
        #

        # generate comparison plots
        fig, axarr = plt.subplots(1,2)

        ftitle = '$f(z)$'
        htitle = '$H_{p}(z)$'

        # plot evaluations at nodes using scatter
        axarr[0].scatter(t,f(t),s=10,color='r', zorder=2)
        axarr[0].scatter(t,hp(t),s=5,color='b', zorder=2)

        # plot function graphs
        fplot = axarr[0].plot(xs, f(xs), 'r-', label=ftitle)
        hplot = axarr[0].plot(xs, hp(xs) , 'b--', label=htitle)

        axarr[0].set_xlabel('x')
        axarr[0].set_ylabel('y')

        h_error_title = '$ | f(z) - H_{p}(z) | $'

        splot_error = axarr[1].plot(xs, abs(f(xs)-hp(xs)), 'b--', label=h_error_title)

        axarr[1].set_xlabel('x')
        axarr[1].set_ylabel('y')
        axarr[1].set_title('Problem 1: $E_{' + str(n) + '}$')


        leg = axarr[0].legend(loc='upper center', shadow=True)


    plt.show()

