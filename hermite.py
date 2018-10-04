#!/usr/bin/env python3

# internal libraries
from functools import reduce

def construct_lagrange_basis_polynomial(i, x_vals):
    """
        Returns a function to evaluate the lagrange basis polynomial li(x).

        li(xi) = 1
        li(xj) = 0 for j != i
        li(x) = other for x not in x_vals

        use:

            ex:
                x_vals = [0,1,2,3]
                l = construct_lagrange_basis_polynomial(2,x_vals)

                # l(2) = 1
                # l(0) = l(1) = l(3) = 0

    """
    xi = x_vals[i]

    # skipping xi
    x_vals = list(x_vals)
    used_x_vals = x_vals[:i] + x_vals[i+1:]

    # for all j != i: (xi - xj)
    demonimators = [xi - xj for xj in used_x_vals]
    denominator_result = reduce(lambda a, b: a*b, demonimators)

    # lagrangian basis polynomial function constructed at index i
    def li(x):

         # for all j != i: (x - xj)
        numerators = [x - xj for xj in used_x_vals]
        numerator_result = reduce(lambda a, b: a*b, numerators)

        return numerator_result/denominator_result

    # return basis polynomial function
    return li

# construct_lagrange_basis_polynomial_derivative
def construct_lagrange_basis_polynomial_derivative(i, x_vals):
    """
        Returns a function to evaluate the derivative of lagrange basis polynomial dli(x).

        dli(x) = li(x)*( sum(1/(x - xi)) for i != j )

    """

    xi = x_vals[i]

    # skipping xi
    x_vals = list(x_vals)
    used_x_vals = x_vals[:i] + x_vals[i+1:]

    # consturct lagrange basis polynomial at i
    li = construct_lagrange_basis_polynomial(i,x_vals)

    # derivative of lagrangian basis polynomial function constructed at index i
    def dli(x):

        # for all j != i: 1/(x - xi)
        multipliers = [1/(x - xi) for xi in used_x_vals]
        multiplier_result = reduce(lambda a, b: a*b, multipliers)

        return li(x)*multiplier_result

    # return derivative of basis polynomial function
    return dli


# construct Ai for hermite interpolation
def construct_Ai(i, x_vals):
    """
        Returns a function to evaluate Ai(x).

        Ai(x) = [1 - 2(x - xi)dli(xi)](li( x )**2)
    """
    xi = x_vals[i]
    li = construct_lagrange_basis_polynomial(i,x_vals)
    dli = construct_lagrange_basis_polynomial_derivative(i, x_vals)

    # evaluate dli at xi
    #dli_xi = dli(xi)

    def Ai(x):
        #return (1 - 2*(x - xi)*dli_xi)*(li(x)**2)
        return (1 - 2*(x - xi)*dli(xi))*(li(x)**2)

    return Ai

# construct Bi for hermite interpolation
def construct_Bi(i, x_vals):
    """
        Returns a function to evaluate Bi(x).

        Bi(x) = (x - xi)*(li(x)**2)
    """
    xi = x_vals[i]
    li = construct_lagrange_basis_polynomial(i,x_vals)

    def Bi(x):
        return (x - xi)*(li(x)**2)

    return Bi

def construct_hermite_interpolating_polynomial(x_values, y_values, dy_values):
    """
        Returns a function to evaluate p(x)

        ph(x) = sum_i2n(yi*Ai(x)) + sum_i2n(dyi*Bi(x))
    """

    def ph(x):

        result = 0
        for i in range(len(x_values)):
            yi, dyi = y_values[i], dy_values[i]
            Ai, Bi = construct_Ai(i, x_values), construct_Bi(i, x_values)

            result += yi*Ai(x) + dyi*Bi(x)

        return result

    return ph

"""
x_vals = [0,1,2,3]
y_vals = [1,5,6,7]
dy_vals = [2,3,4,5]

l2 = construct_lagrange_basis_polynomial(2,x_vals)

A2 = construct_Ai(2, x_vals)

B2 = construct_Bi(2, x_vals)

p = construct_hermite_interpolating_polynomial(x_vals, y_vals, dy_vals)
"""


def hermite(x, y, dy, z):
    """ 
        Usage: p = hermite(x,y,z)

        This routine evaluates the Hermite interpolating                               
        polynomial p(z) 

        [todo]: defined in the Lagrange basis by
        p(z) = y[0]*l1(z) + y[1]*l2(z) + ... + y[n-1]*ln(z)
        Inputs:   x  nodal locations [array of length n],  (we assume that x(i) ~= x(j) for i != j)
                  y  data values [array of length n] 
                  z  evaluation point

        Outputs:  p  p(z)

        LJ Brown
        Math5315 @ SMU
        Fall 2018

        all you need is A and B

    """

    # check inputs
    if (x.size != y.size != dy.size):
        raise ValueError("lagrange error: (x,y) have different sizes")

    ph = construct_hermite_interpolating_polynomial(x, y, dy)

    #return np.array(ph(z))
    return ph(z)


