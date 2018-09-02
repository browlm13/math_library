#!/usr/bin/env python

"""

    Horner's Method

"""


def poly_eval(x, coeffs):
    """
        evalute polynomial at x using horner's method.
        :param x: type float, evaluate the polynomial at x.
        :param coeffs: type list floats, coeffients of polynomial degree n, starting with lowest at index 0.
        Coeffients for nth degree at coeffs[n].
        :returns: float, value of polynomial evaulated at x.
    """

    n = len(coeffs) -2
    p = float(coeffs[-1])

    for k in range(n,-1,-1):
        p = p*x + coeffs[k]

    return p


import numpy as np

ptest = [3, 3, 3 ,4, 4, 4]
x = 5.0

print(np.polynomial.polynomial.polyval(x, ptest))
print(poly_eval(x, ptest))