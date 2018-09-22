#!/usr/bin/env python3

"""
    File name: bisection.py
    Python Version: 3.6

    L.J. Brown
    Math5315 @ SMU
    Fall 2018
"""

__filename__ = "bisection.py"
__author__ = "L.J. Brown"

# external libraries
import numpy as np

def bisect(Ffun, miter, a, b, tol):
    """
    	Bisection root finding method.
    	Where f(a) and f(b) have opposite signs and a < b.
    """
    opposite_signs_y = lambda a, b: (Ffun(a)*Ffun(b) < 0)
    assert opposite_signs_y(a,b) and (a < b)

    for i in range(miter):
        c = (a + b) / 2.
        if opposite_signs_y(a,b): b = c
        else: a = c

        if Ffun(c) <= tol:
        	return c

    return c