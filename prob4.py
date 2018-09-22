#!/usr/bin/env python3

"""
	File name: prob4.py
	Python Version: 3.6

		Fixed point iteration solver.

		Comparisons of aitken acceleration and non aitken acceleration methods for solving x*.

		Method for non aitken acceleration fixed point iteration:

			1.) Convert function g(x) to root finding form. 
				Where f(x*) = 0 is the root finding form of the fixed point problem function g(x*) = x(*).

									f(x*) = g(x*) - x* = 0.

			2.) Compute the coefficients of the newton interpolating polynomial of f(x)
				so the derivative, f'(x), is easy to compute.

	
	L.J. Brown
	Math5315 @ SMU
	Fall 2018
"""

__filename__ = "prob4.py"
__author__ = "L.J. Brown"

# internal libraries
import logging

# external libraries
import numpy as np

# mylib libraries
from fixed_point import *
from aitken_fp import *

# initilize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# testing
if "__main__" in __name__:

	SHOW_OUTPUT = True

	# log SHOW_OUTPUT value
	logger.info("\n\nSHOW_OUTPUT set to %s.\n\n" % SHOW_OUTPUT)

	# test functions:

	# use: 
	#    Gfun_a = test_functions['Gfun_a']
	#    Ffun_a = root_finding_conversion(Gfun_a)

	fixed_point_functions = {

		'Gfun_a' : lambda x: (x**2)/4 -x/2 -1,
		'Gfun_b' : lambda x: np.cos(x),
		'Gfun_c' : lambda x: (x/3 +x)/2,
		'Gfun_d' : lambda x: np.cosh(x)/x -np.arctan(x)
	}

	# intitial guesses:

	initial_guesses = {
		'Gfun_a' : 2,
		'Gfun_b' : 2,
		'Gfun_c' : 2,
		'Gfun_d' : 2
	}

	# For all problems use an absolute solution tolerance of 10−5, 
	# a relative solution tolerance of 10−10, 
	# allow a maximum of 100 iterations.
	maxit = 100
	atol = 10**(-5)
	rtol = 10**(-10)

	# interp_representation_error
	interp_representation_error = lambda Gfun, x_final: abs(Gfun(x_final) - x_final)

	# run trials
	for Gfun_name, Gfun in fixed_point_functions.items():

		Gfun = fixed_point_functions[Gfun_name]
		x = initial_guesses[Gfun_name]

		# find x* using aiken acceleration and method 1

		# method 1 - in 'fixed_point.py'
		x_target = fixed_point(Gfun, x, maxit, rtol, atol, output=SHOW_OUTPUT)

		# display interp_representation_error
		error1 = interp_representation_error(Gfun, x)
		logger.info("\n\n[Method 1]: \nFinal fixed point target error WITHOUT aitken acceleration \nusing newtons method on a newton interpolating polynomial in place of Ffun.\n \
			\n\tFunction name: %s, \n\n\t\t |Gfun(x) - x| = %s.\n" % (Gfun_name, error1))


		# aiken acceleration
		x_target = aitken_fp(Gfun, x, maxit, rtol, atol, output=SHOW_OUTPUT)
	
		# display interp_representation_error
		error2 = interp_representation_error(Gfun, x)
		logger.info("\n\n[Method 2]: \nFinal fixed point target error using acceleration acceleration.\n \
			\n\tFunction name: %s, \n\n\t\t |Gfun(x) - x| = %s.\n" % (Gfun_name, error2))

