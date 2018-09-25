#!/usr/bin/env python3

"""
	File name: fixed_point.py
	Python Version: 3.6

		Fixed point iteration solver.

	L.J. Brown
	Math5315 @ SMU
	Fall 2018
"""

__filename__ = "fixed_point.py"
__author__ = "L.J. Brown"

# internal libraries
import logging

# external libraries
import numpy as np

# initilize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def fixed_point(Gfun, x, maxit, rtol, atol, output=True):

	prev_x = None
	for i in range(maxit):
		next_x = Gfun(x)

		if output:
			output_string = " n = %s, \tx_{n} = %g," % (i, next_x)
			if prev_x is not None:
				# |e_{n}|/|e_{n-1}| ~= |x_{n}-x_{n-1}|/|x_{n-1}-x_{n-2}|
				approximate_error_string = " \t|e_{n}|/|e_{n-1}| ~= %s"
				approximate_error = abs(next_x - x)/abs(x-prev_x)
				output_string += approximate_error_string % approximate_error
			logger.info(output_string)

		if (abs(x - next_x) < atol + rtol*next_x): 
			return next_x

		prev_x = x
		x = next_x

	# log failure to converge
	logger.info("\n\nFailure to converge in fixed_point.\n\n")

	return x
	