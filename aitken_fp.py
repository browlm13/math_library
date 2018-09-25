#!/usr/bin/env python3

"""
	File name: aitken_fp.py
	Python Version: 3.6

		Fixed point iteration solver using aitken acceleration.

	L.J. Brown
	Math5315 @ SMU
	Fall 2018
"""

__filename__ = "aitken_fp.py"
__author__ = "L.J. Brown"

# internal libraries
import logging

# external libraries
import numpy as np

# mylib libraries
from newton_interp import *
from newton2 import *

# initilize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def aitken_fp(Gfun, x, maxit, rtol, atol, output=True):

	# for error convergence output
	prev_x = None
	cur_x = x
	next_x = None

	for i in range(maxit):

		# for error convergence output
		cur_x = x

		x0 = x
		x1 = Gfun(x0)
		x2 = Gfun(x1)

		denominator = (x2 - 2*x1 + x0)

		if denominator != 0:
			r = (x2*x0 - x1*x1)/denominator
			s = r-x
			x = r

			# for error convergence output
			next_x = x

		else:
			s = abs(x0 - x)
			x = x2

			# for error convergence output
			next_x = x

		if output:
			output_string = " n = %s, \tx_{n} = %g," % (i, next_x)
			if prev_x is not None:
				# |e_{n}|/|e_{n-1}| ~= |x_{n}-x_{n-1}|/|x_{n-1}-x_{n-2}|
				approximate_error_string = " \t|e_{n}|/|e_{n-1}| ~= %s"
				approximate_error = abs(next_x - cur_x)/abs(cur_x-prev_x)
				output_string += approximate_error_string % approximate_error
			logger.info(output_string)

			# for error convergence output
			prev_x = x0


		if (abs(s) < atol + rtol*s):
			return x

	# log failure to converge
	logger.info("\n\nFailure to converge in aitken_fp.\n\n")

	return x