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

	for i in range(maxit):
		x0 = x
		x1 = Gfun(x0)
		x2 = Gfun(x1)

		denominator = (x2 - 2*x1 + x0)

		if denominator != 0:
			r = (x2*x0 - x1*x1)/denominator
			s = r-x
			x = r
		else:
			x = x2

		if (abs(s) < atol + rtol*s):
			return x

	# log failure to converge
	logger.info("\n\nFailure to converge in aitken_fp.\n\n")

	return x