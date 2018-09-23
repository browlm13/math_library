#!/usr/bin/env python3

"""
	File name: prob4_analysis.py
	Python Version: 3.6

		Fixed point iteration solver.

		Comparisons of aitken acceleration and non aitken acceleration methods for solving x*.

		a) g(x) = (1/4)x^2 + (-1/2)x -1,  x0 = 2

		plot analysis of y = g(x) vesus y = x
	
	L.J. Brown
	Math5315 @ SMU
	Fall 2018
"""

__filename__ = "prob4_analysis.py"
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


def analysis_plot(fixed_point_functions, order, domain, num_points=300):

	import matplotlib.pylab as plt

	# setting up figure
	num_plots = len(fixed_point_functions)

	fig, axs = plt.subplots(1, num_plots, figsize=(15, 6), facecolor='w', edgecolor='k')
	fig.subplots_adjust(hspace = .75, wspace=.1)
	axs = axs.ravel()

	for i, name in enumerate(order):

		Gfun = fixed_point_functions[name]['method']
		latex = fixed_point_functions[name]['latex']

		x = np.linspace(domain[0],domain[1],num_points)
		y = x # compute data points for plot y = x
		gx = Gfun(x) # compute y data points for plot y = g(x)
		

		# create plot for this function
		axs[i].plot( x, gx, 'r', label=latex)     # y = g(x) in red
		axs[i].plot( x, y, 'b', label='x' ) 	   # y = x in blue
		axs[i].set_title(name, fontsize=12)

		#axs[i].set_aspect('equal')
		axs[i].grid(True, which='both')

		axs[i].axhline(y=0, color='k')
		axs[i].axvline(x=0, color='k')

		# show legned
		axs[i].legend()

	plt.show()


# testing
if "__main__" in __name__:

	SHOW_OUTPUT = True

	# log SHOW_OUTPUT value
	logger.info("\n\nSHOW_OUTPUT set to %s.\n\n" % SHOW_OUTPUT)

	# test functions:

	fixed_point_functions = {
		'a' : {
				'method' : lambda x: (x**2)/4 -x/2 -1,
				'initial guess' : 2,
				'latex' : r"$\frac{x^{2}}{4} -\frac{x}{2} -1$"
				},
		'b' : {
				'method' : lambda x:  np.cos(x),
				'initial guess' : 2,
				'latex' : r"$cos(x)$"
				},

		'c' : {
				'method' : lambda x: (x/3 +x)/2,
				'initial guess' : 2,
				'latex' : r"$\frac{1}{2}\left(\frac{3}{x} + x \right)$"
				},

		'd' : {
				'method' : lambda x: np.cosh(x)/x -np.arctan(x),
				'initial guess' : 2,
				'latex' : r"$\frac{cosh(x)}{x} -arctan(x)$"
				}
	}

	order = ['a', 'b', 'c', 'd']

	# For all problems use an absolute solution tolerance of 10−5, 
	# a relative solution tolerance of 10−10, 
	# allow a maximum of 100 iterations.
	maxit = 100
	atol = 10**(-5)
	rtol = 10**(-10)

	domain = (-5,5)
	#analysis_plot(fixed_point_functions, domain)
	analysis_plot(fixed_point_functions, order, domain)

	# interp_representation_error
	#interp_representation_error = lambda Gfun, x_final: abs(Gfun(x_final) - x_final)

	"""
	# run trials
	for Gfun_name, Gfun in fixed_point_functions.items():

		Gfun = fixed_point_functions[Gfun_name]
		x = initial_guesses[Gfun_name]
	"""

