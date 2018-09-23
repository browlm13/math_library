#!/usr/bin/env python3

"""
	File name: prob4_analysis.py
	Python Version: 3.6

		Fixed point iteration solver.

		Comparisons of aitken acceleration and non aitken acceleration methods for solving x*.

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

def get_iteration_display_points(Gfun, x0, num_iters):

	# first iteration: xaxis (x0, 0) -> y=g (x0, Gfun(x0)) -> y=x (Gfun(x0), Gfun(x0)) -> 
	xs = [x0]
	ys = [0]

	guesses = [x0]

	for i in range(num_iters):
		#y=g (yk, Gfun(xk))
		xk0 = xs[-1]
		yk0 = Gfun(xk0)

		#y=x (Gfun(xk), Gfun(xk))
		xk1 = Gfun(xk0)
		yk1 = Gfun(xk0)

		#y=0, x=xk+1 (Gfun(xk), 0)
		xk2 = Gfun(xk0)
		yk2 = 0

		xs += [xk0, xk1, xk2]
		ys += [yk0, yk1, yk2]

		guesses += [xk2]

	return xs, ys, guesses


def analysis_plot(fixed_point_functions, order, domain, num_fixed_point_iterations=8, num_points=300):

	import matplotlib.pylab as plt
	import matplotlib.colors as colors
	import matplotlib.cm as cmx

	# setting up figure
	num_plots = len(fixed_point_functions)

	fig, axs = plt.subplots(1, num_plots, figsize=(15, 7), facecolor='w', edgecolor='k')
	fig.subplots_adjust(hspace = .25, wspace=.7)
	axs = axs.ravel()

	for i, name in enumerate(order):

		Gfun = fixed_point_functions[name]['method']
		latex = fixed_point_functions[name]['latex']
		x0 = fixed_point_functions[name]['initial guess']

		x = np.linspace(domain[0],domain[1],num_points)
		y = x # compute data points for plot y = x
		gx = Gfun(x) # compute y data points for plot y = g(x)
		

		# create plot for this function
		axs[i].plot( x, gx, 'r', label=latex)     # y = g(x) in red
		axs[i].plot( x, y, 'b', label='x' ) 	   # y = x in blue
		axs[i].set_title(name, fontsize=12)

		axs[i].grid(True, which='both')

		axs[i].axhline(y=0, color='k')
		axs[i].axvline(x=0, color='k')

		axs[i].set_aspect('equal') # equal
		axs[i].set_ylim([domain[0],domain[1]]) # set y lim

		# display fixed point finding iterations
		xs, ys, guesses = get_iteration_display_points(Gfun, x0, num_fixed_point_iterations)
		axs[i].plot(xs, ys, color='m', marker=None, label='fixed point iteration', linestyle='dashed', linewidth=0.5)
	
		cm = plt.get_cmap("hot") 
		cNorm = colors.Normalize(vmin=0, vmax=len(guesses))
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
		for idx,g in enumerate(guesses):
			colorVal = scalarMap.to_rgba(idx)
			alpha=1 #1/(idx+1)
			markersize= 20/(idx+1)
			axs[i].plot(g,0, marker='x', color=colorVal, alpha=alpha,markersize=markersize)

		cNorm = colors.Normalize(vmin=0, vmax=len(xs))
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
		for idx in range(len(xs)):
			colorVal = scalarMap.to_rgba(idx)
			alpha=0.75 #1/(idx+1)
			markersize= 10/(idx+1)
			axs[i].plot(xs[idx],ys[idx], marker='.', color=colorVal, alpha=alpha, markersize=markersize)

		# show legned
		axs[i].legend(loc=9, bbox_to_anchor=(0.5, -0.2))

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

	# plot fixed point iterations
	order = ['a', 'b', 'c', 'd']
	domain = (-2,2)
	num_fixed_point_iterations = 10
	analysis_plot(fixed_point_functions, order, domain, num_fixed_point_iterations)


	# For all problems use an absolute solution tolerance of 10−5, 
	# a relative solution tolerance of 10−10, 
	# allow a maximum of 100 iterations.
	maxit = 100
	atol = 10**(-5)
	rtol = 10**(-10)

	# interp_representation_error
	#interp_representation_error = lambda Gfun, x_final: abs(Gfun(x_final) - x_final)

	"""
	# run trials
	for Gfun_name, Gfun in fixed_point_functions.items():

		Gfun = fixed_point_functions[Gfun_name]
		x = initial_guesses[Gfun_name]
	"""

