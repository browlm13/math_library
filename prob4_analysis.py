#!/usr/bin/env python3

"""
	File name: prob4_analysis.py
	Python Version: 3.6

		Fixed point iteration solver.

		Comparisons of aitken acceleration and non aitken acceleration methods for solving x*.

		plot analysis of y = g(x) vesus y = x

		TODO: add error converge to output during iteration / revamp output in general
				ek+1
				----
				ek	
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
import matplotlib.pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

# mylib libraries
from fixed_point import *
from aitken_fp import *

# initilize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_fp_iteration_display_points(Gfun, x0, atol, rtol, maxits):
	# first iteration: xaxis (x0, 0) -> y=g (x0, Gfun(x0)) -> y=x (Gfun(x0), Gfun(x0)) -> 
	xs = [x0]
	ys = [0]

	guesses = [x0]

	for i in range(maxits):
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

		s = xs[-1] - Gfun(xs[-1])
		if (abs(s) < atol + rtol*s):
			return xs, ys, guesses, i+1, True # return number of iterations and wether or not there was convergence

	return xs, ys, guesses, i+1, False #  return number of iterations and wether or not there was convergence

def get_aa_iteration_display_points(Gfun, x0, atol, rtol, maxits):

	xs = [x0]
	ys = [0]
	guesses = [x0]

	for i in range(maxits):
		x = xs[-1]
		x0 = xs[-1]

		# aitken acceleration
		x1 = Gfun(x0)
		x2 = Gfun(x1)

		denominator = (x2 - 2*x1 + x0)

		if denominator != 0:
			r = (x2*x0 - x1*x1)/denominator
			s = r-x
			x = r
		else:
			s = abs(x0 - x)
			x = x2

		xs += [x, x]
		ys += [0, Gfun(x)]
		guesses += [x]

		if (abs(s) < atol + rtol*s):
			return xs, ys, guesses, i+1, True  # return number of iterations and wether or not there was convergence

	return xs, ys, guesses, i+1, False  # return number of iterations and wether or not there was convergence


def plot_function(Gfun, latex, x0, domain, num_points, axsi, title):

		x = np.linspace(domain[0],domain[1],num_points)
		y = x # compute data points for plot y = x
		gx = Gfun(x) # compute y data points for plot y = g(x)
		

		# create plot for this function
		axsi.plot( x, gx, 'r', label=latex)     # y = g(x) in red
		axsi.plot( x, y, 'b', label='x' ) 	   # y = x in blue

		axsi.set_title(title, fontsize=8)

		axsi.grid(True, which='both')

		axsi.axhline(y=0, color='k')
		axsi.axvline(x=0, color='k')

		axsi.set_aspect('equal') # equal
		axsi.set_ylim([domain[0],domain[1]]) # set y lim

		axsi.legend(loc=9, bbox_to_anchor=(0.5, -0.2))


def plot_iterations(axsi, xs, ys, guesses, its, convergence):

	axsi.plot(xs, ys, color='m', marker=None, linestyle='dashed', linewidth=0.5)

	text = "Convergence (iter:%s)" % its
	text_color = "g"
	if not convergence:
		text = "No Convergence (iter:%s)" % its
		text_color = "r"

	ylim = axsi.get_ylim()
	axsi.text(0, ylim[0],text, color=text_color, fontsize=8, ha='center', va='bottom', alpha=0.8)

	cm = plt.get_cmap("hot") 
	cNorm = colors.Normalize(vmin=0, vmax=len(guesses))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	for idx,g in enumerate(guesses):
		colorVal = scalarMap.to_rgba(idx)
		alpha=1 #1/(idx+1)
		markersize= 20/(idx+1)
		axsi.plot(g,0, marker='x', color=colorVal, alpha=alpha,markersize=markersize)

	cNorm = colors.Normalize(vmin=0, vmax=len(xs))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	for idx in range(len(xs)):
		colorVal = scalarMap.to_rgba(idx)
		alpha=0.75 #1/(idx+1)
		markersize= 10/(idx+1)
		axsi.plot(xs[idx],ys[idx], marker='.', color=colorVal, alpha=alpha, markersize=markersize)


def analysis_plot(fixed_point_functions, order, domain, atol, rtol, maxits, num_points=300):

	# setting up figure
	num_functions = len(fixed_point_functions)

	fig, axs = plt.subplots(2, num_functions, figsize=(15, 7), facecolor='w', edgecolor='k')
	fig.subplots_adjust(hspace = .5, wspace=.7, top=0.94, bottom=0.18)
	axs = axs.ravel()

	for i, name in enumerate(order):

		Gfun = fixed_point_functions[name]['method']
		latex = fixed_point_functions[name]['latex']
		x0 = fixed_point_functions[name]['initial guess']

		# plot non aitken accelerated fixed point iterations
		title = "%s). Without aitken acceleration" % name
		axsi = axs[i]
		plot_function(Gfun, latex, x0, domain, num_points, axsi, title)

		# plot aitken accelerated fixed point itterations
		title = "%s). With aitken acceleration" % name
		axsi = axs[i+num_functions]
		plot_function(Gfun, latex, x0, domain, num_points, axsi, title)


		#
		#	Display iterations without aitken acceleration
		#

		# display fixed point finding iterations
		xs, ys, guesses, its, convergence = get_fp_iteration_display_points(Gfun, x0, atol, rtol, maxits)	
		plot_iterations(axs[i], xs, ys, guesses, its, convergence)


		#
		#	Display iterations with aitken acceleration
		#

		# display fixed point finding iterations
		xs, ys, guesses, its, convergence = get_aa_iteration_display_points(Gfun, x0, atol, rtol, maxits)	
		plot_iterations(axs[i+num_functions], xs, ys, guesses, its, convergence)


	plt.show()


# testing
if "__main__" in __name__:

	SHOW_OUTPUT = True

	# log SHOW_OUTPUT value
	logger.info("\n\nSHOW_OUTPUT set to %s.\n\n" % SHOW_OUTPUT)

	# test functions:
	fixed_point_functions = {
		'A' : {
				'method' : lambda x: (x**2)/4 -x/2 -1,
				'initial guess' : 2,
				'latex' : r"$\frac{x^{2}}{4} -\frac{x}{2} -1$",
				'approx_solution' : -0.605551
				},
		'B' : {
				'method' : lambda x:  np.cos(x),
				'initial guess' : 2,
				'latex' : r"$cos(x)$",
				'approx_solution' : 0.739085
				},

		'C' : {
				'method' : lambda x: (x/3 +x)/2,
				'initial guess' : 2,
				'latex' : r"$\frac{1}{2}\left(\frac{3}{x} + x \right)$",
				'approx_solution' : 0
				},

		'D' : {
				'method' : lambda x: np.cosh(x)/x -np.arctan(x),
				'initial guess' : 2,
				'latex' : r"$\frac{cosh(x)}{x} -arctan(x)$",
				'approx_solution' : 0.881709
				}
	}

	# plot fixed point iterations
	order = ['A', 'B', 'C', 'D']
	domain = (-2,2)
	num_fixed_point_iterations = 10

	# For all problems use an absolute solution tolerance of 10−5, 
	# a relative solution tolerance of 10−10, 
	# allow a maximum of 100 iterations.
	maxit = 100
	atol = 10**(-5)
	rtol = 10**(-10)

	analysis_plot(fixed_point_functions, order, domain, atol, rtol, maxit)

	# interp_representation_error
	interp_representation_error = lambda Gfun, x_final: abs(Gfun(x_final) - x_final)
	
	# run trials
	for name in order:

		Gfun = fixed_point_functions[name]['method']
		x0 = fixed_point_functions[name]['initial guess']

		if SHOW_OUTPUT:
			logger.info("\n\n\tRunning fixed point iterations WITHOUT aitken acceleration on problem %s.\n\n\tFunction %s.\n" % (name,name))

		x_final = fixed_point(Gfun, x0, maxit, rtol, atol, output=SHOW_OUTPUT)

		if SHOW_OUTPUT:
			# display interp_representation_error
			error = interp_representation_error(Gfun, x_final)
			logger.info("\n\n[Method 2]: \nFinal fixed point target error WITHOUT using acceleration acceleration.\n \
				\n\tFunction name: %s, \n\n\t\t |Gfun(x) - x| = %s.\n" % (name, error))

		if SHOW_OUTPUT:
			logger.info("\n\n\tRunning fixed point iterations WITH aitken acceleration on problem %s.\n\n\tFunction %s.\n" % (name,name))

		x_final = aitken_fp(Gfun, x0, maxit, rtol, atol, output=SHOW_OUTPUT)

		if SHOW_OUTPUT:
			# display interp_representation_error
			error = interp_representation_error(Gfun, x_final)
			logger.info("\n\n[Method 2]: \nFinal fixed point target error using acceleration acceleration.\n \
				\n\tFunction name: %s, \n\n\t\t |Gfun(x) - x| = %s.\n" % (name, error))

	# Analysis
	#logger.info("\n\n.\n \
	#\n\tFunction name: %s, \n\n\t\t |Gfun(x) - x| = %s.\n")


