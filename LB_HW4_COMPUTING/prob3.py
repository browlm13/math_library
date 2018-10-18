#!/usr/bin/env python3

"""
File name: prob3.py

	Python Version: 3.6

		Script to compare hermite polynomial interpolation 
		and cubic spline interpolationof uniformly-spaced nodes.

LJ Brown
Fall 2018
"""

__filename__ = "prob3.py"
__version__ = 2.0
__author__ = "L.J. Brown"

#
# imports
#

# internal library
import logging

#
# my library
#

# implimented using class specs
from cubic_spline import *

# implimented using scipy conventions
from LagrangeInterpolatingPolynomial import *
from HermiteInterpolatingPolynomial import *
from CubicSpline import *

# external library
import numpy as np
import matplotlib.pyplot as plt

# initilize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

#
# testing
#

#
#   Define test problems, functions and interpolation methods
#


test_function_1_params = {

	'function_names' : { 
							'python_string' : "arctan(2x^2)", 
							'latex_string' : '$f(x) = tan^{-1}(2x^{2})$'
					  }, 

	'f'  : lambda t: np.arctan(2*(t**2)),                       # arctan(2*(x^2))
	'df' : lambda t: (4*t)/( 4*(t**4) + 1),                     # 4x/( 4(x^4) +1 )
	'ddf': lambda t: ( 4*( -12*(t**4) +1 ) )/( 4*(t**4) +1 )    # [4(-12(x^4) +1)]/( 4(x^4) +1 )
}

test_functions = [test_function_1_params]

# general parameters
L = 3.0
domain = (-L, L)
polynomial_knots = [ 41, 21, 11, 5]
evaluation_xs = np.linspace(-L, L, num=400)

interpolation_methods = {
	
	'lagrange' : {
					'InterpolantClass' : LagrangeInterpolatingPolynomial,
					'params' : lambda node_xs, f : dict(t=node_xs, y=f(node_xs)),
					'color' : 'green',
					'marker' : '1',
					'line' : 'gs'
				},

	'hermite' : {
					'InterpolantClass' : HermiteInterpolatingPolynomial,
					'params' : lambda node_xs, f : dict(t=node_xs, y=f(node_xs), dy=f(node_xs, deriv_order=1)),
					'color' : 'blue',
					'marker' : '2',
					'line' : 'b^'
				},

	'cubic_spline' : {
					'InterpolantClass' : CubicSpline,
					'params' : lambda node_xs, f : dict( t=node_xs, y=f(node_xs), start_bc=(1, f.df(node_xs[0])), end_bc=(2, f.ddf(node_xs[-1]))),
					'color' : 'red',
					'marker' : '3',
					'line' : 'r--'
				}
}

#
#  helper methods
#

class Function:

	def __init__(self, function_params):

		# Test function python string title and latex title
		# for logging and graphs clarity

		self.string_title = function_params['function_names']['python_string']
		self.latex_title = function_params['function_names']['latex_string']

		# Test function, its first, and its second derivatives
		self.f   = function_params['f']
		self.df  = function_params['df']
		self.ddf = function_params['ddf']

	def __call__(self, xs, deriv_order=0):

		if deriv_order == 0:
			results = [ self.f(x) for x in xs ] 

		if deriv_order == 1:
			results = [ self.df(x) for x in xs ] 

		if deriv_order == 2:
			results = [ self.ddf(x) for x in xs ] 

		if len(results) == 1:
			return results[0]
		else:
			return np.array(results).flatten()

		# raise error
		assert 1 == 0

	# try not to use
	def get_derivatives(self):
		return self.f, self.df, self.ddf

def uniform_nodes_for_interpolant(domain, degree, InterpolantClass):

	num_nodes = InterpolantClass.num_nodes_given_degree(degree)
	uniform_nodes = np.linspace( domain[0], domain[1], num=num_nodes)

	return uniform_nodes

if __name__ == "__main__":

	for i, function_params in enumerate(test_functions):

		f = Function(function_params)
		#f, df, ddf = func.get_derivatives()

		fig, axarr = plt.subplots(len(polynomial_knots),2)

		for i, num_knots in enumerate(reversed(polynomial_knots)):

			plot1_title = '$Knots_{' + str(num_knots) + '}$'
			axarr[i,0].set_title(plot1_title, fontsize='small')

			plot2_title = '$Error_{' + str(num_knots) + '}$'
			axarr[i, 1].set_title(plot2_title, fontsize='small')

			# plot true function line
			fplot = axarr[i,0].plot(evaluation_xs, f(evaluation_xs), linewidth= 3, alpha=0.75, color='k', label=f.latex_title)


			for method_name in interpolation_methods.keys():
				
				InterpolantClass = interpolation_methods[method_name]['InterpolantClass']


				# get number of nodes needed for degree generically
				node_xs = uniform_nodes_for_interpolant(domain, num_knots, InterpolantClass)

				# build generic interpolant
				class_params = interpolation_methods[method_name]['params'](node_xs, f) # get class parameters dictonary
				interpolant = InterpolantClass(**class_params)


				#
				# plotting
				#

				# get class color
				color = interpolation_methods[method_name]['color']

				#
				# Direct overlay
				#

				# plot line graph
				iplot = axarr[i,0].plot(evaluation_xs, interpolant(evaluation_xs),  color=color, linestyle='--', label=interpolant.latex_title)
				
				#
				# Error Plots
				#

				iplot_error = axarr[i,1].plot(evaluation_xs, abs(f(evaluation_xs)-interpolant(evaluation_xs)), color=color, linestyle='-', label=interpolant.latex_error_title)

		#
		# log analysis
		#

		analysis_string = """

			The degrees of the hermite and lagrange 
		interpolating polynomials increase with the 
		number of knots to interpolate. The runge 
		phenominon begins to take effect on high d-
		egree polynomials especially when the knots 
		to interpolate are evenly spaced. This can 
		be observed as unexpected oscillations in 
		the interpolation at its tail ends. These 
		oscillations increase in amplitude with 
		the degree of the interpolating polynomia-
		ls.
		\n
			We see the runge phenominon begining to 
		effect the error of the lagrange interpola-
		-ting polynomial with as few as 11 knots. 
		The tail ends of the lagrange interpolant go 
		wild when we reach 21 knots. At 21 knots the 
		lagrange interpolant error is so high it ob-
		-scures the error of the other two interpol-
		-ants. 
		\n
			The hermite interpolating polynomials 
		degree also increases with the number of kn-
		-ots to interpolate. In this set up the her-
		-mite polynomial is of the same degree as t-
		-he lagrange polynomial in each graph. The 
		advantage we see with the hermite interpola-
		-ting polynomial comes from matching the fi-
		-rst derivative at each knot it interpolates. 
		This additional interpolation condition is 
		able to reduce the negitive effects of the 
		runge phenominon, but only up until a point. 
		At 41 knots the hermite interpolant also suf-
		-fers massivley at its tail ends.
		\n
			The cubic spline on the other hand is 
		always composed of 3rd degree polynomials 
		regarless of the number of knots it is 
		interpolating. The splines interpolation 
		is able to greatly reduce oscillation be-
		-tween data points. It ensures continuity 
		of both the first and second derivatives 
		at each knot and so it ensures smoothness 
		at each point where the polynomials join.
		As the number of knots increase, the spa-
		-cing between the knots decreases, and 
		the cubic spline becomes increasingly ac-
		-curate. When interpolating only a few 
		knots seperated by relativly large dista-
		-nces the cubic spline interpolation is 
		at its weakest. This is the oposite of 
		the other two interpolation methods.

		"""

		logger.info("\n\n" + analysis_string)

		#leg = axarr[0,0].legend(loc='lower center', shadow=True, bbox_to_anchor=(0, 0), fancybox=True)
		leg = axarr[1,0].legend(bbox_to_anchor=(1.1, 0, 0, 0), loc=3,
       		ncol=1,  shadow=True, fontsize='small', fancybox=True) #mode="expand", borderaxespad=0,  shadow=True,  fancybox=True)

		leg = axarr[2,1].legend(bbox_to_anchor=(-0.85, 0, 0, 0), loc=3,
       		ncol=1,  shadow=True, fontsize='small', fancybox=True) #mode="expand", borderaxespad=0,  shadow=True,  fancybox=True)

		plt.subplots_adjust(top=0.95, bottom=0.10, left=0.10, right=0.95, hspace=0.65,
                    wspace=0.95)

		# display plots
		plt.show()

