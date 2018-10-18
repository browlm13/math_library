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
from lagrange import *
from cubic_spline import *
from hermite import *

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
#   Define test problems and functions
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

def uniform_nodes_for_interpolant(domain, degree, InterpolantClass):

	num_nodes = InterpolantClass.num_nodes_given_degree(degree)
	uniform_nodes = np.linspace( domain[0], domain[1], num=num_nodes)

	return uniform_nodes

if __name__ == "__main__":

	for i, function_params in enumerate(test_functions):

		f = Function(function_params)
		#f, df, ddf = func.get_derivatives()

		fig, axarr = plt.subplots(len(polynomial_knots),2)

		for i, num_knots in enumerate(polynomial_knots):

			plot1_title = '$Knots_{' + str(num_knots) + '}$'
			axarr[i,0].set_title(plot1_title, fontsize='small')

			plot2_title = '$Error_{' + str(num_knots) + '}$'
			axarr[i, 1].set_title(plot2_title, fontsize='small')

			#axarr[i,0].set_xlabel('x')
			#axarr[i,0].set_ylabel('y')
			#axarr[i,1].set_xlabel('x')
			#axarr[i,1].set_ylabel('y')


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


		#leg = axarr[0,0].legend(loc='lower center', shadow=True, bbox_to_anchor=(0, 0), fancybox=True)
		leg = axarr[1,0].legend(bbox_to_anchor=(1.1, 0, 0, 0), loc=3,
       		ncol=1,  shadow=True, fontsize='small', fancybox=True) #mode="expand", borderaxespad=0,  shadow=True,  fancybox=True)

		leg = axarr[2,1].legend(bbox_to_anchor=(-0.85, 0, 0, 0), loc=3,
       		ncol=1,  shadow=True, fontsize='small', fancybox=True) #mode="expand", borderaxespad=0,  shadow=True,  fancybox=True)

		plt.subplots_adjust(top=0.95, bottom=0.10, left=0.10, right=0.95, hspace=0.65,
                    wspace=0.95)

		# display plots
		plt.show()

""""
problems:

Runge phenominon - especially for lagrange

Flat regions for lagrange and hermite, but not spline.
I think it becomes difficult when the degree of the polynomial is high, this is not true for cubic spline interpolation

   
If there are a large number of data points use piecewise polynomial interpolation
Use piecewise polynomials instead of using a single polynomial (of high degree) to interpolate these points

A spline of degree m is a piecewise polynomial (of degree m) with the maximum possible 
smoothness at each of the points where the polynomials join

optimality theorem?

Polynomials are not well suited to interpolation problems because of Runge's phenomenon -- 
the short version is that if you try to interpolate values using a polynomial, 
you get unexpected oscillations in the interpolation, and using higher degree polynomials can make things worse: 

Say you have the red function above, and you try to interpolate using a fifth order polynomial 
(the blue line) -- that looks pretty rough. So you say, hmmm, let me try a 9th order polynomial. 
And what you find is that it actually makes things worse -- there are weird oscillations between the 
interpolation points. 

A second reason when used in modelling is that cubic splines have local control. 
When you modify a control point for a spline, it only affects the area around the control points 
(for a cubic spline, it is limited by the convex hull of the two prior and two subsequent control points). 
Polynomials: you change one interpolation point and you can change the curve anywhere along the curve 
(e.g. in the above picture, you take the minimum at x=-0.95 and change it, the local minimum at x=0.95 
	will also shift.

for a lot of data sets it can be beneficial. The problem with having lots of data, 
especially if it’s roughly equally spaced apart, is that polynomial interpolation suffers from Runge’s Phenomena. 
This means that as you add more data, the derivatives at each of the data points tend to grow. 
This results in large oscillations between data points that typically don’t tend to be “right”.

Spline interpolation tends to greatly reduce oscillation between data points.
 Part of this is in the derivation of the splines. For typical 3rd order splines, 
 you can actually derive the spline equations using the 1D Beam Equations. Given that, 
 you can use some physical intuition to understand why the splines may have less error between data. 
 You can view the spline result as trying to bend a beam to touch each of the data points. 
 Since it’s a beam, it will deflect a minimal amount between the data points, doing just enough to 
 get to the next data point.

 Thus, spline interpolation often does not have as much error as polynomial interpolation. 
 However, it is possible to space data out in such a way that the polynomial fit will actually perform well.
  But splines do a pretty good job for a generic data set (assuming no outliers), so it’s a bit more reliable.

  One can use, Lagrange's Interpolation formula when there is an uneven spread in the domain set.
Example: When sets such as {0,1,4,8,11,12,17} or {1,2,4,6,9,10} are the domains of functions.
"""
