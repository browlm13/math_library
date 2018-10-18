#!/usr/bin/env python3

"""
File name: CubicSpline.py

	Python Version: 3.6

		==========================
		Cubic Spline Interpolation.
		==========================

			----------------------------------------

				'n+1' data points to interpolate:

			----------------------------------------

				Given numpy arrays t and y:

					t = [t_0, ..., t_n] ~ 'knots', n+1 data points (paramater 1)
					y = [y_0, ..., y_n] ~ 'f(t)',  n+1 data points (parameter 2)

						* where y_i = f(t_i)

									Table:
									------

					t ||  t_0   |  t_1   |  ...  |  t_n   |  
					---------------------------------------
					y || f(t_0) | f(t_1) |  ...  | f(t_n) | 

			------------------------------------------------------------------

				Defining the interpolating spline function S of degree k,
				using 'n' polynomials (S_0, ..., S_n-1) of at most degree k

						* where k = 3 for the cubic spline S

			------------------------------------------------------------------

					The 'n+1' data points are interpolated using 'n' cubic 
				splines (S_0, ..., S_n-1). Let S be the cubic spline function 
				constructed to interpolate the table data above.

					On each of the n subintervals: [t_0, t_1], [t_1, t_2], ..., [t_n-1, t_n], 
				let S_i be the cubic polynomial that represents S on [t_i, t_i+1].

					* S ~  
							Spline function of degree '3'.
					   S is a peicewise polynomial that interpolates 
					   all 'n+1' data points. 

					* S_0, ..., S_n-1 ~ 

							'n' polynomials of at most degree '3' 
						that together form the interpolting cubic spline S. 
						Each of the n cubic polynomials has k+1 = 4 coefficients,
						a_i0, ..., a_i3,

						S_i =  (a_i0)  +(a_i1)*(x)  +(a_i2)*(x^2)  +(a_i3)*(x^3)


								{  
								   S_0( x )   	for x in [t_0, t_1),
								   S_1( x )   	for x in [t_1, t_2),

					S( x ) =	   .
								   .
								   .

								   S_n-1( x ) 	for x in [t_n-1, t_n]  
								}

						* S has a total of '4n' coefficients, 4 from each of 
						the n polynomials that form S.

						* (k+1)n coefficients in general

				------------------------------------------------------------------------

					Our goals:

					1) Interpolation Condition:
					   Construct S to interpolate all n+1 datapoints

								S(t_i) = y_i, for all t_i

					2) Continuity Conditions:
					   Construct S to have continous derivatives of all orders up k-1

							denoted,
									 " S in C^k-1 "

				-------------------------------------------------------------------------

					---------------------------
					1) Interpolation Condition:
					---------------------------
					- - - - - - - - - - - - - - - - - - - - - 
					Condition 0 : S(t_i) = y_i for all t_i :
					- - - - - - - - - - - - - - - - - - - - - 

				   To construct S so that it interpolate all n+1 data points we 
				   assign the constraint,

						S(t_i) = y_i for all t_i

					There are 2 interpolation conditions on 
					each of the 'n' subinterval,
							
							S(t_i) = y_i,
							S(t_i+1) = y_i+1

					* This gives rise to '2n' constraints.

					--------------------------------------------
					2) Continuity Conditions (Conditions 1,2,3):
					--------------------------------------------

						To construct a cubic spline S so that,

										S in C^2 

								( k-1 = 2 so S must have continous 
									derivatives of all orders up to 2 )

						we must satisfy 3 conditions:

						(Condition 1) : 

								S in C^0, 
								" S must be continous "

						(Condition 2) : 

								S in C^1, 
								" S's first derivative, S', must be continous "

						(Condition 3) : 

								S in C^2, 
								" S's second derivative, S'', must be continous "


					- - - - - - - - - - - - - 
					 Condition 1 - S in C^0 :
					- - - - - - - - - - - - - 

						The polynomials making up S are all continous on the 
					subintervals they are defined, 

								S_i(x) continous for all x in [t_i-1, t_i)

						So we will focus on the interior knots, t_1, ..., t_n-1. If S is continous here then 
					Condition 1 will be met.

						The polynomials S_i-1 and S_i interpolate the 
					same value y_i at the point t_i making S continous at all n-1 interior knots,
				
									S_i-1( t_i ) = y_i = S_i( t_i )   eq (2)

					Then S is continous by eq (2) above. 

									S in C^0 ~ (Condition 1 met) 

					* This continuity does not give rise to any 
					additional constraints because it has already 
					been counted in Conition 0 (the interpolation condition)

					- - - - - - - - - - - - - 
					 Condition 2 - S in C^1 :
					- - - - - - - - - - - - - 

						The first derivatives of the polynomials making up S are all continous on the open
					subintervals they are defined because they are degree 3 polynomials,

									S_i'(x) continous for all x in (t_i-1, t_i)

						To ensure that S' is continous we just need to focus on the 
					'n-1' interior knots, t_1, ..., t_n-1, and add the constraints that,
									
									S_i-1'(t_i) = S_i'(t_i)  eq (3)

					Then S' is continous by eq (3) above. 

					* This gives an additional 'n-1' constraints

					- - - - - - - - - - - - - 
					 Condition 3 - S in C^2 :
					- - - - - - - - - - - - - 

						The second derivatives of the polynomials making up S are also all continous 
					on the open subintervals they are defined because they are degree 3 polynomials,

									S_i''(x) contious for all x in (t_i-1, t_i)

						Again in order to ensure that S'' is continous we just focus on the 'n-1' 
					interior knots, t_1, ..., t_n-1, and add the constraints that,
									
									S_i-1''(t_i) = S_i''(t_i)  eq (4)

					Then S'' is continous by eq (4) above. 

					* This again gives an additional 'n-1' constraints

				------------------------------------------

					Underdetermined System:

						As it stands we have, 

							4n coefficients with 4n -2 constraints

				-------------------------------------------

					By constructing S to satisfy conditions 0, 1, 2 and 3 
				we arrive at a total of '4n -2' constraints,
						
							2n + n-1 + n-1 = 4n -2

						--------------------------
						Conditions | # Constraints
						--------------------------
							0 	   -> 	2n  (interpolation condtion)
						(X) 1 	   -> 	0   (n-1, already counted)
							2 	   -> 	n-1
						+	3 	   -> 	n-1
						-------------------------
					total constraints:  4n -2

					Altogether the conditions above lead to 4n -2 constraints 
				for determining 4n coeffiencients, leaving two degrees of freedom 
				that can be chosen inteligently to assist with the particular 
				interpolation problem at hand.
				
					[Side Note]

							In general for a spline of degree k 
						and n+1 data points to interpolate:

							* 2n + (k-1)*(n-1) total constraints

									n-1 additional conditions are added 
									each time k increments, 
									(1 for each interior knot) 
									k-1 times so S in C^k-1

							* (k+1)*n total coefficients 

									(k+1)*n coefficients need to be 
									determined for a spline of degree k

				------------------------------------------------------------

					Boundry Conditions:

					chosing additional constraints for the underdetermined
					system of polynomial equations.

						4n (-2) constraints vs. 4n coeffiencients

					* 2 degrees of freedom for the cubic spline (k=3)

				------------------------------------------------------------


				--------------------------------------------------------------

				Deriving the equation for S_i(x) on the interval [t_i, t_i+1].

				--------------------------------------------------------------

					S'' is continous at each interior knot, 
				and S_i is a cubic polynomial on [t_i, t_i+1].

					let S_i''(t_i) = z_i,
					let h_i = t_i+1 - t_i

				Then S_i''(x) is a linear function satisfying,

					S_i''(t_i)   = z_i,
					S_i''(t_i+1) = z_i+1

				Then S_i''(x) is given by the straight line between z_i and z_i+1,

								z_i 			   z_i+1
					S_i''(x) =  ---( t_i+1 -x )  + -----( x - t_i ) 
								h_i 				h_i

				After integrating twice and imposing the interpolation conditions 
				S_i(t_i) = y_i and S_i(t_i+1) = y_i+1 you get an equation for S_i(x),

						------
						eq (*):
						------

									 z_i 			        z_i+1
						S_i(x) =  + ------( t_i+1 -x )^3  + ------( x - t_i )^3  +
									 6*h_i 				     6*h_i

									 [ 	y_i+1	 (z_i+1 * h_i) ]
								   + | 	-----  -  -----------  |( x - t_i ) +
									 [   h_i  		   6       ]

									 [ 	 y_i	 (z_i * h_i) ]
								   + | 	-----  - ----------- |( t_i+1 - x )
									 [   h_i  		  6      ]


				Once z_0, ..., z_n are found, eq (*) can be used in 
				conjunction with the definition of S(x) to evaluate 
				the cubic spline at any point.


		Determine z_1, ..., z_n-1:

			The continuity conditions for S' at the interior knots imply,

				S_i-1'(t_i) = S_i'(t_i)

			S_i'(x) can be obtained by differentiating eq (*), 
			after substituting x=t_i the simplified equation for S_i'(t_i) becomes,

									 h_i			      h_i        y_i       y_i+1
				S_i'(t_i) =  - (z_i) -----    - (z_i+1)  -----    - -----    + -----
									  3  				   6 		 h_i 		h_i

			Similarily for S_i-1'(t_i),


										h_i-1			 h_i-1       y_i-1         y_i
				S_i-1'(t_i) =   (z_i-1) -----    + (z_i) -----    - -------    + ------
										  6  			   3 	     h_i-1 		  h_i-1

			Setting these two equations equal to one another and the result can be written as,


																	6					 6
			(h_i-1)(z_i-1) + 2(h_i + h_i-1)(z_i) + (h_i)(z_i+1) = ----(y_i+1 - y_i)  - -----(y_i - y_i-1)   eq($)
																   h_i 				   h_i-1

			let, 
				u_i = 2(h_i + h_i-1),
				b_i = [ 6*(y_i+1 - y_i) ]/(h_i),
				v_i = b_i - b_i-1

			then eq ($) becomes,

				(h_i-1)(z_i-1) + (u_i)(z_i) + (h_i)(z_i+1) = v_i,     for i = 1, ..., n-1

				or,

					< h_i-1,  u_i,  h_i >  dot  < z_i-1,  z_i,  z_i+1 >  =  v_i

			written as a system of equations,


				C is an (n-1) x (n-1) tridagonal matrix,

				C z_c = v_c


				|u1 h1 						   | |  z1   |     |  v1   |
				|h1 u2 h2 					   | |  z2   |     |  v2   |
				|   h2 u3 h3  				   | |  z3   |     |  v3   |
				|       ... 				   | |   .   |     |  .    |
				|         ...  				   | |   .   |  =  |  .    |
				|           ...  			   | |   .   |     |  .    |
				|            h_n-3 u_n-2 h_n-2 | | z_n-2 |     | v_n-2 |
				|                  h_n-2 u_n-1 | | z_n-1 |     | v_n-1 |




			This equation is used only for i = 1, ..., n-1, giving a system of n-1 
			linear equations for the n+1 unkown z_0, ..., z_n. 

			z_0 and z_n can be selected arbitrarily because of the two 
			degrees of freedom at your disposal.


		   Boundary conditions (2.0):

			start_bc = (deriv_order, deriv_value)
			end_bc = (deriv_order, deriv_value)

			deriv_order ~ derivative order, 1 or 2
			deriv_value ~ value of derivative

			#
			# start_bc
			#

			if deriv_order == 1:
				S'(t_0) = deriv_value

				2*h_0 * (z_0) + h_0 * (z_1) = b_0 -6 * deriv_value

				or first row of matrix A becomes:

					A[0,:] = [2*h_0, h_0, 0, ..., 0] 

				and first element of vector r becomes:

					r[0] = b_0 -6 * deriv_value

			if deriv_order == 2:
				S''(t_0) = deriv_value

				(z_0) = deriv_value

				or first row of matrix A becomes:

					A[0,:] = [1, 0, ..., 0] 

				and first element of vector r becomes:

					r[0] = deriv_value

		   #
		   # end_bc
		   #

			if deriv_order == 1:
				S'(t_n) = deriv_value

				2*h_n-1 * (z_n-1) + h_n-1 * (z_n) = b_n-1 -6 * deriv_value

				or last row of matrix A becomes:

					A[n,:] = [0, ..., 0, 2*h_n-1, h_n-1] 

				and first element of vector r becomes:

					r[n] = b_n-1 -6 * deriv_value

			if deriv_order == 2:
				S''(t_n) = deriv_value

				(z_n) = deriv_value

				or last row of matrix A becomes:

					A[n,:] = [0, ..., 0, 1] 

				and last element of vector r becomes:

					r[n] = deriv_value


			----------
			HW Problem:
			-----------

			The given boundry conditions for problem 3 become:

			S'(t_0) = alpha
			S''(t_n) = beta

			start_bc = (1, alpha)
			end_bc = (2, beta)


			1.) 2*h_0 * (z_0) + h_0 * (z_1) = b_0 -6 * alpha

				first row:

					[2*h_0,  h_0, 0, ..., 0 ] * [z_0, z_1, ..., z_n].T = [ (b_0 -6*alpha), v_1, ..., v_n-1, beta]

			2.) (z_n) = beta

				last row:

					[0, ..., 1] * [z_0, ..., z_n].T  = [ (b_0 -6*alpha), v_1, ..., v_n-1, beta]
			

LJ Brown
Fall 2018
"""

__filename__ = "CubicSpline.py"
__author__ = "L.J. Brown"

# internal libraries
import logging

# external libraries
import numpy as np

# initilize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class CubicSpline:

	"""
		use example:

			# test function and its derivatives:
			f = lambda t: np.arctan(2*(t**2)) 						# arctan(2*(x^2))
			df = lambda t: (4*t)/( 4*(t**4) + 1) 					# 4x/( 4(x^4) +1 )
			ddf = lambda t: ( 4*( -12*(t**4) +1 ) )/( 4*(t**4) +1 ) # [4(-12(x^4) +1)]/( 4(x^4) +1 )

			# general parameters
			L = 3.0 # range (-L, L) for nodes
			n = 41  # number of nodes

			# generate nodes
			t = np.linspace(-L, L, num=n)
			y = f(t)

			# boundry conditions from HW problem 3
			alpha = df(t[0])
			beta = ddf(t[-1])

			# construct cubic spline class passing nodes and boundry conditions
			cs = CubicSpline(t, y, start_bc=(1,alpha), end_bc=(2,beta))

			# evaluation points to test cubic spline after construction
			xs = np.linspace(-L, L, 800)
		
			# evaluate
			yhats = cs(xs) 

	"""

	def __init__(self, t, y, start_bc=(2, 0.0), end_bc=(2, 0.0)):
		self.t = t
		self.y = y 
		self.start_bc = start_bc
		self.end_bc = end_bc

		# calculate coefficients
		self.z = cubic_spline_coefficients_bc(self.t, self.y, start_bc=self.start_bc, end_bc=self.end_bc)

		n = len(self.t)
		self.degree = n

		self.latex_title = self.get_latex_title()
		self.latex_error_title = self.get_latex_error_title()

	def __call__(self, xs):
		yhats = cubic_spline_evaluate(self.t, self.y, self.z, xs)
		return yhats

	@staticmethod
	def num_nodes_given_degree(degree):
		return degree

	@staticmethod
	def get_latex_title():
		title = '$Spline_{3}(x)$'
		return title
		
	@staticmethod
	def get_latex_error_title():
		error_title = '$ | f(x) - Spline_{3}(x) | $'
		return error_title


def cubic_spline_coefficients_bc(t, y, start_bc=(2, 0.0), end_bc=(2, 0.0)):
	"""
		Determine coeffiencients, z, of cubic spline. Takes nodes (t, f(t)) ~ where f(t) = y.

		Solve system:

				Az = r, 

				where, 
					A is an (n+1)x(n+1) matrix,
					z is an unkown vector of length n+1,
					r is a vector of length n+1

			Boundry conditions,

			start_bc = (deriv_order, deriv_value)
			end_bc = (deriv_order, deriv_value)

			deriv_order ~ derivative order, 1 or 2
			deriv_value ~ value of derivative

			#
			# first row determined by start_bc (starting boundry conditions)
			#

			# deriv_order == 1 then (S'(t_0) = deriv_value)

			|2h0 h0     						   | |  z0   |     | b0 -6*deriv_value |
									  or
			# deriv_order == 2 then (S''(t_0) = deriv_value)
			| 1          						   | | z_n   |     | deriv_value |

			#
			#	tridiagonal system for interior nodes
			#

			| h0 u1 h1 						   	   | |  z1   |     |  v1   		 |
			|    h1 u2 h2 					   	   | |  z2   |     |  v2   		 |
			|       h2 u3 h3  				       | |  z3   |     |  v3   		 |
			|          ... 				           | |   .   |     |  .    		 |
			|             ...  				       | |   .   |  =  |  .    		 |
			|              ...  			       | |   .   |     |  .    		 |
			|               h_n-3 u_n-2 h_n-2      | | z_n-2 |     | v_n-2 		 |
			|                     h_n-2 u_n-1 h_n-1| | z_n-1 |     | v_n-1 		 |
			
			#
			# last row determined by end_bc (ending boundry conditions)
			#

			# deriv_order == 1 then (S'(t_n) = deriv_value)

			|    					2*h_n-1, h_n-1 | |  zn   |     | bn-1 -6*deriv_value |

									  or
			# deriv_order == 2 then (S''(t_n) = deriv_value)
			|           						 1 | | z_n   |     | deriv_value |

			where,

				ui = 2(hi + hi-1),
				bi = [ 6*(yi+1 - yi) ]/(hi),
				vi = bi - bi-1


		:param t: numpy array of nodes first element.
		:param y: numpy array of nodes second emelemnt ( true function evalutated at t, f(t))
		:param start_bc: tuple, (deriv_order, deriv_value) specify boundry conditions for the start node of the cubic spline, where deriv_order is 1 or 2 for cubic spline and deriv_value is its chosen value (default natural second derivative equal to 0)
		:param end_bc: tuple, (deriv_order, deriv_value) specify boundry conditions for the end node of the cubic spline, where deriv_order is 1 or 2 for cubic spline and deriv_value is its chosen value (default natural second derivative equal to 0)
		:returns: numpy array of coefficients for each cubic polynomial (4)*** for each interval, order corresponds to node ordering.
	"""

	# len(t) = n+1
	n = len(t) -1

	# build tridiagonal matrix initially filled with zeros
	A = np.zeros(shape=(n+1, n+1))

	# build vector r from equation Az = r
	r = np.zeros(shape=(n+1,))

	# h_i = t_i+1 - t_i
	h = lambda i: t[i+1] - t[i] 

	# u_i = 2(h_i + h_i-1)
	u = lambda i: 2*(h(i) + h(i-1))

	# b_i = (6/h_i) * (y_i+1 - y_i)
	b = lambda i: ((6)*(y[i+1] - y[i]))/(h(i))

	# v_i = b_i - b_i-1
	v = lambda i: b(i) - b(i-1)

	#
	# starting soundry conditions
	#

	# if condition for first node is S'(t_0) = deriv_value
	if start_bc[0] == 1:

		# A row 0:   [2h0, h0, 0, ... , 0]
		A[0,0], A[0,1] = 2*h(0), h(0)

		# first element in r = b(0) -6 * deriv_value
		r[0] = b(0) -6 * start_bc[1]

	# if condition for first node is S''(t_0) = deriv_value
	if start_bc[0] == 2:
		# A row 0:   [1, 0, ... , 0]
		A[0,0] = 1

		# first element in r = deriv_value
		r[0] = start_bc[1]

	#
	# interior nodes system
	#

	# fill in equations for intererior nodes rows 1->n-1
	for i in range(1,n):

		# fill tridagonal system
		A[i, i-1] = h(i-1)
		A[i, i] = u(i)
		A[i, i+1] = h(i)

		# fill result vector neglecting first 2 and last 2 elements
		r[i] = v(i)
	
	#
	# ending boundry conditions
	#

	# if condition for first node is S'(t_n) = deriv_value
	if end_bc[0] == 1:

		# A row n:   [0, ..., 0, 2hn-1, hn]
		A[n,n-1], A[n,n] = 2*h(n-1), h(n)

		# last element in r = b(n-1) -6 * deriv_value
		r[n] = b(n-1) -6 * end_bc[1]

	# if condition for first node is S''(t_n) = deriv_value
	if end_bc[0] == 2:
		
		# A row n:   [0, ..., 1]
		A[n,n] = 1

		# last element in r:
		r[n] = end_bc[1]

	# solve system for coefficients, z
	z = np.linalg.solve(A, r)

	return z


def cubic_spline_evaluate(t, y, z, xs):
	"""

				{  
				   S_0( x )   	for x in [t_0, t_1),
				   S_1( x )   	for x in [t_1, t_2),

	S( x ) =	   .
				   .
				   .

				   S_n-1( x ) 	for x in [t_n-1, t_n]  
				}

	where S_i is defined,

					 z_i 			        z_i+1
		S_i(x) =  + ------( t_i+1 -x )^3  + ------( x - t_i )^3  +
					 6*h_i 				     6*h_i

					 [ 	y_i+1	 (z_i+1 * h_i) ]
				   + | 	-----  -  -----------  |( x - t_i ) +
					 [   h_i  		   6       ]

					 [ 	 y_i	 (z_i * h_i) ]
				   + | 	-----  - ----------- |( t_i+1 - x )
					 [   h_i  		  6      ]

	:param t: numpy array of n+1 values
	:param y: numpy array of n+1 values of the function to interpolate evaluated at the points t.
	:param z: numpy array of n+1 coefficients for evaluating cubic spline, or second derivatives at the knots.
	:param xs: numpy array of floats containing points to evaluate cubic spline at.
	:returns: numpy array of floats of cubic spline evaluated at all points in xs.

	"""

	# len(t) = n+1
	n = len(t)-1

	# useful interval width function takes interval indexs as parameter
	h = lambda i: t[i+1] - t[i]

	# evalutation method given interval and x value
	def S(i,x):

		term1 = (z[i] * (t[i+1] - x)**3)/(6*h(i)) 
		term2 = (z[i+1] * (x - t[i])**3)/(6*h(i)) 
		term3 = ((y[i+1]/h(i)) - ((z[i+1]*h(i))/6)) * (x - t[i])
		term4 = ((y[i]/h(i)) - ((z[i]*h(i))/6)) * (t[i+1] -x)

		return term1 + term2 + term3 + term4

	# evaluate using polynomial in correct interval
	def evaluate(x):
		for i in range(1, n):
			if x < t[i]: 
				return S(i-1,x) 
		return S(n-1, x)

	# evaluate all points in xs and add to results vector
	results = [evaluate(x) for x in xs]

	return np.array(results)


def cubic_spline_coefficients(t, y, alpha, beta):
	"""
		Determine coeffiencients, z, of cubic spline. Takes nodes (t, f(t)) ~ where f(t) = y.

		Solve system:

				Az = r, 

				where, 
					A is an (n+1)x(n+1) matrix,
					z is an unkown vector of length n+1,
					r is a vector of length n+1

			|2h0 h0     						   | |  z0   |     | b0 -6*alpha |
			| h0 u1 h1 						   	   | |  z1   |     |  v1   		 |
			|    h1 u2 h2 					   	   | |  z2   |     |  v2   		 |
			|       h2 u3 h3  				       | |  z3   |     |  v3   		 |
			|          ... 				           | |   .   |     |  .    		 |
			|             ...  				       | |   .   |  =  |  .    		 |
			|              ...  			       | |   .   |     |  .    		 |
			|               h_n-3 u_n-2 h_n-2      | | z_n-2 |     | v_n-2 		 |
			|                     h_n-2 u_n-1 h_n-1| | z_n-1 |     | v_n-1 		 |
			|           						1  | | z_n   |     | beta  		 |

			where,

				ui = 2(hi + hi-1),
				bi = [ 6*(yi+1 - yi) ]/(hi),
				vi = bi - bi-1

			Boundry equations,

				S'(t_0) = alpha
				2*h_0 * (z_0) + h_0 * (z_1) = b_0 -6 * alpha

				row 0: |2h0 h0     						   | |  z0   |   =  | b0 -6*alpha |

				S''(t_n) = beta
				(z_n) = beta

				row n: |           						1  | | z_n   |   =  | beta |

		:param t: numpy array of nodes first element.
		:param y: numpy array of nodes second emelemnt ( true function evalutated at t, f(t))
		:returns: numpy array of coefficients for each cubic polynomial (4)*** for each interval, order corresponds to node ordering.
	"""

	# calculate coefficients
	z = cubic_spline_coefficients_bc(t, y, start_bc=(1,alpha), end_bc=(2,beta))

	return z
	
#
#	Testing
#

import matplotlib.pyplot as plt

if __name__ == '__main__':

	#
	# testing
	#

	# test function and its derivatives:
	f = lambda t: np.arctan(2*(t**2)) 						# arctan(2*(x^2))
	df = lambda t: (4*t)/( 4*(t**4) + 1) 					# 4x/( 4(x^4) +1 )
	ddf = lambda t: ( 4*( -12*(t**4) +1 ) )/( 4*(t**4) +1 ) # [4(-12(x^4) +1)]/( 4(x^4) +1 )

	# general parameters
	L = 3.0 # range (-L, L) for nodes
	nvals = [ 41, 21, 11, 5] # number of nodes

	for n in nvals:

		# generate nodes
		t = np.linspace(-L, L, num=n)
		y = f(t)

		# boundry conditions from HW problem 3
		alpha = df(t[0])
		beta = ddf(t[-1])

		# construct cubic spline class passing nodes and boundry conditions
		cs = CubicSpline(t, y, start_bc=(1,alpha), end_bc=(2,beta))

		# evaluation points to test cubic spline after construction
		xs = np.linspace(-L, L, 400)

		#
		# plot results
		#

		# generate comparison plots
		fig, axarr = plt.subplots(1,2)

		ftitle = '$f(z)$'
		stitle = cs.latex_error_title()

		# plot evaluations at nodes using scatter
		axarr[0].scatter(t,f(t),s=10,color='r', zorder=2)
		axarr[0].scatter(t,cs(t),s=5,color='b', zorder=2)

		# plot function graphs
		fplot = axarr[0].plot(xs, f(xs), 'r-', label=ftitle)
		splot = axarr[0].plot(xs, cs(xs) , 'b--', label=stitle)

		axarr[0].set_xlabel('x')
		axarr[0].set_ylabel('y')

		s_error_title = cs.latex_error_title()

		splot_error = axarr[1].plot(xs, abs(f(xs)-cs(xs)), 'b--', label=s_error_title)

		axarr[1].set_xlabel('x')
		axarr[1].set_ylabel('y')
		axarr[1].set_title('Problem 1: $E_{' + str(n) + '}$')


		leg = axarr[0].legend(loc='upper center', shadow=True)
		leg = axarr[1].legend(loc='upper center', shadow=True)


	plt.show()

