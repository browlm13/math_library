#!/usr/bin/env python3

"""
File name: hermite.py

    Python Version: 3.6

    	Cubic Spline Interpolation.

    	Boundary conditions:

    		S'(t_0) = alpha
    		S''(t_n) = beta

    		two equations for boundry conditions:

    		1.) 2*h_0 * (z_0) + h_0 * (z_1) = b_0 -6 * alpha

				first row:

					[2*h_0,  h_0, 0, ..., 0 ] * [z_0, z_1, ..., z_n].T = [ (b_0 -6*alpha), v_1, ..., v_n-1, beta]

    		2.) (z_n) = beta

				last row:

					[0, ..., 1] * [z_0, ..., z_n].T  = [ (b_0 -6*alpha), v_1, ..., v_n-1, beta]
    		

LJ Brown
Fall 2018
"""

__filename__ = "hermite.py"
__author__ = "L.J. Brown"

# internal libraries
import logging
from functools import reduce

# external libraries
import numpy as np

# initilize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def cubic_spline_coefficients(t, y, alpha, beta):

	n = len(t) 

	# build tridiagonal matrix initially filled with zeros
	A = np.zeros(shape=(n, n))

	# build vector r from equation Az = r
	r = np.zeros(shape=(n,))

	# h_i = t_i+1 - t_i
	h = lambda i: t[i+1] - t[i] 

	# d_i = 2(h_i-1 + h_i)
	d = lambda i: 2*(h(i-1) + h(i))

	# b_i = (6/h_i) * (y_i+1 - y_i)
	b = lambda i: (6/h(i))*(y[i+1] - y[i])

	# v_i = b_i - b_i-1
	v = lambda i: b(i) - b(i-1)

	# fill in equations for intererior nodes rows 1->n-1
	for i in range(1,n-1):
		A[i, i-1] = h(i-1)
		A[i, i] = d(i)
		A[i, i+1] = h(i)

		r[i] = v(i)

	# fill in boundry equations

	# S'(t_0) = alpha
	# 2*h_0 * (z_0) + h_0 * (z_1) = b_0 -6 * alpha
	A[0,0], A[0,1] = 2*h(0), h(0)
	r[0] = b(0) -6 * alpha

	# S''(t_n) = beta
	# (z_n) = beta
	A[n-1,n-1] = 1
	r[n-1] = beta

	# solve system for coefficients, z
	z = np.linalg.solve(A, r)

	print(A)
	print(r)
	print(z)

	return z

def cubic_spline_evaluate(t, y, z, x):

	# h_i = t_i+1 - t_i
	h = lambda i: t[i+1] - t[i] 

	n = len(t)-1

	S = lambda i,x: (z[i]/(6*h(i)))*(t[i+1] - x)**3 + (z[i+1]/(6*h(i)))*(x - t[i])**3 + ((y[i+1]/h(i)) - (z[i+1]*h(i)/6))*(x - t[i]) + ((y[i]/h(i)) - (z[i]*h(i)/6))*(t[i+1] - x)
	for i in range(1, n-2):
		if x <= t[i]:
			return S(i,x)
	return S(n-1, x)


import matplotlib.pyplot as plt

if __name__ == '__main__':

	# test
	f = lambda t: np.arctan(2*(t**2))

		# general parameters
	L = -3.0
	n = 5
	t = np.linspace(-L, L, 5)
	y = f(t)
	alpha = 1
	beta = 3
	
	z = cubic_spline_coefficients(t,y,alpha,beta)

	interpolation_ys = []
	for x in t:
		interpolation_ys.append(cubic_spline_evaluate(t,y,z,x))

	xs = np.linspace(-L, L, 201)     # evaluation points
	answers = []
	for x in xs:
		answers.append(cubic_spline_evaluate(t,y,z,x))

	# generate comparison plots
	fig, axarr = plt.subplots(1,2)


	ftitle = '$f(z)$'
	stitle = '$S(z)$'

	axarr[0].scatter(t,y,s=10,color='r', zorder=2)
	axarr[0].scatter(t,interpolation_ys,s=10,color='g', zorder=2)

	fplot = axarr[0].plot(xs, f(xs), 'r-', label=ftitle)
	lplot = axarr[0].plot(xs, answers, 'g--', label=stitle)

	axarr[0].set_xlabel('x')
	axarr[0].set_ylabel('y')

	leg = axarr[0].legend(loc='upper center', shadow=True)


	plt.show()

