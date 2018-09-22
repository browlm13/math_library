#!/usr/bin/env python3

"""
Script to compare the performance of Newton’s method 
and the Picard iteration on the root-finding problem

Laurence Brown
Math5315 @ SMU
Fall 2018
"""

import numpy as np
import time
from newton import *
from picard import *

# general parameters for all tests
maxit = 50
Satols = np.array([1e-6])
Srtols = np.array([1e-6])
Ratols = np.array([1e-10])
Rrtols = np.array([1e-10])
ntols = len(Satols)


# Testing the performance of Newton’s method and the Picard iteration
# on the root-finding problem defined in question 3.

print("\nProblem 3:")


x0 = np.array([1,0,0])                     # initial guess
x1 = np.array([0.95,0,0.01])               # second initial guess
x2 = np.array([1.05,0.01,0])               # third initial guess

"""
x0 = np.random.uniform(-100,100,3)                     # initial guess
x1 = np.random.uniform(-100,100,3)               # second initial guess
x2 = np.random.uniform(-100,100,3)               # third initial guess
"""

initial_guesses = [x0,x1,x2]

def to_homogeneous(x):
  """ 
    converts numpy array x to homogeneous coordinates 
    by adding an extra row to the vector with a value 1.

          use: h = to_homogeneous(x)

    :param x: numpy array, vector to convert to homogeneous coordinates.
    :returns: homogeneous numpy array
  """
  h = np.ones(shape=(x.shape[0]+1))
  h[:x.shape[0]] = x[:]
  return h

def compute_J(x, *Fs):
  """
    Compute the jacobian matrix for a given vector x, 
    using the passed list of matrices to compute the columns of J.

      use: J = compute_J(x, Fx1, Fx2, Fx3)

    :param x: vector of length n.
    :param Fs: list of matrices of dimension nxn+1, where column j of J is: Jij = Fij dot homogeneous(x).
    :returns: Jacobian Matrix
  """
  Jis = []
  for Fxi in Fs: 
    Jis.append(np.dot(Fxi,to_homogeneous(x)))
  return np.vstack(tuple(Jis)).T

def J(x):
  """
   Compute Jacobian Matrix for Question 4 given an x vector.
   :param x: vector of length n.
   :returns: Jacobian Matrix (nxn).
  """

  # specific Jacobian matrix components for dfi/dx1 'first column in J'
  # Fx1 dot x = Ji0
  Fx1 = np.array([ [0, 0, 0,  1.004],
                   [0, 0, 0, -0.004],
                   [0, 0, 0,      0] ])

  # Fx2 dot x = Ji1
  Fx2 = np.array([ [0,  0, -1000, 0],
                   [0,  60, 1000, 1],
                   [0, -60,    0, 0] ])

  # Fx3 dot x = Ji2
  Fx3 = np.array([ [0,  -1000, 0, 0],
                   [0,   1000, 0, 0],
                   [0,      0, 0, 1] ])


  return compute_J(x, Fx1, Fx2, Fx3)

def f(x):
  """
    Compute the function evaluted at the vector x.
      use: F = f(x)

    :param x: vector x of length 3.
    :returns: vector of length 3, or function from question 4 evaluted at x.
  """
  f0 = lambda x: x[0]*(1 + 0.004) + x[1]*x[2]*(-1000) -1
  f1 = lambda x: x[0]*(-0.004) + x[1] + x[1]**2*(30) +  x[1]*x[2]*1000
  f2 = lambda x: x[1]**2*(-30) + x[2]

  return np.array([f0(x), f1(x), f2(x)])

for x_initial in initial_guesses:
  input("Press Enter to continue...")
  print("Initial guess x: %s\n\n" % x_initial)

  # call solvers on this problem
  for i in range(ntols):
    print("\n  Srtol = ", Srtols[i], ", Satol = ", Satols[i],", Rrtol = ",
          Rrtols[i], ", Ratol = ", Ratols[i], ":")

    print("  Newton:")
    tstart = time.time()
    x, its = newton(f, J, x_initial, maxit, Srtols[i], Satols[i], Rrtols[i], Ratols[i], True)
    ttot = time.time() - tstart
    print("   final solution = ", x, ", residual = ", np.abs(f(x)), ", time = ", ttot, ", its = ", its)

    print("  Picard:")
    tstart = time.time()
    x, its = picard(f, J(x_initial), x_initial, maxit, Srtols[i], Satols[i], Rrtols[i], Ratols[i], True)
    ttot = time.time() - tstart
    print("   final solution = ", x, ", residual = ", np.abs(f(x)), ", time = ", ttot, ", its = ", its)
