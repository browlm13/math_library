# newton.py
#
# Daniel R. Reynolds
# Math5315 @ SMU
# Fall 2017

def newton2(Ffun, Jfun, x, maxit, rtol, atol, output):
    """
    Usage: x, its = newton(Ffun, Jfun, x, maxit, rtol, atol, output)

    This routine uses Newton's method to approximate a root of
    the nonlinear system of equations F(x)=0.  The iteration ceases 
    when the following condition is met:

       ||xnew - xold|| < atol + rtol*||xnew||

    inputs:   Ffun     nonlinear function name/handle
              Jfun     Jacobian function name/handle
              x        initial guess at solution
              maxit    maximum allowed number of iterations
              rtol    relative solution tolerance
              atol    absolute solution tolerance
              output   flag (true/false) to output iteration history/plot
    outputs:  x        approximate solution
              its      number of iterations taken
    """

    # imports
    import numpy as np
    import math
    import sys

    # check input arguments
    if (int(maxit) < 1):
        sys.stdout.write("newton: maxit = %i < 1. Resetting to 10\n" % (int(maxit)))
        maxit = 10
    if (rtol < 1e-15):
        sys.stdout.write("newton: Srtol = %g < %g. Resetting to %g\n" % (rtol, 1e-15, 1e-15))
        rtol = 1e-15
    if (atol < 0):
        sys.stdout.write("newton: Satol = %g < 0. Resetting to %g\n" % (atol, 1e-15))
        atol = 1e-15

    # evaluate function
    F = Ffun(x)
    f0norm = np.linalg.norm(F)
        
    # perform iteration
    for its in range(1,maxit+1):

        # evaluate derivative
        J = Jfun(x)

        # compute Newton update, new guess at solution, new residual
        if (np.isscalar(x)):       # if problem is scalar-valued
            h = F/J
        else:                      # if problem is vector-valued
            h = np.linalg.solve(J, F)
        x = x - h
        F = Ffun(x)

        # check for convergence and output diagnostics
        hnorm = np.linalg.norm(h)
        xnorm = np.linalg.norm(x)
        fnorm = np.linalg.norm(F)
        if (output):
            sys.stdout.write("Method1:   iter %3i, \t||h|| = %g, \thtol = %g, \t||f|| = %g\n" % 
                             (its, hnorm, atol+rtol*xnorm, fnorm))

        if (hnorm < atol + rtol*xnorm):
            break

    return [x, its]

