def lagrange(x, y, z):
    """ Usage: p = lagrange(x,y,z)                               """
    """                                                          """
    """ This routine evaluates the Lagrange interpolating        """
    """ polynomial p(z) defined in the Lagrange basis by         """
    """   p(z) = y[0]*l1(z) + y[1]*l2(z) + ... + y[n-1]*ln(z).   """
    """                                                          """
    """ Inputs:   x  nodal locations [array of length n],        """
    """              (we assume that x(i) ~= x(j) for i != j)    """
    """           y  data values [array of length n]             """
    """           z  evaluation point                            """
    """ Outputs:  p  p(z)                                        """
    """                                                          """
    """ Daniel R. Reynolds                                       """
    """ Math5315 @ SMU                                           """
    """ Fall 2018                                                """

    # check inputs
    if (x.size != y.size):
        raise ValueError("lagrange error: (x,y) have different sizes")

    n = x.size

    # initialize output
    p = 0.0

    # iterate over Lagrange basis functions
    for k in range(n):
   
        # initialize l (the kth Lagrange basis function)
        l = 1.0
   
        # iterate over data to construct l(z)
        for j in range(n):
            # exclude the k-th data point
            if (j != k):
                l *= (z-x[j]) / (x[k]-x[j])
   
        # add contribution from this basis function (and data value) into p
        p += y[k]*l
   
    return p
