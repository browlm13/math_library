"""


	eurlers number


	lim n->inf of P(n)

	where 

	P(n) = sum from i=0 to i=n of n!/n**i /i! / (n-i)!

"""


import math
import numpy as np


def p1(n):

	terms = []

	for i in range(n+1):

		t = math.factorial(n) / n**i / math.factorial(n-i) / math.factorial(i)
		terms += [t]

	return sum(terms)

def p2(n):

	d1 = [1/n**i for i in range(n+1)]
	d2 = [1 / math.factorial(n-i) / math.factorial(i) for i in range(n+1)]

	d = np.dot(d1, d2)
	return math.factorial(n) * d


def p3(n):

	d0 = 1/np.array([n**i for i in range(n+1)])

	d1a = 1/np.array([math.factorial(i) for i in range(n+1)])
	d1b = np.flip(d1a)
	d1 = np.multiply(d1a, d1b)

	d = np.dot(d0, d1)

	return math.factorial(n) * d

from scipy.special import comb
def p4(n):

	t1 = 1/np.array([n**i for i in range(n+1)])
	t2 = np.array([comb(n,i) for i in range(n+1)])
	t = np.dot(t1,t2)

	return t

def p5(n):

	ts = np.array([comb(n,i)/n**i for i in range(n+1)])
	return np.sum(ts)

def p6(n):

	t1 = np.array([comb(n,i) for i in range(n+1)])
	t2 = np.array([n**i for i in range(n+1)])

	return np.sum(np.divide(t1,t2))

# newton's series Expansion for e
def p7(n):

	ts = np.array([1/math.factorial(i) for i in range(n+1)])

	return np.sum(ts)

# brothers' formulae
def brothers(n):
	ts = np.array([(2*i + 2)/math.factorial(2*i + 1) for i in range(n+1) ])
	return np.sum(ts)


def p8(n):
	ts = np.array([(i + 2)/math.factorial(i + 1) for i in range(0,n+1,2) ])
	return np.sum(ts)


# continue combining terms in newtons method (combining 3 terms)
def p9(n):

	ts = np.array([(i**2 + 4*i + 5)/math.factorial(i + 2) for i in range(0,n+1,3) ])
	return np.sum(ts)

# check
#for i in range(171):
#	print("P(%s) = %s" % (i, p9(i)))

eTrue = 2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260
for i in range(16):

	err8 = abs(p8(i)- eTrue)
	err9 = abs(p9(i) - eTrue)

	print("\ni = %s,  err8 = %s,  err9 = %s" % (i, err8, err9)) 
	print("i = %s,  err8 - err9 = %s" % (i, err8-err9))

