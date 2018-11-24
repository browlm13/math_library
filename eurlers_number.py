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


# check
for i in range(171):
	print("P(%s) = %s" % (i, p5(i)))

