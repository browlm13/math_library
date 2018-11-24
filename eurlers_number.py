"""


	eurlers number


	lim n->inf of P(n)

	where 

	P(n) = sum from i=0 to i=n of n!/n**i /i! / (n-i)!

"""


import math

def p(n):

	terms = []

	for i in range(n+1):

		t = math.factorial(n) / n**i / math.factorial(n-i) / math.factorial(i)
		terms += [t]

	return sum(terms)


# check
for i in range(171):
	print("P(%s) = %s" % (i, p(i)))
