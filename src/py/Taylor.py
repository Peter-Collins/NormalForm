"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: Taylor

PURPOSE:

Experimental code to implement cached lazy-evaluated Taylor series.

NOTES:

The code here is for the 1-variable case!

"""

from Powers import Powers
from Polynomial import Polynomial
from Utility import factorial

class Lazy:

    def _ith_term(self, i):
        raise NotImplementedError

    def __getitem__(self, i):
        x = self._ith_term(i)
        for m, c in x.powers_and_coefficients():
            assert m.degree() == i
        return x

    def __add__(self, other):
        assert isinstance(other, Lazy)
        return Sum(self, other)

    def __mul__(self, other):
        assert isinstance(other, Lazy)
        return Product(self, other)

class Sum(Lazy):

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def _ith_term(self, i):
        return self.a[i] + self.b[i]

class Product(Lazy):

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def _ith_term(self, i):
        res = self.a[0]*self.b[i]
        for j in xrange(1, i+1):
            res += self.a[j]*self.b[i-j]
        return res

class Cached(Lazy):

    def __init__(self, obj):
        self.obj = obj
        self.values = {}

    def __getitem__(self, i):
        if self.values.has_key(i):
            return self.values[i]
        else:
            x = self.obj[i]
            self.values[i] = x
            return x

class Sine(Lazy):

    def __init__(self, n_vars, index):
        self.n_vars = n_vars
        self.index = index

    def _ith_term(self, i):
        if i%2:
            powers = [0,]*self.n_vars
            powers[self.index] = i
            if i%4 == 3:
                sign = -1
            else:
                sign = 1
            return (float(sign)/float(factorial(i)))*Polynomial.Monomial(tuple(powers))
        else:
            return Polynomial(self.n_vars)

class Cosine(Lazy):

    def __init__(self, n_vars, index):
        self.n_vars = n_vars
        self.index = index

    def _ith_term(self, i):
        if i%2:
            return Polynomial(self.n_vars)
        else:
            powers = [0,]*self.n_vars
            powers[self.index] = i
            if i%4 == 2:
                sign = -1
            else:
                sign = 1
            return (float(sign)/float(factorial(i)))*Polynomial.Monomial(tuple(powers))

tanh_numerators = [1,
                   1,
                   2,
                   17,
                   62,
                   1382,
                   21844,
                   929569,
                   6404582,
                   443861162,
                   18888466084,
                   113927491862,
                   58870668456604,
                   8374643517010684,
                   689005380505609448,
                   129848163681107301953,
                   1736640792209901647222,
                   418781231495293038913922]

#http://www.research.att.com/~njas/sequences/
#Bernouilli numerators

a027641 = [1,-1,1,0,-1,0,1,0,-1,0,5,0,-691,0,7,0,-3617,0,43867,0,
           -174611,0,854513,0,-236364091,0,8553103,0,-23749461029,0,
           8615841276005,0,-7709321041217,0,2577687858367,0,
           -26315271553053477373,0,2929993913841559,0,
           -261082718496449122051]

#http://www.research.att.com/~njas/sequences/
#Bernouilli denominators

a027642 = [1,2,6,1,30,1,42,1,30,1,66,1,2730,1,6,1,510,1,798,1,330,1,
           138,1,2730,1,6,1,870,1,14322,1,510,1,6,1,1919190,1,6,1,
           13530,1,1806,1,690,1,282,1,46410,1,66,1,1590,1,798,1,870,1,
           354,1,56786730,1]
           
def bernoulli(i):
    return float(a027641[i])/float(a027642[i])

class Tanh(Lazy):

    def __init__(self, n_vars, index):
        self.n_vars = n_vars
        self.index = index

    def _ith_coeff(self, i):
        if (i%2) == 0:
            return 0.0
        else:
            b = float(bernoulli(i+1))
            p = float(2 ** (i+1))
            f = float(factorial(i+1))
            return b*p*(p-1)/f

    def _ith_term(self, i):
        if i > 0:
            sign = 1-(i%2)*2
            powers = [0,]*self.n_vars
            powers[self.index] = i
            m = Polynomial.Monomial(tuple(powers))
            return self._ith_coeff(i)*m
        else:
            return Polynomial(self.n_vars)
