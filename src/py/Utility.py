"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE:

PURPOSE:

Miscellaneous utility functions, such as binomial coefficients.

NOTES:

"""

def pb(a, b):
    return a.poisson_bracket(b)

def binomial(n, k):
    """

    ADB: This is an efficient iterative implementation, from the paper
    [Manolopoulos] "Binomial Coefficient Computation: Recursion or
    Iteration?", SIGCSE Bulletin, Vol. 34, No. 4, 2002 December.  It
    is compared, in the unit tests, against a reference
    implementation, against the binomial identities, and against
    sample values.

    """
    if n < 0: raise ValueError
    if k < 0 or k > n: raise ValueError
    t = 1
    if k < (n-k):
        for i in xrange(n, (n-k), -1):
            t *= i
            t /= (n-i+1)
    else:
        for i in xrange(n, k, -1):
            t *= i
            t /= (n-i+1)
    return t

def factorial(x):
    r = 1
    for i in xrange(2, x+1):
        r *= i
    return r
