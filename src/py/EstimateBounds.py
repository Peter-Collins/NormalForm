# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

def fact(n):
    """

    Factorial of n, with safety checks.

    """
    assert isinstance(n, long) or isinstance(n, int)
    assert n>=0
    return _fact(n)

def _fact(n):
    """

    Factorial of n, without any checks.

    """
    prod = 1
    for k in xrange(1, n+1):
        prod *= k
    return prod

def binom(n, k):
    """

    Binomial coefficient by the inefficient formula.

    """
    assert isinstance(n, long) or isinstance(n, int)
    assert isinstance(k, long) or isinstance(k, int)
    assert n>=0
    assert k>=0
    assert k<=n
    return _fact(n)/(_fact(k)*_fact(n-k))

def p(k, n):
    """
    
    Number of (order-independent) partitions of natural number n,
    using addends of size at least k.

    """
    assert isinstance(k, long) or isinstance(k, int)
    assert isinstance(n, long) or isinstance(n, int)
    assert k >= 0
    assert n >= 0
    return _p(k, n)

def _p(k, n):
    """

    Version without checks.
    
    """
    if k > n:
        return 0
    if k == n:
        return 1
    return _p(k+1, n) + _p(k, n-k)

def parts(n):
    """

    Order-independent partitions of natural number n.

    """
    return p(1, n)

if __name__ == '__main__':
    #factorials
    assert fact(0) == 1
    assert fact(1) == 1
    assert fact(2) == 2
    for i in xrange(1, 100):
        assert fact(i) == i*fact(i-1)

    #partitions (order-independent) using addends of at least k
    assert p(1,4) == 5
    assert p(2,8) == 7
    assert p(3,12) == 9
    assert p(4,16) == 11
    assert p(5,20) == 13
    assert p(6,24) == 16

    #partitions (order-independent)
    assert parts(1) == 1 #1
    assert parts(2) == 2 #2, 1+1
    assert parts(3) == 3 #3, 2+1, 1+1+1
    assert parts(4) == 5 #4, 3+1, 2+2, 2+1+1, 1+1+1+1
    assert parts(5) == 7
    assert parts(6) == 11
    assert parts(7) == 15
    assert parts(8) == 22
    assert parts(9) == 30
    assert parts(10) == 42
    #assert parts(100) == 190569292
    #assert parts(1000) == 24061467864032622473692149727991
    #assert parts(10000) == 36167251325636293988820471890953695495016030339315650422081868605887952568754066420592310556052906916435144

