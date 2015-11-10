//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: Utility implementation
//
//----------------------------------------------------------------------

//standard headers
#include <assert.h>

//library headers

//project headers
#include "Types.h"
#include "Utility.h"

//efficiently raise a value to a power
Coeff coeffPow(const Coeff c, const Power p) {
  assert(isNonnegativePower(p));
  Coeff result(CoeffOne);
  if (p > 0) {
    Power mask(1);
    while (mask <= p) {
      mask <<= 1;
    }
    //position of p's most-significant bit
    mask >>= 1;
    while (mask > 0) {
      result *= result;
      if (p & mask) {
	result *= c;
      }
      mask >>= 1;
    }
  }
  return result;
}

// efficiently raise any integer to a non-negative power
int intPow(const int i, const size_t p) {
  assert(p >= 0);
  size_t result = 1;
  if (p > 0) {
    size_t mask(1);
    while (mask <= p) {
      mask <<= 1;
    }
    //position of p's most-significant bit
    mask >>= 1;
    while (mask > 0) {
      result *= result;
      if (p & mask) {
	result *= i;
      }
      mask >>= 1;
    }
  }
  return result;
}

size_t factorial(const int i) {
  //(1)(2)(3)...(i-2)(i-1)(i)
  if (i < 0)
    throw UtilityNegativeArgumentError();
  size_t result = 1;
  for (size_t j = 2; j <= i; ++j)
    result *= j;
  //worry a little about overflow
  if (!(i <= result))
    throw TypesOverflowError();
  return result;
}

size_t factorial(const int i, const int j) {
  //utility for use in binomial
  //
  //(i)(i+1)(i+2)...(j-2)(j-1)(j)
  if (i < 0)
    throw UtilityNegativeArgumentError();
  if (j < 0)
    throw UtilityNegativeArgumentError();
  //we allow j<i for convenience in binomial
  size_t result = 1;
  for (size_t k = i; k <= j; ++k)
    result *= k;
  //worry a little about overflow
  if (i <= j)
    if (!((i <= result) && (j <= result)))
      throw TypesOverflowError();
  return result;
}

size_t binomial(const int n, const int k) {
  // This algorithm appeared in "Binomial Coefficient Computation:
  // Recursion or Iteration?" Yannis Manolopoulos, Data Engineering
  // Laboratory, Department of Informatics, Aristotle University
  // Thessaloniki, 54006 Greece <manolopo@delab.csd.auth.gr>.
  //
  // 1. It is iterative, thus avoiding time overhead for function
  // calls and space overhead for stacks, 2. It has optimal
  // complexity, that is O(min(k,n-k)), 3. It is robust, as it
  // performs multiplications and division alternatively, thus
  // avoiding data type overflows.
  if (n < 0)
    throw UtilityNegativeArgumentError();
  if (k < 0)
    throw UtilityNegativeArgumentError();
  if (k > n)
    throw UtilityOutOfRangeError();
  size_t t = 1;
  if (k < n-k)
    for (size_t i = n; i > n-k; --i) {
      t *= i;
      t /= (n-i+1);
    }
  else
    for (size_t i = n; i > k; --i) {
      t *= i;
      t /= (n-i+1);
    }
  return t;
}

/*
size_t binomialWithOverflow(const size_t n, const size_t k) {
  //this implementation will overflow for large n, k.
  //we should use Reals for internal computation then?
  assert(0 <= n);
  assert((0 <= k) && (k <= n));
  const size_t min = (k <= n-k) ? k : (n - k);
  const size_t max = (k <= n-k) ? (n - k) : k;
  //std::cerr << factorial(max+1, n) << std::endl;
  //std::cerr << factorial(min) << std::endl;
  return Real(factorial(max+1, n))/Real(factorial(min));
}
*/
