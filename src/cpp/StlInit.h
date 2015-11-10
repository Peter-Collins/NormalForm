//Author: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.

#ifndef STL_INIT_H
#define STL_INIT_H

//system
#include <vector>
#include <assert.h>

//library

//project
// gcc 4 needs forward declarations 
template < typename T >
class InitVector;
template < typename T >
InitVector< T > init(std::vector< T >& v); //to construct

template < typename T >
class InitVector {

  // It is the purpose of this little bit of code to allow us to fill
  // values into STL vectors in a less painful way.  We can write code
  // such as:
  //
  //   std::vector< double > v;
  //   init(v) = 1.0, 2.2, -7.9, 56.3;
  //
  // instead of:
  //
  //   std::vector< double > v;
  //   v.push_back(1.0);
  //   v.push_back(2.2);
  //   v.push_back(-7.9);
  //   v.push_back(56.3);
  //
  // we can even extend such a vector, via code like:
  //
  //   init(v) += 77.9, 32.1;
  //
  // This is accomplished by overloading the comma, i.e., operator,(),
  // the assignment, and the in-place addition operators of a class
  // InitVector.
  //
  // InitVector instances are created only via the friend function
  // init, and such entities as copy constructors are made private to
  // effectively disable them (- Keep It Simple, Stupid!).
  //
  // This method could probably be used for other STL containers,
  // although perhaps we would need 2 inner classes to cope with
  // assigning (key, value)-pairs for maps?

 private:

  InitVector< T >(void); //default - private, use init
  InitVector< T >(const InitVector< T >& other); //copy - private! disabled
  InitVector< T > operator=(const InitVector< T >& other); //assign - private!
  InitVector< T > operator=(InitVector< T >& other); //assign - private!
  explicit InitVector< T >(std::vector< T >& v); //from existing vector

 protected:

  std::vector< T >& vec_; //must be a reference; we mutate the original!
  size_t count_;

 public:

  ~InitVector< T >(void); //destruct
  InitVector< T >& operator=(const T& x);
  InitVector< T >& operator+=(const T& x);
  InitVector< T >& operator,(const T& x);
  friend InitVector< T > init< T >(std::vector< T >& v); //to construct
};

template < typename T >
InitVector< T > init(std::vector< T >& v) {
  //creation method
  return InitVector< T >(v);
}

template < typename T >
InitVector< T >::InitVector(void)
  : vec_(std::vector< T >()), count_(0) {
  //default constuction; must make a new anonymous container - useful?
  assert(vec_.size() == count_);
}

template < typename T >
InitVector< T >::~InitVector(void) {
  //destruction; nothing to do here
}

template < typename T >
InitVector< T >::InitVector(std::vector< T >& v)
  : vec_(v), count_(v.size()) {
  //construction from a container - important that we record its size!
  assert(vec_.size() == count_);
}

template < typename T >
InitVector< T >& InitVector< T >::operator=(const T& x) {
  //initialize the container from scratch; empty it if necessary
  assert(vec_.size() == count_);
  vec_.clear();
  count_ = 0;
  assert(vec_.size() == 0);
  vec_.push_back(x);
  assert(vec_[count_] == x);
  ++count_;
  assert(vec_.size() == count_);
  return *this;
}

template < typename T >
InitVector< T >& InitVector< T >::operator+=(const T& x) {
  //append to an existing container
  assert(vec_.size() == count_);
  vec_.push_back(x);
  assert(vec_[count_] == x);
  ++count_;
  assert(vec_.size() == count_);
  return *this;
}

template < typename T >
InitVector< T >& InitVector< T >::operator,(const T& x) {
  //append to the container (within a comma-separated list of items)
  assert(vec_.size() == count_);
  vec_.push_back(x);
  assert(vec_[count_] == x);
  ++count_;
  assert(vec_.size() == count_);
  return *this;
}

#endif //STL_INIT_H
