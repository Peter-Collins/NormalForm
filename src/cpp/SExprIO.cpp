//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: SExprIO implementation.
//
//----------------------------------------------------------------------

//standard headers
#include <iostream>
#include <string>
#include <assert.h>

//library headers

//project headers
#include "Types.h"
#include "SExprIO.h"

//template function implementations
template< typename T >
std::vector< T > seWriteVector(std::ostream& os, const std::vector< T >& vec) {
  seWriteOBrace(os);
  for (size_t i = 0; i < vec.size(); ++i) {
    if (i >= 1)
      os << " ";
    os << vec[i];
  }
  seWriteCBrace(os);
}

template< typename T >
std::vector< T > seReadVector(std::istream& is, const size_t len) {
  seReadOBrace(is);
  T item;
  std::vector< T > result;
  for (size_t i = 0; i < len; ++i) {
    assert(is);
    is >> item;
    result.push_back(item);
  }
  seReadCBrace(is);
  return result;
}

void seWriteOBrace(std::ostream& os) {
  os << "(";
}

void seWriteCBrace(std::ostream& os) {
  os << ")";
}

void seWriteStart(std::ostream& os, const std::string& s) {
  seWriteOBrace(os);
  os << s;
}

void seWriteEnd(std::ostream& os, const std::string& s) {
  seWriteCBrace(os);
}

void seWriteString(std::ostream& os, const std::string& s) {
  os << "\"" << s << "\"";
}

void seReadOBrace(std::istream& is) {
  is >> std::ws;
  assert(is);
  char c;
  is.get(c);
  assert(c == '(');
}

void seReadCBrace(std::istream& is) {
  is >> std::ws;
  assert(is);
  char c;
  is.get(c);
  assert(c == ')');
}

void seReadStart(std::istream& is, const std::string& s) {
  seReadOBrace(is);
  assert(is);
  std::string buffer;
  is >> buffer;
  assert(buffer == s);
}

void seReadEnd(std::istream& is, const std::string& s) {
  //for now, we do not allow any ending comments!
  seReadCBrace(is);
}

std::string seReadString(std::istream& is) {
  is >> std::ws;
  assert(is);
  char c;
  is.get(c);
  assert(c == '\"');
  //
  char buffer[256];
  is.get(buffer, 256, '\"');
  //
  assert(is);
  is.get(c);
  assert(c == '\"');
  //
  return buffer;
}
