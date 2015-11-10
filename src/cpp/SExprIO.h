#ifndef S_EXPR_IO_H
#define S_EXPR_IO_H

//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: SExprIO
//
// PURPOSE:
//
// A few experimental routines for parsing of simple s-expressions.
//
// NOTES:
//
// This was designed to get polynomial reading from the python side of
// the code up and running quickly.  This is quite an ugly piece of code
// and it would be a good idea to replace it as soon as time pressures
// allow.
//
//----------------------------------------------------------------------

//standard headers
#include <iostream>
#include <vector>
#include <assert.h>

//library headers

//project headers

void seWriteOBrace(std::ostream& os);
void seWriteCBrace(std::ostream& os);
void seWriteStart(std::ostream& os, const std::string& s);
void seWriteEnd(std::ostream& os, const std::string& s);
void seWriteString(std::ostream& os, const std::string& s);

void seReadOBrace(std::istream& is);
void seReadCBrace(std::istream& is);
void seReadStart(std::istream& is, const std::string& s);
void seReadEnd(std::istream& is, const std::string& s);
std::string seReadString(std::istream& is);

//template declarations

template< typename T >
extern
std::vector< T > seWriteVector(std::ostream& os, const std::vector< T >& vec);

template< typename T >
extern
std::vector< T > seReadVector(std::istream& is, const size_t len);

#endif //S_EXPR_IO_H
