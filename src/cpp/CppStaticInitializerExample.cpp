//----------------------------------------------------------------------
//
// AUTHOR: Dr. Andrew David Burbanks, 2005.
// This software is Copyright (C) 2004-2008  Bristol University
// and is released under the GNU General Public License version 2.
//
// MODULE: CppStaticInitializerExample implementation
//
// NOTES:
//
// This was going to be an experient with how various compiler
// standards deal with initialisation of static member variables, but
// was later found not to be needed.
//
//----------------------------------------------------------------------

class dummy {
public:
  static const int var = -3;
};
