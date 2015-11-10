"""

AUTHOR: Dr. Andrew David Burbanks, 2005.
This software is Copyright (C) 2004-2008  Bristol University
and is released under the GNU General Public License version 2.

MODULE: SExprIO

PURPOSE:

Basic routines for writing and reading simple S-expressions.

NOTES:

These are rather primitive, but they work for the subset of
S-expressions that I have decided to use.

"""

def strip_begin(s, head=None):
    #strip beginning brace and (optional) keyword from sexp
    s = s.strip()
    if head:
        assert s.startswith('('+head), s
        r = s[len('('+head):]
    else:
        assert s.startswith('('), s
        r = s[1:]
    return r.strip()
def strip_end(s, head=None):
    #strip end brace and (optional) keyword from sexp
    s = s.strip()
    if s.endswith(')'):
        r = s[:-1]
    else:
        assert s.endswith(');;'+head), s
        r = s[:-len(');;'+head)]
    return r.strip()
def strip_braces(s, head=None):
    #strip braces and keyword enclosing sexp
    return strip_end(strip_begin(s, head=head), head=head)
def elts(s):
    #get space-separated elements from a string
    return s.strip().split()

