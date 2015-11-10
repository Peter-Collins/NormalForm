# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest

from ShelvedDict import ShelvedDict

class EmptyDict(unittest.TestCase):
    def test_empty(self):
        s = ShelvedDict()
        d = {}
        self.assert_(s==d, 'Empty shelf != empty dict.')
        self.assert_(d==s, 'Empty shelf != empty dict.')

class DictTestExamples(unittest.TestCase):

    def setUp(self):
        self.examples = [{},
                         {1:1, 2:2, 3:3},
                         {1:'a', 2:'b', 3:'c'}]

    def tearDown(self):
        pass

    def _cond(self, a, b):
        raise NotImplementedError

    def test_examples(self):
        for d in self.examples:
            s = ShelvedDict()
            for key, item in d.iteritems():
                s[key] = item
            self.assert_(self._cond(s, d))

class EqualsLikeDict(DictTestExamples):
    def _cond(self, a, b):
        return (a == b) and (b == a)

class NotIsLikeDict(DictTestExamples):
    def _cond(self, a, b):
        return (not (a is b)) and (not (b is a))

class EqualsSelf(DictTestExamples):
    def _cond(self, a, b):
        return (a == a) and (b == b)

def suite():
    suites = []
    suites.append(unittest.makeSuite(EmptyDict))
    suites.append(unittest.makeSuite(EqualsLikeDict))
    suites.append(unittest.makeSuite(NotIsLikeDict))
    suites.append(unittest.makeSuite(EqualsSelf))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')

