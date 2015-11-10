# This software is Copyright (C) 2004-2008  Bristol University
# and is released under the GNU General Public License version 2.

import unittest

class ClassAttributesAsDefaults(unittest.TestCase):

    def setUp(self):
        """

        Note: we must encapsulate the class definitions here, since
        our test code modifies class attributes as well as instance
        attributes; thus, if we were to define our classes once, we
        run into the problem of those definitions changing between
        test methods.

        """
        class ClassWithDefaultAttribute:
            _hello = 'class hello'
        self.ClassWithDefaultAttribute = ClassWithDefaultAttribute
        class ClassWithOptionalOverride:
            _hello = 'class hello'
            def __init__(self, hello=None):
                if hello:
                    self._hello = hello
        self.ClassWithOptionalOverride = ClassWithOptionalOverride
        class ClassWithOverrideMethod:
            _hello = 'class hello'
            def change_self(self):
                self._hello = 'new self hello'
            def change_class(self):
                ClassWithOverrideMethod._hello = 'new class hello'
        self.ClassWithOverrideMethod = ClassWithOverrideMethod

    def tearDown(self):
        pass

    def test_default(self):
        ins = self.ClassWithDefaultAttribute()
        self.assert_(ins._hello == 'class hello')

    def test_default_in_override(self):
        ins = self.ClassWithOptionalOverride()
        self.assert_(ins._hello == 'class hello')

    def test_overridden(self):
        ins = self.ClassWithOptionalOverride('instance hello')
        self.assert_(ins._hello == 'instance hello')

    def test_setting_instance_attributes_shadows_class(self):
        a = self.ClassWithOptionalOverride('a')
        a._hello = 'c'
        self.assert_(a._hello == 'c')
        b = self.ClassWithOptionalOverride()
        self.assert_(b._hello == 'class hello')

    def test_override_default(self):
        a = self.ClassWithOverrideMethod()
        self.assert_(a._hello == 'class hello')

    def test_override_self_method_overrides(self):
        a = self.ClassWithOverrideMethod()
        a.change_self()
        self.assert_(a._hello == 'new self hello')

    def test_override_class_method_overrides_default(self):
        a = self.ClassWithOverrideMethod()
        a.change_class()
        self.assert_(a._hello == 'new class hello')

    def test_override_self_shadows_class(self):
        a = self.ClassWithOverrideMethod()
        a.change_self()
        a.change_class()
        self.assert_(a._hello == 'new self hello')

    def test_override_self_method_does_not_override_other(self):
        a = self.ClassWithOverrideMethod()
        a.change_self()
        b = self.ClassWithOverrideMethod()
        self.assert_(a._hello == 'new self hello')
        self.assert_(b._hello == 'class hello')

def suite():
    suites = []
    suites.append(unittest.makeSuite(ClassAttributesAsDefaults))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
    
