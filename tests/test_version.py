import unittest
import setuptools_scm

import GMatElastoPlastic as GMat


class Test_main(unittest.TestCase):
    """ """

    def test_version_dependencies(self):

        deps = GMat.version_dependencies()
        deps = [i.split("=")[0] for i in deps]

        self.assertTrue("gmatelastoplastic" in deps)
        self.assertTrue("gmatelastic" in deps)
        self.assertTrue("gmattensor" in deps)
        self.assertTrue("xtensor" in deps)
        self.assertTrue("xtensor-python" in deps)
        self.assertTrue("xtl" in deps)


if __name__ == "__main__":

    unittest.main()
