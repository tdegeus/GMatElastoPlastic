import unittest

import GMatElastoPlastic.Cartesian3d as GMat
import h5py
import numpy as np


class Test(unittest.TestCase):
    def test_main(self):

        with h5py.File("Cartesian3d_random.hdf5", "r") as data:

            mat = GMat.LinearHardening2d(
                K=data["K"][...],
                G=data["G"][...],
                sigy0=data["sigy0"][...],
                H=data["H"][...],
            )

            for i in range(20):

                mat.increment()
                mat.Eps = data[f"/data/{i:d}/Eps"][...]
                self.assertTrue(np.allclose(mat.Sig, data[f"/data/{i:d}/Stress"][...]))
                self.assertTrue(np.allclose(mat.C, data[f"/data/{i:d}/Tangent"][...]))
                self.assertTrue(np.allclose(mat.epsp, data[f"/data/{i:d}/epsp"][...]))


if __name__ == "__main__":

    unittest.main()
