import unittest

import GMatElastoPlastic.Cartesian3d as GMat
import GMatTensor.Cartesian3d as tensor
import numpy as np


class Test_main(unittest.TestCase):
    """ """

    def test_elastic(self):

        shape = [2, 3]
        K = np.random.random(shape)
        G = np.random.random(shape)
        sigy0 = 1000 * np.ones(shape)
        H = np.random.random(shape)
        mat = GMat.LinearHardening2d(K, G, sigy0, H)

        gamma = np.random.random(shape)
        epsm = np.random.random(shape)

        Eps = np.zeros(shape + [3, 3])
        Eps[..., 0, 0] = epsm
        Eps[..., 1, 1] = epsm
        Eps[..., 2, 2] = epsm
        Eps[..., 0, 1] = gamma
        Eps[..., 1, 0] = gamma
        mat.Eps = Eps

        Sig = np.zeros(shape + [3, 3])
        Sig[..., 0, 0] = 3 * K * epsm
        Sig[..., 1, 1] = 3 * K * epsm
        Sig[..., 2, 2] = 3 * K * epsm
        Sig[..., 0, 1] = 2 * G * gamma
        Sig[..., 1, 0] = 2 * G * gamma

        self.assertTrue(np.allclose(GMat.Epseq(mat.Eps), 2 / np.sqrt(3) * gamma))
        self.assertTrue(np.allclose(GMat.Sigeq(mat.Sig), 2 * np.sqrt(3) * G * gamma))
        self.assertTrue(np.allclose(mat.Sig, Sig))
        self.assertTrue(np.allclose(tensor.A4_ddot_B2(mat.C, mat.Eps), Sig))
        self.assertTrue(np.allclose(mat.epsp, np.zeros(shape)))

    def test_elastoplastic(self):

        shape = [2, 3]
        K = np.random.random(shape)
        G = np.random.random(shape)
        sigy0 = np.random.random(shape)
        H = np.random.random(shape)

        G = np.ones(shape)
        sigy0 = np.ones(shape)
        H = np.ones(shape)

        mat = GMat.LinearHardening2d(K, G, sigy0, H)

        gamma = np.random.random(shape)
        epsm = np.random.random(shape)

        gamma = np.ones(shape)

        Eps = np.zeros(shape + [3, 3])
        Eps[..., 0, 0] = epsm
        Eps[..., 1, 1] = epsm
        Eps[..., 2, 2] = epsm
        Eps[..., 0, 1] = gamma
        Eps[..., 1, 0] = gamma
        mat.Eps = Eps

        Sig = np.zeros(shape + [3, 3])
        Sig[..., 0, 0] = 3 * K * epsm
        Sig[..., 1, 1] = 3 * K * epsm
        Sig[..., 2, 2] = 3 * K * epsm

        sigeq = 2 * np.sqrt(3) * G * gamma
        2 / np.sqrt(3) * gamma
        epsp = (sigeq - sigy0) / (3 * G + H)
        sigeq = np.where(sigeq > sigy0, sigy0 + H * epsp, sigeq)
        tau = sigeq / np.sqrt(3)

        Sig[..., 0, 1] = tau
        Sig[..., 1, 0] = tau

        self.assertTrue(np.allclose(mat.Sig, Sig))
        self.assertTrue(np.allclose(mat.epsp, epsp))

    def test_elastic_tangent(self):

        shape = [2, 3]
        Eps = np.random.random(shape + [3, 3])
        Eps = tensor.A4_ddot_B2(tensor.Array2d(shape).I4s, Eps)

        mat = GMat.LinearHardening2d(
            K=np.random.random(shape),
            G=np.random.random(shape),
            sigy0=1000 * np.ones(shape),
            H=np.random.random(shape),
        )

        mat.Eps = Eps
        self.assertTrue(np.allclose(tensor.A4_ddot_B2(mat.C, mat.Eps), mat.Sig))

    def test_plastic_tangent(self):

        shape = [2, 3]
        Eps = np.random.random(shape + [3, 3])
        Eps = tensor.A4_ddot_B2(tensor.Array2d(shape).I4s, Eps)

        mat = GMat.LinearHardening2d(
            K=np.random.random(shape),
            G=np.random.random(shape),
            sigy0=np.zeros(shape),
            H=np.random.random(shape),
        )

        mat.Eps = Eps
        self.assertTrue(np.allclose(tensor.A4_ddot_B2(mat.C, mat.Eps), mat.Sig))


if __name__ == "__main__":

    unittest.main()
