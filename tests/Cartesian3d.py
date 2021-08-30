import unittest
import numpy as np
import GMatElastoPlastic.Cartesian3d as GMat

class Test_main(unittest.TestCase):

    def test_Elastic(self):

        K = 12.3
        G = 45.6

        gamma = 0.02
        epsm = 0.12

        Eps = np.array(
            [[epsm, gamma, 0.0],
             [gamma, epsm, 0.0],
             [0.0, 0.0, epsm]])

        Sig = np.array(
            [[3.0 * K * epsm, 2.0 * G * gamma, 0.0],
             [2.0 * G * gamma, 3.0 * K * epsm, 0.0],
             [0.0, 0.0, 3.0 * K * epsm]])

        self.assertTrue(np.isclose(float(GMat.Epseq(Eps)), 2.0 / np.sqrt(3.0) * gamma))

        mat = GMat.Elastic(K, G)
        mat.setStrain(Eps)

        self.assertTrue(np.allclose(mat.Stress(), Sig))

    def test_Array2d(self):

        K = 12.3
        G = 45.6

        gamma = 0.02
        epsm = 0.12

        Eps = np.array(
            [[epsm, gamma, 0.0],
             [gamma, epsm, 0.0],
             [0.0, 0.0, epsm]])

        Sig = np.array(
            [[3.0 * K * epsm, 2.0 * G *gamma, 0.0],
             [2.0 * G *gamma, 3.0 * K * epsm, 0.0],
             [0.0, 0.0, 3.0 * K * epsm]])

        nelem = 3
        nip = 2
        mat = GMat.Array2d([nelem, nip])
        ndim = 3

        I = np.ones([nelem, nip], dtype='int')
        mat.setElastic(I, K, G)

        eps = np.zeros((nelem, nip, ndim, ndim))
        sig = np.zeros((nelem, nip, ndim, ndim))

        for e in range(nelem):
            for q in range(nip):
                fac = float((e + 1) * nip + (q + 1))
                eps[e, q, :, :] = fac * Eps
                sig[e, q, :, :] = fac * Sig

        mat.setStrain(eps)

        self.assertTrue(np.allclose(mat.Stress(), sig))

if __name__ == '__main__':

    unittest.main()
