import h5py
import numpy as np
import GMatElastoPlastic.Cartesian3d as GMat

with h5py.File('Cartesian3d_random.hdf5', 'w') as data:

    nelem = 1000
    nip = 4
    iden = 2.0 * np.random.random([nelem, nip])
    iden = np.where(iden < 1.0, 0.0, iden)
    iden = np.where(iden >= 1.0, 1.0, iden)
    iden = iden.astype(np.int)

    shape = np.array([nelem, nip], np.int)

    data['/shape'] = shape

    mat = GMat.Array2d(shape)

    I = np.where(iden == 0, 1, 0).astype(np.int)
    n = np.sum(I)
    idx = np.zeros(I.size, np.int)
    idx[np.argwhere(I.ravel() == 1).ravel()] = np.arange(n)
    idx = idx.reshape(I.shape)
    K = np.ones(n)
    G = np.ones(n)
    sigy0 = 0.1 * np.ones(n)
    H = np.ones(n)

    data['/LinearHardening/I'] = I
    data['/LinearHardening/idx'] = idx
    data['/LinearHardening/K'] = K
    data['/LinearHardening/G'] = G
    data['/LinearHardening/sigy0'] = sigy0
    data['/LinearHardening/H'] = H

    mat.setLinearHardening(I, idx, K, G, sigy0, H)

    I = np.where(iden == 1, 1, 0).astype(np.int)
    n = np.sum(I)
    idx = np.zeros(I.size, np.int)
    idx[np.argwhere(I.ravel() == 1).ravel()] = np.arange(n)
    idx = idx.reshape(I.shape)
    K = np.ones(n)
    G = np.ones(n)

    data['/Elastic/I'] = I
    data['/Elastic/idx'] = idx
    data['/Elastic/K'] = K
    data['/Elastic/G'] = G

    mat.setElastic(I, idx, K, G)

    for i in range(20):

        mat.increment()

        GradU = 20 * np.random.random([nelem, nip, 3, 3])

        data['/random/{0:d}/GradU'.format(i)] = GradU

        Eps = np.einsum('...ijkl,...lk->...ij', mat.I4s(), GradU)
        mat.setStrain(Eps)

        data['/random/{0:d}/Stress'.format(i)] = mat.Stress()
        data['/random/{0:d}/Tangent'.format(i)] = mat.Tangent()
        data['/random/{0:d}/Epsp'.format(i)] = mat.Epsp()

