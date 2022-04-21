import GMatElastoPlastic.Cartesian3d as GMat
import GMatTensor.Cartesian3d as tensor
import h5py
import numpy as np

with h5py.File("Cartesian3d_random.hdf5", "w") as data:

    shape = [1000, 4]

    mat = GMat.LinearHardening2d(
        K=np.random.random(shape),
        G=np.random.random(shape),
        sigy0=np.random.random(shape),
        H=np.random.random(shape),
    )

    data["K"] = mat.K
    data["G"] = mat.G
    data["sigy0"] = mat.sigy0
    data["H"] = mat.H

    for i in range(20):

        mat.increment()
        mat.Eps = tensor.Sym(20 * np.random.random(shape + [3, 3]))

        data[f"/data/{i:d}/Eps"] = mat.Eps
        data[f"/data/{i:d}/Stress"] = mat.Sig
        data[f"/data/{i:d}/Tangent"] = mat.C
        data[f"/data/{i:d}/epsp"] = mat.epsp
