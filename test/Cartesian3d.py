
import GMatElastoPlastic.Cartesian3d as GMat
import numpy as np

# ==================================================================================================

def EQ(a,b):
  assert np.abs(a-b) < 1.e-12

def ALLEQ(a, b):
  assert np.allclose(a, b)

# ==================================================================================================

# material model
# - parameters
K = 12.3
G = 45.6
# - model
mat = GMat.Elastic(K,G)

# simple shear + volumetric deformation
# - parameters
gamma = 0.02
epsm  = 0.12
# - strain
Eps = [[epsm , gamma, 0.0 ],
       [gamma, epsm , 0.0 ],
       [0.0  , 0.0  , epsm]]
# - stress
Sig = mat.Stress(Eps)
# - analytical solution
EQ(Sig[0,0], 3.0 * K * epsm)
EQ(Sig[1,1], 3.0 * K * epsm)
EQ(Sig[2,2], 3.0 * K * epsm)
EQ(Sig[0,1], 2.0 * G * gamma)
EQ(Sig[1,0], 2.0 * G * gamma)
EQ(Sig[0,2], 0)
EQ(Sig[1,2], 0)
EQ(Sig[2,0], 0)
EQ(Sig[2,1], 0)

# ==================================================================================================

# parameters
K = 12.3
G = 45.6

# allocate matrix
nelem = 2
nip = 2
mat = GMat.Matrix(nelem, nip)

# all rows: elastic
I = np.ones([nelem, nip], dtype='int')
mat.setElastic(I,K,G)

# simple shear + volumetric deformation
# - parameters
gamma = 0.02;
epsm  = 0.12;
# - strain
Eps = np.zeros((nelem, nip, 3, 3))
Eps[:,:,0,0] = epsm
Eps[:,:,1,1] = epsm
Eps[:,:,2,2] = epsm
Eps[:,:,0,1] = gamma
Eps[:,:,1,0] = gamma
# - stress
Sig = mat.Stress(Eps)

# - analytical solution
EQ(Sig[0,0,0,0], 3.0 * K * epsm);  EQ(Sig[0,1,0,0], 3.0 * K * epsm)
EQ(Sig[0,0,1,1], 3.0 * K * epsm);  EQ(Sig[0,1,1,1], 3.0 * K * epsm)
EQ(Sig[0,0,2,2], 3.0 * K * epsm);  EQ(Sig[0,1,2,2], 3.0 * K * epsm)
EQ(Sig[0,0,0,1], 2.0 * G * gamma); EQ(Sig[0,1,0,1], 2.0 * G * gamma)
EQ(Sig[0,0,0,1], 2.0 * G * gamma); EQ(Sig[0,1,1,0], 2.0 * G * gamma)
EQ(Sig[1,0,0,0], 3.0 * K * epsm);  EQ(Sig[1,1,0,0], 3.0 * K * epsm)
EQ(Sig[1,0,1,1], 3.0 * K * epsm);  EQ(Sig[1,1,1,1], 3.0 * K * epsm)
EQ(Sig[1,0,2,2], 3.0 * K * epsm);  EQ(Sig[1,1,2,2], 3.0 * K * epsm)
EQ(Sig[1,0,0,1], 2.0 * G * gamma); EQ(Sig[1,1,0,1], 2.0 * G * gamma)
EQ(Sig[1,0,0,1], 2.0 * G * gamma); EQ(Sig[1,1,1,0], 2.0 * G * gamma)
ALLEQ(Sig[:,:,0,2], 0)
ALLEQ(Sig[:,:,1,2], 0)
ALLEQ(Sig[:,:,2,0], 0)
ALLEQ(Sig[:,:,2,1], 0)

# ==================================================================================================

print('All checks passed')
