
import GMatElastoPlastic.Cartesian3d as GMat
import numpy as np

def EQ(a,b):
  assert np.abs(a-b) < 1.e-12

def ALLEQ(a, b):
  assert np.allclose(a, b)

K = 12.3
G = 45.6

gamma = 0.02
epsm  = 0.12

Eps = np.array(
    [[epsm, gamma, 0.0],
     [gamma, epsm, 0.0],
     [0.0, 0.0, epsm]])

# Elastic

mat = GMat.Elastic(K, G)

Sig = mat.Stress(Eps)

EQ(Sig[0,0], 3.0 * K * epsm)
EQ(Sig[1,1], 3.0 * K * epsm)
EQ(Sig[2,2], 3.0 * K * epsm)
EQ(Sig[0,1], 2.0 * G * gamma)
EQ(Sig[1,0], 2.0 * G * gamma)
EQ(Sig[0,2], 0)
EQ(Sig[1,2], 0)
EQ(Sig[2,0], 0)
EQ(Sig[2,1], 0)

# Matrix

nelem = 2
nip = 2
mat = GMat.Matrix(nelem, nip)

# all rows: elastic
I = np.ones([nelem, nip], dtype='int')
mat.setElastic(I, K, G)

eps = np.zeros((nelem, nip, 3, 3))
for i in range(3):
    for j in range(3):
        eps[:, :, i, j] = Eps[i, j]

sig = mat.Stress(eps)

for e in range(nelem):
    for q in range(nip):

        EQ(sig[e,q,0,0], 3.0 * K * epsm)
        EQ(sig[e,q,1,1], 3.0 * K * epsm)
        EQ(sig[e,q,2,2], 3.0 * K * epsm)
        EQ(sig[e,q,0,1], 2.0 * G * gamma)
        EQ(sig[e,q,0,1], 2.0 * G * gamma)

ALLEQ(sig[:,:,0,2], 0.0)
ALLEQ(sig[:,:,1,2], 0.0)
ALLEQ(sig[:,:,2,0], 0.0)
ALLEQ(sig[:,:,2,1], 0.0)

print('All checks passed')
