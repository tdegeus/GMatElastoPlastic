import GMatElastoPlastic as gmat
import matplotlib.pyplot as plt
import numpy as np

try:
    plt.style.use(['goose', 'goose-latex'])
except:
    pass

# tensor operations

def ddot42(A4, B2):
    return np.einsum('ijkl,lk->ij', A4, B2)


def ddot22(A2, B2):
    return np.einsum('ij,ji', A2, B2)

I4d = gmat.Cartesian3d.I4d()

# define material model

mat = gmat.Cartesian3d.LinearHardening(10.0, 1.0, 0.1, 0.2)

# pre-loading

ninc = 301

epseq = np.zeros(ninc)
sigeq = np.zeros(ninc)

for igamma, gamma in enumerate(np.linspace(0.0, 0.1, ninc)):

    mat.increment()

    Eps = np.array([
        [0.0, gamma,   0.0],
        [gamma,   0.0,   0.0],
        [0.0,   0.0,   0.0],
    ])

    Sig = mat.Stress(Eps)

    Epsd = ddot42(I4d, Eps)
    Sigd = ddot42(I4d, Sig)

    epseq[igamma] = np.sqrt(2./3.*ddot22(Epsd, Epsd))
    sigeq[igamma] = np.sqrt(3./2.*ddot22(Sigd, Sigd))

# plot result

fig, ax = plt.subplots()

ax.plot(epseq, mat.sigy0() * np.ones(epseq.shape), ls='--', lw=1, c='k')

ax.plot(epseq, sigeq, c='k')

ax.set_xlabel(r'$\varepsilon_\mathrm{eq}$')
ax.set_ylabel(r'$\sigma_\mathrm{eq}$')

plt.savefig('stress-strain.pdf')
plt.show()
