import GMatElastoPlastic as gmat
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(['goose', 'goose-latex'])

# tensor operations

def ddot42(A4, B2):
    return np.einsum('ijkl,lk->ij', A4, B2)

def ddot22(A2, B2):
    return np.einsum('ij,ji', A2, B2)

def norm(A2):
    return np.abs(np.einsum('ij,ji', A2, A2))

# define material model

mat = gmat.Cartesian3d.LinearHardening(10.0, 1.0, 0.1, 0.2)

# pre-loading

ninc = 301

for igamma, gamma in enumerate(np.linspace(0.0, 0.1, ninc)):

    mat.increment()

    Eps0 = np.array([
        [0.0, gamma, 0.0],
        [gamma, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ])

    Sig0, C4 = mat.Tangent(Eps0)

# consistency check

x = np.logspace(-16, 0, 100)
y = np.zeros(x.shape)

for i in range(len(x)):

    dEps = np.random.random((3, 3)) * x[i]
    dEps = 0.5 * (dEps + dEps.T)

    Sig = mat.Stress(Eps0 + dEps)

    dSig = Sig - Sig0

    y[i] = norm(dSig - ddot42(C4, dEps)) / norm(dSig)

# plot result

fig, ax = plt.subplots()

ax.plot(x, y, color='r', label=r'measurement')

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-18, 1e0])
ax.set_ylim([1e-18, 1e0])

ax.set_xlabel(r'$|| \delta \bm{\varepsilon} ||$')
ax.set_ylabel(r'$\eta$')

gplt.plot_powerlaw(-2, 0.0, 1.0, 0.5, axis=ax, units='relative', color='k', linewidth=1,
                   label=r'rounding error: $|| \delta \bm{\varepsilon} ||^{-2}$')

gplt.plot_powerlaw(+2, 0.5, 0.0, 0.5, axis=ax, units='relative', color='k', linewidth=1, linestyle='--',
                   label=r'linearisation error: $|| \delta \bm{\varepsilon} ||^{+2}$')

ax.legend()

plt.savefig('consistency.pdf')
plt.show()
