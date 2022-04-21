import GMatElastoPlastic.Cartesian3d as GMat
import GMatTensor.Cartesian3d as tensor
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose", "goose-latex"])


def norm(A2):
    return tensor.A2_ddot_B2(A2, A2)


# define material model

mat = GMat.LinearHardening0d(K=10.0, G=1.0, sigy0=0.1, H=0.2)

# pre-loading

ninc = 301

for igamma, gamma in enumerate(np.linspace(0.0, 0.1, ninc)):

    mat.increment()

    Eps0 = np.array(
        [
            [0.0, gamma, 0.0],
            [gamma, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )

    mat.Eps = Eps0

Sig0 = np.copy(mat.Sig)
C4 = np.copy(mat.C)

# consistency check

x = np.logspace(-16, 0, 100)
y = np.zeros(x.shape)

for i in range(len(x)):

    dEps = np.random.random((3, 3)) * x[i]
    dEps = 0.5 * (dEps + dEps.T)

    mat.Eps = Eps0 + dEps
    dSig = mat.Sig - Sig0

    y[i] = norm(dSig - tensor.A4_ddot_B2(C4, dEps)) / norm(dSig)

# plot result

fig, ax = plt.subplots()

ax.plot(x, y, color="r", label=r"measurement")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlim([1e-18, 1e0])
ax.set_ylim([1e-18, 1e0])

ax.set_xlabel(r"$|| \delta \bm{\varepsilon} ||$")
ax.set_ylabel(r"$\eta$")

gplt.plot_powerlaw(
    -2,
    0.0,
    1.0,
    0.5,
    axis=ax,
    units="relative",
    color="k",
    linewidth=1,
    label=r"rounding error: $|| \delta \bm{\varepsilon} ||^{-2}$",
)

gplt.plot_powerlaw(
    +2,
    0.5,
    0.0,
    0.5,
    axis=ax,
    units="relative",
    color="k",
    linewidth=1,
    linestyle="--",
    label=r"linearisation error: $|| \delta \bm{\varepsilon} ||^{+2}$",
)

ax.legend()

fig.savefig("consistency.pdf")
plt.close(fig)
