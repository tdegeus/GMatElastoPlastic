import GMatElastoPlastic.Cartesian3d as GMat
import matplotlib.pyplot as plt
import numpy as np

try:
    plt.style.use(["goose", "goose-latex"])
except FileNotFoundError:
    pass

# define material model

mat = GMat.LinearHardening0d(K=10.0, G=1.0, sigy0=0.1, H=0.2)

# pre-loading

ninc = 301

epseq = np.zeros(ninc)
sigeq = np.zeros(ninc)

for igamma, gamma in enumerate(np.linspace(0.0, 0.1, ninc)):

    mat.increment()

    mat.Eps = np.array(
        [
            [0.0, gamma, 0.0],
            [gamma, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )

    epseq[igamma] = GMat.Epseq(mat.Eps)
    sigeq[igamma] = GMat.Sigeq(mat.Sig)

# plot result

fig, ax = plt.subplots()

ax.plot(epseq, mat.sigy0 * np.ones(epseq.shape), ls="--", lw=1, c="k")
ax.plot(epseq, sigeq, c="k")

ax.set_xlabel(r"$\varepsilon_\mathrm{eq}$")
ax.set_ylabel(r"$\sigma_\mathrm{eq}$")

fig.savefig("stress-strain.pdf")
plt.show()
