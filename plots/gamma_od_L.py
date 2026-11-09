
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Folder z wynikami
FOLDER = Path(
    r"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data"
)

# PARAMETRY
E_VALUES = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
W = 2.0
R = 100
PBC = False
L_MAX = 1000000

# Zakres osi X
X_MIN = 0
X_MAX = 100000

# Linia pionowa
X_LINE = 10000
SHOW_LINE = True

# Folder na wykresy
FOLDER_WYKRESY = FOLDER.parent / "plots"
FOLDER_WYKRESY.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 6))

for E in E_VALUES:
    filename = f"gamma_1D_E{E}_R{R}_W{W}_PBC{str(PBC).lower()}_Lx{L_MAX}.data"
    filepath = FOLDER / filename

    if not filepath.exists():
        print(f"Nie znaleziono pliku: {filename}")
        continue

    data = np.loadtxt(filepath)

    L = data[:, 0]
    gamma = data[:, 1]

    ax.plot(
        L, gamma, "-o",
        markersize=3,
        linewidth=1.2,
        label=rf"$E={E:g}$"
    )

# Linia pionowa
if SHOW_LINE:
    ax.axvline(
        x=X_LINE,
        color="black",
        linestyle="--",
        linewidth=1.5,
        label=rf"$L={X_LINE}$"
    )

ax.set_xlabel(r"$L$", fontsize=14)
ax.set_ylabel(r"$\langle \gamma_L \rangle_R$", fontsize=14)
ax.set_title(rf"$W={W},\ R={R},\ t=1$", fontsize=14)

# Zakres osi X
ax.set_xlim(X_MIN, X_MAX)

ax.grid(True, alpha=0.3)
ax.legend(fontsize=11)

fig.tight_layout()

# Zapis wykresu
filename_plot = f"gamma_vs_L_W{W}_R{R}_Lmax{L_MAX}.png"
fig.savefig(FOLDER_WYKRESY / filename_plot, dpi=300)

plt.show()
