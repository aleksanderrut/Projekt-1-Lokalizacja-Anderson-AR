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
L_MAX = 10000000

# Skala logarytmiczna
LOG_X = True
LOG_Y = True 

# Zakres osi X
X_MIN = 0
X_MAX = 10000000

# Linia pionowa
X_LINE = 100000
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

    # Usuwamy wartości niedozwolone dla skali logarytmicznej
    mask = np.isfinite(L) & np.isfinite(gamma)

    if LOG_X:
        mask &= L > 0

    if LOG_Y:
        mask &= gamma > 0

    L = L[mask]
    gamma = gamma[mask]

    ax.plot(
        L, gamma, "-o",
        markersize=3,
        linewidth=1.2,
        label=rf"$E={E:g}$"
    )

# Linia pionowa
if SHOW_LINE and (not LOG_X or X_LINE > 0):
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

# Skala osi
if LOG_X:
    ax.set_xscale("log")

if LOG_Y:
    ax.set_yscale("log")

# Zakres osi X
if LOG_X and X_MIN <= 0:
    ax.set_xlim(right=X_MAX)
else:
    ax.set_xlim(X_MIN, X_MAX)

ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=11)

fig.tight_layout()

# Zapis wykresu
filename_plot = (
    f"gamma_vs_L_W{W}_R{R}_Lmax{L_MAX}"
    f"_logX{LOG_X}_logY{LOG_Y}.png"
)

fig.savefig(FOLDER_WYKRESY / filename_plot, dpi=300)

plt.show()
