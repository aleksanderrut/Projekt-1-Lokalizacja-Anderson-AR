from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Folder z wynikami
FOLDER = Path(
    r"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data"
)

# PARAMETRY
W_VALUES = [1.0, 1.5, 2.0, 5.0]
R = 100
PBC = False
L = 100000

# Skala logarytmiczna
LOG_X = False
LOG_Y = False

# Zakres osi (None = automatyczny)
X_MIN = None
X_MAX = None
Y_MIN = None
Y_MAX = None

# Folder na wykresy
FOLDER_WYKRESY = FOLDER.parent / "plots"
FOLDER_WYKRESY.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 6))

# Każda krzywa odpowiada innej wartości W
for W in W_VALUES:

    filename = (
        f"gamma_vs_E_R{R}_W{W}_"
        f"PBC{str(PBC).lower()}_Lx{L}.data"
    )

    filepath = FOLDER / filename

    if not filepath.exists():
        print(f"Nie znaleziono pliku: {filename}")
        continue

    # Pierwsza kolumna: E
    # Druga kolumna: średni wykładnik Lapunowa
    data = np.loadtxt(filepath, ndmin=2)

    E = data[:, 0]
    gamma_mean = data[:, 1]

    # Usunięcie niedozwolonych wartości
    mask = np.isfinite(E) & np.isfinite(gamma_mean)

    if LOG_X:
        mask &= E > 0

    # Długość lokalizacji wymaga gamma > 0
    mask &= gamma_mean > 0

    E = E[mask]
    gamma_mean = gamma_mean[mask]

    # Obliczenie długości lokalizacji
    xi = 1.0 / gamma_mean

    # Sortowanie energii
    indeksy_sort = np.argsort(E)
    E = E[indeksy_sort]
    xi = xi[indeksy_sort]

    if len(E) == 0:
        print(f"Brak poprawnych danych dla W={W}")
        continue

    print(f"W={W}: narysowano {len(E)} punktów")

    ax.plot(
        E,
        xi,
        "-",
        linewidth=1.5,
        label=rf"$W={W:g}$"
    )

# Opisy osi
ax.set_xlabel(r"$E$", fontsize=14)
ax.set_ylabel(r"$\xi_L(E)=1/\langle\gamma_L(E)\rangle_R$", fontsize=14)

# Tytuł
ax.set_title(
    rf"$t=1,\quad L={L},\quad R={R}$",
    fontsize=14
)

# Skale osi
if LOG_X:
    ax.set_xscale("log")

if LOG_Y:
    ax.set_yscale("log")

# Zakres osi
ax.set_xlim(X_MIN, X_MAX)
ax.set_ylim(Y_MIN, Y_MAX)

ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=11)

fig.tight_layout()

# Zapis wykresu
filename_plot = (
    f"xi_vs_E_L{L}_R{R}"
    f"_logX{LOG_X}_logY{LOG_Y}.png"
)

fig.savefig(FOLDER_WYKRESY / filename_plot, dpi=300)

plt.show()
