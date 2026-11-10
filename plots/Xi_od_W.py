
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Folder z wynikami
FOLDER = Path(
    r"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data"
)

# PARAMETRY
E_VALUES = [0.0, 1.0, 1.5]
W_RANGE = "0.1-10.0"
R = 100
PBC = False
L = 100000

# Skala logarytmiczna
LOG_X = True
LOG_Y = True

# Funkcja referencyjna f(W) = W^a + b
SHOW_REFERENCE = True
A_REF = -2.0
B_REF = 25.0

# Zakres funkcji referencyjnej
W_REF_MIN = 0.1
W_REF_MAX = 10.0

# Zakres osi (None = automatyczny)
X_MIN = None
X_MAX = None
Y_MIN = None
Y_MAX = None

# Folder na wykresy
FOLDER_WYKRESY = FOLDER.parent / "plots"
FOLDER_WYKRESY.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 6))

# Każda krzywa odpowiada innej energii E
for E in E_VALUES:

    filename = (
        f"gamma_vs_W_E{E}_R{R}_W{W_RANGE}_"
        f"PBC{str(PBC).lower()}_Lx{L}.data"
    )

    filepath = FOLDER / filename

    if not filepath.exists():
        print(f"Nie znaleziono pliku: {filename}")
        continue

    # Pierwsza kolumna: W
    # Druga kolumna: średni wykładnik Lapunowa
    data = np.loadtxt(filepath, ndmin=2)

    W = data[:, 0]
    gamma_mean = data[:, 1]

    # Usunięcie niedozwolonych wartości
    mask = np.isfinite(W) & np.isfinite(gamma_mean)

    if LOG_X:
        mask &= W > 0

    # Długość lokalizacji wymaga gamma > 0
    mask &= gamma_mean > 0

    W = W[mask]
    gamma_mean = gamma_mean[mask]

    if len(W) == 0:
        print(f"Brak poprawnych danych dla E={E}")
        continue

    # Obliczenie długości lokalizacji
    xi = 1.0 / gamma_mean

    # Sortowanie wartości W
    indeksy_sort = np.argsort(W)
    W = W[indeksy_sort]
    xi = xi[indeksy_sort]

    print(f"E={E}: narysowano {len(W)} punktów")

    ax.plot(
        W,
        xi,
        "-",
        linewidth=1.5,
        label=rf"$E={E:g}$"
    )

# ============================================================
# FUNKCJA REFERENCYJNA
# ============================================================

if SHOW_REFERENCE:

    if LOG_X:
        W_ref = np.logspace(
            np.log10(W_REF_MIN),
            np.log10(W_REF_MAX),
            500
        )
    else:
        W_ref = np.linspace(
            W_REF_MIN,
            W_REF_MAX,
            500
        )

    # f(W) = W^a + b
    y_ref = B_REF * W_ref**A_REF

    # Usunięcie niedozwolonych wartości
    mask = np.isfinite(y_ref)

    if LOG_Y:
        mask &= y_ref > 0

    # Automatyczny wzór w legendzie
    label=rf"$f(W)={B_REF:g}\cdot W^{{{A_REF:g}}}$"

    ax.plot(
        W_ref[mask],
        y_ref[mask],
        "--",
        color="gray",
        linewidth=2,
        label=rf"$f(W)={B_REF:g}\cdot W^{{{A_REF:g}}}$"
    )

# Opisy osi
ax.set_xlabel(r"$W$", fontsize=14)
ax.set_ylabel(
    r"$\xi_L(W)=1/\langle\gamma_L(W)\rangle_R$",
    fontsize=14
)

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
    f"xi_vs_W_L{L}_R{R}"
    f"_logX{LOG_X}_logY{LOG_Y}.png"
)

fig.savefig(FOLDER_WYKRESY / filename_plot, dpi=300)

plt.show()
