from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Folder z wynikami
FOLDER = Path(
    r"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data"
)

# PARAMETRY
L_VALUES = [377, 610, 987, 1597, 2584]
LAMBDA_RANGE = "lambda_0.2-4.0"
PBC = True
t = 1.0

# Wybór wykresu
MODE = "IPR"       # Średnia IPR
# MODE = "SIGMA"       # Względna wariancja IPR

# Maksimum każdej krzywej w legendzie
SHOW_MAX_IN_LEGEND = False

# Skala logarytmiczna
LOG_X = False
LOG_Y = False

# Zakres osi (None = automatyczny)
X_MIN = None
X_MAX = None
Y_MIN = None
Y_MAX = None

# Pionowa linia przy lambda/t = 2
SHOW_LINE = True
X_LINE = 2.0

# Folder na wykresy
FOLDER_WYKRESY = FOLDER.parent / "plots"
FOLDER_WYKRESY.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 6))

# Każda krzywa odpowiada innemu rozmiarowi L
for L in L_VALUES:

    filename = (
        f"IPR_sigma_rel_AA_W{LAMBDA_RANGE}_"
        f"PBC{str(PBC).lower()}_Lx{L}.data"
    )

    filepath = FOLDER / filename

    if not filepath.exists():
        print(f"Nie znaleziono pliku: {filename}")
        continue

    # Kolumna 1: lambda
    # Kolumna 2: średnia IPR
    # Kolumna 3: względna wariancja IPR
    data = np.loadtxt(filepath, ndmin=2)

    lam = data[:, 0]
    IPR_mean = data[:, 1]
    sigma_rel = data[:, 2]

    # Wybór wielkości do wyświetlenia
    if MODE == "IPR":
        y = IPR_mean
    elif MODE == "SIGMA":
        y = sigma_rel
    else:
        raise ValueError("Nieznany MODE")

    # Usunięcie niedozwolonych wartości
    mask = np.isfinite(lam) & np.isfinite(y)

    if LOG_X:
        mask &= lam > 0

    if LOG_Y:
        mask &= y > 0

    lam = lam[mask]
    y = y[mask]

    if len(lam) == 0:
        print(f"Brak poprawnych danych dla L={L}")
        continue

    # Sortowanie wartości lambda
    indeksy_sort = np.argsort(lam)
    lam = lam[indeksy_sort]
    y = y[indeksy_sort]

    print(f"L={L}: narysowano {len(lam)} punktów")

    # Położenie maksimum na osi X
    indeks_max = np.argmax(y)
    lam_max = lam[indeks_max]

    # Podpis w legendzie
    if SHOW_MAX_IN_LEGEND:
        label = rf"$L={L}$, maks dla $\lambda/t={lam_max:.3f}$"
    else:
        label = rf"$L={L}$"

    ax.plot(
        lam,
        y,
        "-o",
        markersize=4,
        linewidth=1.5,
        label=label
    )

# Pionowa linia przy lambda/t = 2
if SHOW_LINE:
    ax.axvline(
        x=X_LINE,
        color="gray",
        linestyle="--",
        linewidth=1.5,
        label=rf"$\lambda/t={X_LINE:g}$"
    )

# Opisy osi
ax.set_xlabel(r"$\lambda/t$", fontsize=14)

if MODE == "IPR":
    ax.set_ylabel(
        r"$\langle \mathrm{IPR} \rangle$ - Średnia IPR",
        fontsize=14
    )

elif MODE == "SIGMA":
    ax.set_ylabel(
        r"$\sigma_{\mathrm{rel}}^2$ - Względna wariancja IPR",
        fontsize=14
    )

# Tytuł
ax.set_title(
    rf"$t={t:g},\quad \mathrm{{PBC}}={str(PBC).lower()},\quad$"
    r"Środek widma energetycznego: "
    r"indeksy stanów od $\approx 0.35L$ do $\approx 0.65L$",
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
ax.legend(fontsize=10)

fig.tight_layout()

# Zapis wykresu
filename_plot = (
    f"AA_{MODE}_vs_lambda_PBC{str(PBC).lower()}"
    f"_logX{LOG_X}_logY{LOG_Y}.png"
)

fig.savefig(FOLDER_WYKRESY / filename_plot, dpi=300)

plt.show()
