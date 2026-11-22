
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# ============================================================
# PARAMETRY
# ============================================================

plik = Path(
    r"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data\IPR_matrix_3D_R100_W0_PBCtrue_Lx10_Ly10_Lz10.data"
)

# Zakresy osi (None = automatyczny)
X_MIN = None
X_MAX = None

Y_MIN = None
Y_MAX = None

# Zakres skali kolorów
VMIN = -3.0
VMAX = -0.5

# Skala kolorystyczna
CMAP = "inferno"

# ============================================================
# ODCZYT PARAMETRÓW Z NAZWY PLIKU
# ============================================================

nazwa = plik.stem

dimension = re.search(r"IPR_matrix_(\dD)", nazwa).group(1)
R = int(re.search(r"_R(\d+)", nazwa).group(1))
PBC = re.search(r"_PBC(true|false)", nazwa).group(1)

Lx = int(re.search(r"_Lx(\d+)", nazwa).group(1))
Ly = int(re.search(r"_Ly(\d+)", nazwa).group(1))
Lz = int(re.search(r"_Lz(\d+)", nazwa).group(1))

print("Wymiar:", dimension)
print("R:", R)
print("PBC:", PBC)
print("Lx:", Lx)
print("Ly:", Ly)
print("Lz:", Lz)

# ============================================================
# WCZYTANIE MACIERZY IPR
# ============================================================

IPR_matrix = np.loadtxt(plik)

N_W, N_binow = IPR_matrix.shape

print("Rozmiar macierzy:", IPR_matrix.shape)

# ============================================================
# OŚ W I OŚ ENERGII
# ============================================================

# Parametry zgodne z kodem Julii
W_min = 1.0
W_max = 20.0
delta_W = 0.2

E_min = -16.0
E_max = 16.0

# Wartości W odpowiadające kolejnym wierszom
W_values = np.arange(
    W_min,
    W_max + delta_W / 2,
    delta_W
)

# Środki binów energetycznych
delta_E = (E_max - E_min) / N_binow

energie = E_min + (
    np.arange(N_binow) + 0.5
) * delta_E

if len(W_values) != N_W:
    raise ValueError(
        f"Liczba wartości W ({len(W_values)}) "
        f"nie zgadza się z liczbą wierszy macierzy ({N_W})."
    )

# ============================================================
# PRZYGOTOWANIE HEATMAPY
# ============================================================

# Pominięcie NaN, inf i wartości niedodatnich
maska = ~np.isfinite(IPR_matrix) | (IPR_matrix <= 0)

print("Liczba NaN:", np.isnan(IPR_matrix).sum())
print("Liczba inf:", np.isinf(IPR_matrix).sum())
print("Liczba pomijanych punktów:", maska.sum())

IPR_clean = IPR_matrix.copy()
IPR_clean[maska] = np.nan

# Logarytm dziesiętny
with np.errstate(divide="ignore", invalid="ignore"):
    IPR_log = np.log10(IPR_clean)

IPR_plot = np.ma.masked_invalid(IPR_log)

if IPR_plot.count() == 0:
    raise ValueError(
        "Brak poprawnych wartości IPR do narysowania."
    )

print("Minimum log10(IPR):", float(IPR_plot.min()))
print("Maksimum log10(IPR):", float(IPR_plot.max()))

# Mapa kolorów inferno
cmap = plt.get_cmap(CMAP).copy()
cmap.set_bad(color="white")

# ============================================================
# WYKRES
# ============================================================

fig, ax = plt.subplots(figsize=(10, 7))

im = ax.imshow(
    IPR_plot,
    origin="lower",
    aspect="auto",
    interpolation="nearest",
    extent=[
        E_min,
        E_max,
        W_min - delta_W / 2,
        W_max + delta_W / 2
    ],
    cmap=cmap,
    vmin=VMIN,
    vmax=VMAX
)

ax.set_xlim(X_MIN, X_MAX)
ax.set_ylim(Y_MIN, Y_MAX)

ax.set_xlabel("Energia E", fontsize=16)
ax.set_ylabel("Nieporządek W", fontsize=16)

ax.set_title(
    f"Participation Ratio {dimension}\n"
    f"{Lx} x {Ly} x {Lz}, "
    f"PBC = {PBC}, R = {R}",
    fontsize=17
)

ax.tick_params(axis="both", labelsize=13)

# ============================================================
# COLORBAR
# ============================================================

cbar = fig.colorbar(im, ax=ax)

cbar.set_label(
    r"$\log_{10}(\langle IPR\rangle)$",
    fontsize=15
)

# Podpisy od -3 do -0.5, co 0.5
ticks = np.arange(-3.0, -0.49, 0.5)

cbar.set_ticks(ticks)

# Bez niepotrzebnych miejsc po przecinku
cbar.ax.yaxis.set_major_formatter(
    FuncFormatter(lambda x, pos: f"{x:g}")
)

cbar.ax.tick_params(labelsize=12)

plt.tight_layout()

# ============================================================
# ZAPIS
# ============================================================

folder_wykresy = plik.parent.parent / "plots"
folder_wykresy.mkdir(parents=True, exist_ok=True)

plot_name = f"Heatmap_{nazwa}.png"
plot_path = folder_wykresy / plot_name

plt.savefig(
    plot_path,
    dpi=300,
    bbox_inches="tight"
)

print("Wykres zapisany do:")
print(plot_path)

plt.show()
