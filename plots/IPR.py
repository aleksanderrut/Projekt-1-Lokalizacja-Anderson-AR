
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# =========================
# PARAMETRY
# =========================

folder = Path(
    r"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data"
)

# "binned" - IPR uśrednione w binach energetycznych
# "mean"   - IPR uśrednione po realizacjach
tryb = "binned"

# (L_x, W, PBC, R)
PARAMETRY = [
    (50,   8.0, True, 500),
    (100,  8.0, True, 500),
    (200,  8.0, True, 500),
    (500,  8.0, True, 500),
    (1000, 8.0, True, 500),
]

# Zakres osi
X_MIN = -6
X_MAX = 6

Y_MIN = None
Y_MAX = None


# =========================
# WYBÓR TRYBU
# =========================

if tryb == "binned":
    prefix = "IPR_1D"
elif tryb == "mean":
    prefix = "IPR_mean_1D"
else:
    raise ValueError("Nieznany tryb!")


# =========================
# WSPÓLNE PARAMETRY DO TYTUŁU
# =========================

_, W_title, PBC_title, R_title = PARAMETRY[0]

# Sprawdzenie, czy porównujemy tylko różne rozmiary L_x
if any(
    (W, PBC, R) != (W_title, PBC_title, R_title)
    for L_x, W, PBC, R in PARAMETRY
):
    raise ValueError(
        "Aby wyświetlać W, PBC i R w tytule, "
        "muszą być takie same dla wszystkich krzywych."
    )


# =========================
# WYKRES
# =========================

plt.figure(figsize=(9, 6))

for L_x, W, PBC, R in PARAMETRY:

    # =========================
    # WYBÓR PLIKU
    # =========================

    nazwa_podstawowa = (
        f"{prefix}_R{R}_W{W}_PBC{str(PBC).lower()}"
        f"_Lx{L_x}"
    )

    # Pliki bez seeda i z seedem
    pliki = []

    plik_bez_seeda = folder / f"{nazwa_podstawowa}.data"

    if plik_bez_seeda.exists():
        pliki.append(plik_bez_seeda)

    pliki.extend(
        folder.glob(f"{nazwa_podstawowa}_seed*.data")
    )

    if not pliki:
        print(f"Nie znaleziono pliku: {nazwa_podstawowa}")
        continue

    # Wybór najnowszego pliku
    plik = max(pliki, key=lambda p: p.stat().st_mtime)

    print("Wczytuję:", plik.name)

    # =========================
    # WCZYTANIE DANYCH
    # =========================

    data = np.loadtxt(plik, ndmin=2)

    energia = data[:, 0]
    IPR_mean = data[:, 1]

    # Pomijanie pustych binów
    if tryb == "binned" and data.shape[1] >= 3:
        liczba_stanow = data[:, 2]
        mask = liczba_stanow > 0
    else:
        mask = np.ones(len(energia), dtype=bool)

    mask &= np.isfinite(energia)
    mask &= np.isfinite(IPR_mean)
    mask &= IPR_mean > 0

    energia = energia[mask]
    IPR_mean = IPR_mean[mask]

    if len(energia) == 0:
        print("Brak poprawnych danych:", plik.name)
        continue

    # Participation Ratio
    PR = 1.0 / IPR_mean

    # Sortowanie względem energii
    indeksy = np.argsort(energia)
    energia = energia[indeksy]
    PR = PR[indeksy]

    # =========================
    # RYSOWANIE KRZYWEJ
    # =========================

    plt.plot(energia, PR, marker="o", markersize=3, linewidth=1.5, label=f"L = {L_x}")

# =========================
# WYGLĄD WYKRESU

plt.xlim(X_MIN, X_MAX)
plt.ylim(Y_MIN, Y_MAX)

plt.xlabel("Energia E", fontsize=16)
plt.ylabel(r"$1/\langle IPR\rangle$", fontsize=16)

plt.title(f"Participation Ratio 1D\n" f"W = {W_title}, PBC = {PBC_title}, R = {R_title}", fontsize=18)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
plt.grid(alpha=0.3)
plt.tight_layout()


# =========================
# ZAPIS

script_folder = Path(__file__).resolve().parent
plot_name = "fig3.png"
plot_path = script_folder / plot_name
plt.savefig(plot_path, dpi=300)
print("Wykres zapisany do:")
print(plot_path)

plt.show()
