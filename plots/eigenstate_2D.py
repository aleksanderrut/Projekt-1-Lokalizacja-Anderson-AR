
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt


# ==========================================
# USTAWIENIA
# ==========================================

folder = Path(
    r"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data"
)

filename = "eigenvectors_2D_W0.5_PBCtrue_Lx50_Ly50_seed1234.data"
filepath = folder / filename

# Numery stanów do przedstawienia
states = [1, 1250, 2500]

# Skala kolorów
cmap = "inferno"


# ==========================================
# ODCZYT PARAMETRÓW Z NAZWY PLIKU
# ==========================================

Lx = int(re.search(r"_Lx(\d+)", filename).group(1))
Ly = int(re.search(r"_Ly(\d+)", filename).group(1))
W = re.search(r"_W([^_]+)", filename).group(1)
PBC = re.search(r"_PBC([^_]+)", filename).group(1)
seed = re.search(r"_seed(\d+)", filename).group(1)


# ==========================================
# WCZYTANIE WEKTORÓW WŁASNYCH
# ==========================================

eigenvectors = np.loadtxt(filepath)

N = Lx * Ly

if eigenvectors.shape[0] != N:
    raise ValueError(
        f"Niepoprawny rozmiar wektorów: "
        f"{eigenvectors.shape[0]} zamiast {N}"
    )

if max(states) > eigenvectors.shape[1] or min(states) < 1:
    raise ValueError("Numer stanu poza zakresem!")


# ==========================================
# RYSOWANIE HEATMAP
# ==========================================

fig, axes = plt.subplots(
    1,
    len(states),
    figsize=(6 * len(states), 6),
    constrained_layout=True
)

axes = np.atleast_1d(axes)

for ax, state in zip(axes, states):

    # Wektor własny: kolumna nr state - 1
    psi = eigenvectors[:, state - 1]

    # Kwadrat modułu funkcji falowej
    probability = np.abs(psi) ** 2

    # Zamiana wektora 2500 na macierz 50 x 50
    probability_2D = probability.reshape(Ly, Lx)

    # Heatmapa
    im = ax.imshow(
        probability_2D,
        origin="lower",
        cmap=cmap,
        interpolation="nearest",
        extent=[0.5, Lx + 0.5, 0.5, Ly + 0.5],
        aspect="equal"
    )

    ax.set_title(
        f"Stan n = {state}",
        fontsize=16
    )

    ax.set_xlabel("x", fontsize=14)
    ax.set_ylabel("y", fontsize=14)

    ax.tick_params(labelsize=12)

    cbar = fig.colorbar(
        im,
        ax=ax,
        fraction=0.046,
        pad=0.04
    )

    cbar.set_label(
        r"$|\psi(x,y)|^2$",
        fontsize=14
    )


fig.suptitle(
    "Kwadraty modułów wektorów własnych, 2D\n"
    f"{Lx} × {Ly}, W = {W}, PBC = {PBC}, seed = {seed}",
    fontsize=18
)


# ==========================================
# ZAPIS
# ==========================================

script_folder = Path(__file__).resolve().parent

base_name = Path(
    filename.replace("eigenvectors_", "")
).stem

states_text = "_".join(str(s) for s in states)

plot_name = f"modul2_stany_{states_text}_{base_name}.png"
plot_path = script_folder / plot_name

plt.savefig(plot_path, dpi=300, bbox_inches="tight")

print("Wykres zapisany do:")
print(plot_path)

plt.show()
