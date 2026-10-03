from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt


folder = Path(
    r"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data"
)

filename = "eigenvectors_1D_W0.5_PBCfalse_Lx100_seed1234.data"
filepath = folder / filename

states = [1, 50, 100]

energies = {
    1: 2.0799,
    50: -0.0196,
    100: -2.0507
}

eigenvectors = np.loadtxt(filepath)

# ODCZYT PARAMETROW Z NAZWY PLIKU
dimension_match = re.search(r"eigenvectors_(\dD)", filename)
W_match = re.search(r"_W([^_]+)", filename)
Lx_match = re.search(r"_Lx(\d+)", filename)
Ly_match = re.search(r"_Ly(\d+)", filename)
Lz_match = re.search(r"_Lz(\d+)", filename)
PBC_match = re.search(r"_PBC([^_\.]+)", filename)
seed_match = re.search(r"_seed(\d+)", filename)

dimension = dimension_match.group(1) if dimension_match else "?"
W = W_match.group(1) if W_match else "?"
Lx = Lx_match.group(1) if Lx_match else None
Ly = Ly_match.group(1) if Ly_match else None
Lz = Lz_match.group(1) if Lz_match else None
PBC = PBC_match.group(1) if PBC_match else "?"
seed = seed_match.group(1) if seed_match else "?"

if dimension == "1D":
    size_text = f"L = {Lx}"

elif dimension == "2D":
    size_text = f"{Lx} × {Ly}"

elif dimension == "3D":
    size_text = f"{Lx} × {Ly} × {Lz}"

else:
    size_text = "?"


# WYKRES
plt.figure(figsize=(9, 6))

for state in states:

    psi = eigenvectors[:, state - 1]
    probability = np.abs(psi) ** 2
    n = np.arange(1, len(probability) + 1)

    plt.plot(
        n,
        probability,
        marker="o",
        markersize=4,
        linewidth=1.5,
        label=f"E = {energies[state]:.4f}"
    )

plt.xlabel("Numer węzła", fontsize=16)
plt.ylabel(r"$|\psi(i)|^2$", fontsize=16)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.title(
    f"Kwadraty modułów wektorów własnych, {dimension}\n"
    f"{size_text}, W = {W}, PBC = {PBC}, seed = {seed}",
    fontsize=18
)

plt.legend(fontsize=13)
plt.grid(alpha=0.3)
plt.tight_layout()

# ZAPIS
script_folder = Path(__file__).resolve().parent
base_name = filename.replace("eigenvectors_", "")
base_name = Path(base_name).stem
states_text = "_".join(str(state) for state in states)
plot_name = f"modul2_stany_{states_text}_{base_name}.png"
plot_path = script_folder / plot_name

plt.savefig(plot_path, dpi=300)

print("Wykres zapisany do:")
print(plot_path)

plt.show()