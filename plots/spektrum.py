from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt


folder = Path(
    r"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data"
)

filename = "eigenvalues_1D_W0.5_PBCfalse_Lx100_seed1234.data"
filepath = folder / filename

eigenvalues = np.loadtxt(filepath)
n = np.arange(len(eigenvalues))


# ODCZYT PARAMETROW Z NAZWY PLIKU
dimension_match = re.search(r"eigenvalues_(\dD)", filename)
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
plt.scatter(n, eigenvalues, s=20)

# dodatkowe poziomy energii i odpowiadające im numery wartości własnych
energy_levels = [
    (-2.0507335971431253, 1),
    ( 2.0799342423272753, 100),
    ( -0.019594035122677766, 50)
]

for E, number in energy_levels:

    plt.axhline(
        y=E,
        linestyle="--",
        linewidth=1
    )

    plt.text(
        len(eigenvalues) * 0.98,
        E,
        f"E = {E:.4f} ({number})",
        ha="right",
        va="bottom",
        fontsize=14
    )

plt.xlabel("Numer wartości własnej n", fontsize=16)
plt.ylabel("Energia $E_n$", fontsize=16)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.title(
    f"Spektrum {dimension}\n"
    f"{size_text}, W = {W}, PBC = {PBC}, seed = {seed}",
    fontsize=18
)

plt.grid(alpha=0.3)
plt.tight_layout()

# ZAPIS
script_folder = Path(__file__).resolve().parent
plot_name = filename.replace("eigenvalues_", "spektrum_")
plot_name = Path(plot_name).stem + ".png"
plot_path = script_folder / plot_name

plt.savefig(plot_path, dpi=300)

print("Wykres zapisany do:")
print(plot_path)

plt.show()