from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============================================================
# FOLDER Z DANYMI
# ============================================================
FOLDER = Path(
    r"C:\Users\aleks\Desktop\PWR\Kwantowe układy wielu cząstek z elementami kwantowej fizyki statystycznej\Projekt 1 Lokalizacja Andersona\data"
)

# ============================================================
# PARAMETRY
# ============================================================
L_VALUES = [377, 610, 987, 1597, 2584]
PBC = True
t = 1.0

# Środkowa część widma
SRODEK_MIN = 0.35
SRODEK_MAX = 0.65

# Wybór modelu
# MODEL = "BOTH"        # Anderson i Aubry-André
# MODEL = "AA"        # Tylko Aubry-André
MODEL = "ANDERSON"  # Tylko Anderson

# Wybór współczynnika
MODE = "IPR"          # Średni IPR
#MODE = "SIGMA"      # Względna wariancja IPR

# Skala logarytmiczna
LOG_X = False
LOG_Y = False

# Zakres osi
X_MIN = 0.0
X_MAX = 8.0

Y_MIN = None
Y_MAX = None

# Zapis wykresu
SAVE_PLOT = True

FOLDER_WYKRESY = FOLDER.parent / "plots"
FOLDER_WYKRESY.mkdir(parents=True, exist_ok=True)

# ============================================================
# NAZWY PLIKÓW
# ============================================================
def file_anderson(L):
    return FOLDER / (
        f"IPR_sigma_anderson_Ww_0.2-4.0_"
        f"PBC{str(PBC).lower()}_Lx{L}.data"
    )

def file_aa(L):
    return FOLDER / (
        f"IPR_sigma_rel_AA_Wlambda_0.2-4.0_"
        f"PBC{str(PBC).lower()}_Lx{L}.data"
    )

# ============================================================
# USTAWIENIA WYKRESU
# ============================================================
fig, ax = plt.subplots(figsize=(10, 6))

# Jeden kolor dla danego L w obu modelach
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# Kolumna danych
if MODE == "IPR":
    kolumna = 1
    ylabel = r"$\langle \mathrm{IPR} \rangle$"
    title_main = "Średni IPR"

elif MODE == "SIGMA":
    kolumna = 2
    ylabel = r"$\sigma_{\mathrm{rel}}$"
    title_main = r"Względna wariancja IPR, $\sigma_{\mathrm{rel}}$"

else:
    raise ValueError("Nieznany MODE!")

# ============================================================
# STYLE LINII
# ============================================================
# Aubry-André: linia ciągła
STYLE_AA = "-"
MARKER_AA = None

# Anderson: linia kreska-kropka z markerami
STYLE_ANDERSON = "-."
MARKER_ANDERSON = "o"

# ============================================================
# RYSOWANIE
# ============================================================
for color, L in zip(colors, L_VALUES):

    # --------------------------------------------------------
    # MODEL ANDERSONA
    # --------------------------------------------------------
    if MODEL in ("BOTH", "ANDERSON"):

        plik_and = file_anderson(L)

        if plik_and.exists():

            dane_and = np.loadtxt(plik_and)

            x_and = dane_and[:, 0] / t
            y_and = dane_and[:, kolumna]

            ax.plot(
                x_and,
                y_and,
                linestyle=STYLE_ANDERSON,
                linewidth=1.8,
                marker=MARKER_ANDERSON,
                markersize=4.5,
                color=color,
                label=f"L={L}"
            )

        else:
            print(f"Brak pliku Andersona: {plik_and}")

    # --------------------------------------------------------
    # MODEL AUBRY-ANDRÉ
    # --------------------------------------------------------
    if MODEL in ("BOTH", "AA"):

        plik_aa = file_aa(L)

        if plik_aa.exists():

            dane_aa = np.loadtxt(plik_aa)

            x_aa = dane_aa[:, 0] / t
            y_aa = dane_aa[:, kolumna]

            ax.plot(
                x_aa,
                y_aa,
                linestyle=STYLE_AA,
                linewidth=2.2,
                marker=MARKER_AA,
                color=color,
                label=f"L={L}"
            )

        else:
            print(f"Brak pliku Aubry-André: {plik_aa}")

# ============================================================
# TYTUŁ I OPISY OSI
# ============================================================
if MODEL == "BOTH":
    title_model = "Anderson i Aubry-André"
    xlabel = r"$W/t$ (Anderson) / $\lambda/t$ (Aubry-André)"

elif MODEL == "AA":
    title_model = "Aubry-André"
    xlabel = r"$\lambda/t$"

elif MODEL == "ANDERSON":
    title_model = "Anderson"
    xlabel = r"$W/t$"

else:
    raise ValueError("Nieznany MODEL!")

# Tytuł jak na poprzednim wykresie
ax.set_title(
    rf"$t={t:g}$,   $\mathrm{{PBC}}={str(PBC).lower()}$,   "
    rf"Środek widma energetycznego: indeksy stanów od "
    rf"$\approx {SRODEK_MIN}L$ do $\approx {SRODEK_MAX}L$",
    fontsize=16,
    pad=10
)

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)

# ============================================================
# SKALE I ZAKRESY OSI
# ============================================================
if LOG_X:
    ax.set_xscale("log")

    dodatnie_x = [
        x for line in ax.lines
        for x in line.get_xdata()
        if x > 0
    ]

    if dodatnie_x:
        ax.set_xlim(
            max(X_MIN, min(dodatnie_x)),
            X_MAX
        )
else:
    ax.set_xlim(X_MIN, X_MAX)

if LOG_Y:
    ax.set_yscale("log")

if Y_MIN is not None or Y_MAX is not None:
    ax.set_ylim(Y_MIN, Y_MAX)

ax.grid(True, alpha=0.3)

# ============================================================
# LEGENDY
# ============================================================

# Legenda kolorów — rozmiary L
legend_L = [
    Line2D(
        [0], [0],
        color=color,
        linewidth=2,
        linestyle="-",
        label=f"L={L}"
    )
    for color, L in zip(colors, L_VALUES)
]

leg1 = ax.legend(
    handles=legend_L,
    title="Rozmiar układu",
    loc="upper left",
    fontsize=10
)

ax.add_artist(leg1)

# Legenda stylów linii — modele
if MODEL == "BOTH":

    legend_models = [
        Line2D(
            [0], [0],
            color="black",
            linewidth=1.8,
            linestyle="-.",
            marker="o",
            markersize=5,
            label="Anderson"
        ),
        Line2D(
            [0], [0],
            color="black",
            linewidth=2.2,
            linestyle="-",
            label="Aubry-André"
        )
    ]

    ax.legend(
        handles=legend_models,
        title="Model",
        loc="upper right",
        fontsize=10
    )

plt.tight_layout()

# ============================================================
# ZAPIS WYKRESU
# ============================================================
if SAVE_PLOT:

    outname = f"{MODE}_{MODEL}_porownanie.png"

    outfile = FOLDER_WYKRESY / outname

    plt.savefig(
        outfile,
        dpi=300,
        bbox_inches="tight"
    )

    print(f"Zapisano wykres: {outfile}")

plt.show()
