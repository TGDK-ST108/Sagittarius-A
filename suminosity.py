import math
import pandas as pd
import matplotlib.pyplot as plt

# Constants
PMZ = 0.0102
M_sun = 1.98847e30
G = 6.67430e-11
c = 299792458
AU_m = 1.495978707e11
pc_to_AU = 206265.0

# Base formulas
def L_edd_erg_s(M_solar):
    return 1.26e38 * M_solar

def R_sub_AU(L_erg_s):
    return 0.5 * pc_to_AU * math.sqrt(L_erg_s / 1e45)

def R_s_AU(M_solar):
    M_kg = M_solar * M_sun
    Rs_m = 2 * G * M_kg / (c**2)
    return Rs_m / AU_m

def ESM_norm(PMZ, L, M):
    Ledd = L_edd_erg_s(M)
    Rs = R_s_AU(M)
    Rsub = R_sub_AU(L)
    return PMZ * math.sqrt(L / Ledd) * (Rsub / Rs), Ledd, Rsub, Rs

# Object data
objects = [
    ("Sagittarius A* (Quiescent)", 4.297e6, 1e36),
    ("Sagittarius A* (Flare)", 4.297e6, 1.5e39),
    ("Sun (Reference)", 1.0, 3.828e33),
]

rows = []
for name, M, L in objects:
    esm, Ledd, Rsub, Rs = ESM_norm(PMZ, L, M)
    rows.append({
        "Object": name,
        "Mass (Msun)": M,
        "Luminosity (erg/s)": L,
        "Eddington Luminosity (erg/s)": Ledd,
        "L/L_Edd": L / Ledd,
        "R_sub (AU)": Rsub,
        "R_s (AU)": Rs,
        "ESM_norm (dimensionless)": esm
    })

df = pd.DataFrame(rows)

# Plot ESM_norm vs Luminosity (log scale)
plt.figure(figsize=(8,6))
plt.loglog(df["Luminosity (erg/s)"], df["ESM_norm (dimensionless)"], marker='o', linestyle='-', color='blue')
for i, row in df.iterrows():
    plt.text(row["Luminosity (erg/s)"], row["ESM_norm (dimensionless)"]*1.1, row["Object"], fontsize=8, ha='center')
plt.xlabel("Luminosity (erg/s)")
plt.ylabel("Normalized Enso Sublimation Metric (dimensionless)")
plt.title("Normalized Enso Sublimation Metric vs Luminosity")
plt.grid(True, which="both", ls="--", lw=0.5)
plt.tight_layout()
plt.show()

# Save CSV for records
df.to_csv("/mnt/data/enso_sublimation_metric_normalized.csv", index=False)

df
