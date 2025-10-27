# Calculate Enso Sublimation Metric (ESM) for Sagittarius A* using PMZ=0.0102
# Assumptions:
# - Dust sublimation radius scaling: R_sub(pc) ≈ 0.5 * sqrt(L/1e45 erg/s)
# - Sgr A* mass ~ 4.297e6 Msun
# - Two luminosity states: quiescent L=1e36 erg/s, bright flare L=1.5e39 erg/s
# - Eddington luminosity: L_Edd = 1.26e38 * (M/Msun) erg/s
#
# ESM definition (proposed):
#   ESM = PMZ * R_sub(AU)
# where R_sub(AU) = 0.5 * 206265 * sqrt(L/1e45)

import math
import pandas as pd
from caas_jupyter_tools import display_dataframe_to_user

PMZ = 0.0102
M_sun = 1.98847e30  # kg
G = 6.67430e-11     # m^3/kg/s^2
c = 299792458       # m/s
AU_m = 1.495978707e11  # m
pc_to_AU = 206265.0

M = 4.297e6  # solar masses
L_quiescent = 1e36  # erg/s
L_flare = 1.5e39    # erg/s

def R_sub_AU(L_erg_s):
    return 0.5 * pc_to_AU * math.sqrt(L_erg_s / 1e45)

def L_edd_erg_s(M_solar):
    return 1.26e38 * M_solar

def schwarzschild_radius_AU(M_solar):
    # Rs = 2GM/c^2
    M_kg = M_solar * M_sun
    Rs_m = 2 * G * M_kg / (c**2)
    return Rs_m / AU_m

rows = []
for label, L in [("Quiescent", L_quiescent), ("Bright flare", L_flare)]:
    Rsub = R_sub_AU(L)
    Ledd = L_edd_erg_s(M)
    edd_ratio = L / Ledd
    ESM = PMZ * Rsub
    rows.append({
        "State": label,
        "Luminosity L (erg/s)": L,
        "Eddington ratio (L/L_Edd)": edd_ratio,
        "Dust sublimation radius (AU)": Rsub,
        "ESM = PMZ * R_sub (dimensionless, AU-scaled)": ESM
    })

Rs_AU = schwarzschild_radius_AU(M)

df = pd.DataFrame(rows)
df_meta = pd.DataFrame([{"Parameter": "Sagittarius A* mass (Msun)", "Value": M},
                        {"Parameter": "Schwarzschild radius (AU)", "Value": Rs_AU},
                        {"Parameter": "PMZ factor", "Value": PMZ}])

display_dataframe_to_user("Enso Sublimation Metric for Sagittarius A*", df)
display_dataframe_to_user("Constants & Derived Parameters", df_meta)

# Save as CSVs for the user to download if desired
df.to_csv("/mnt/data/enso_sublimation_metric_sgra.csv", index=False)
df_meta.to_csv("/mnt/data/enso_sublimation_metric_constants.csv", index=False)

print("Files saved:\n - /mnt/data/enso_sublimation_metric_sgra.csv\n - /mnt/data/enso_sublimation_metric_constants.csv")
