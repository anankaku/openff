import numpy as np
import pandas as pd
from pathlib import Path

output_dir = Path("processed_dft")
output_dir.mkdir(exist_ok=True)

temperature = 300.0
hartree_to_kcal = 627.509474
R = 0.0019872041


def wrap_angle(angle):
    return ((angle + 180.0) % 360.0) - 180.0


def read_scan_table(filename):
    data = []

    with open(filename, "r", errors="ignore") as f:
        for line in f:
            parts = line.split()

            if len(parts) == 2:
                try:
                    phi_raw = float(parts[0])
                    energy = float(parts[1])
                    data.append((phi_raw, energy))
                except ValueError:
                    pass

    df = pd.DataFrame(data, columns=["phi_raw", "energy_hartree"])

    if df.empty:
        raise ValueError(f"No scan data found in {filename}")

    return df


scan_map = pd.read_csv("scan.csv")

all_scans = []

for _, row in scan_map.iterrows():
    filename = row["file"]
    psi_raw = float(row["psi_deg"])
    psi_deg = wrap_angle(psi_raw)

    df = read_scan_table(filename)

    # 37 points: last point closes the cycle, so keep first 36
    df = df.iloc[:36].copy()

    df["step"] = np.arange(len(df))
    df["phi_deg"] = df["phi_raw"].apply(wrap_angle)
    df["psi_raw"] = psi_raw
    df["psi_deg"] = psi_deg
    df["source_file"] = filename

    out_1d = output_dir / f"{Path(filename).stem}_processed.csv"
    df.to_csv(out_1d, index=False)

    all_scans.append(df)

df_all = pd.concat(all_scans, ignore_index=True)

df_all = df_all.drop_duplicates(
    subset=["phi_deg", "psi_deg"],
    keep="first"
).reset_index(drop=True)

df_all["deltaE_kcal"] = (
    df_all["energy_hartree"] - df_all["energy_hartree"].min()
) * hartree_to_kcal

df_all["boltzmann_weight"] = np.exp(
    -df_all["deltaE_kcal"] / (R * temperature)
)

df_all["population"] = (
    df_all["boltzmann_weight"] / df_all["boltzmann_weight"].sum()
)

df_all.to_csv(output_dir / "dft_2d_phi_psi_population.csv", index=False)

print("Done.")
print("Number of input scans:", len(scan_map))
print("Number of final points:", len(df_all))
print("Expected if 36 scans x 36 points:", 36 * 36)
print(df_all.head())