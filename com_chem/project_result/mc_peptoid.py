import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# =========================
# User settings
# =========================
BASE_DIR = Path(__file__).resolve().parent

CSV_FILES = {
    "SAR": BASE_DIR / "sar_omega.csv",
    "PMP": BASE_DIR / "pmp_omega.csv",
}

TEMPERATURES = [200, 300, 500]   # K
N_STEPS = 50000
BURN_IN = 5000
SAMPLE_EVERY = 10
R_KCAL = 0.0019872041  # kcal/mol/K

# After shifting:
# first angle point -> 0 deg (treated as trans)
TRANS_CENTER = 0
CIS_CENTER = 180
REGION_HALF_WIDTH = 30  # degrees


# =========================
# Helper functions
# =========================
def angular_distance(a, b):
    diff = abs((a - b) % 360)
    return min(diff, 360 - diff)


def classify_angle(angle):
    if angular_distance(angle, TRANS_CENTER) <= REGION_HALF_WIDTH:
        return "trans"
    elif angular_distance(angle, CIS_CENTER) <= REGION_HALF_WIDTH:
        return "cis"
    else:
        return "intermediate"


def summarize_populations(sampled_angles):
    labels = [classify_angle(a) for a in sampled_angles]
    counts = pd.Series(labels).value_counts()
    total = len(labels)
    return {
        "trans": 100 * counts.get("trans", 0) / total,
        "cis": 100 * counts.get("cis", 0) / total,
        "intermediate": 100 * counts.get("intermediate", 0) / total,
    }


def load_scan_data(csv_file):
    """
    Load scan data from CSV.

    Required columns:
      angle_deg
    and one of:
      rel_kcal
      energy_hartree
      energy_kcal

    Special handling:
    - angle does not need to start at 0
    - first angle will be shifted to 0
    - all angles mapped into [0, 360)
    """
    csv_file = Path(csv_file)

    if not csv_file.exists():
        raise FileNotFoundError(f"Cannot find file: {csv_file}")

    df = pd.read_csv(csv_file)
    df.columns = [str(c).strip().lower() for c in df.columns]

    if "angle_deg" not in df.columns:
        raise ValueError(
            f"{csv_file} must contain a column named 'angle_deg'. "
            f"Found columns: {df.columns.tolist()}"
        )

    df["angle_deg"] = pd.to_numeric(df["angle_deg"], errors="coerce")
    df = df.dropna(subset=["angle_deg"]).copy()

    angle0 = df["angle_deg"].iloc[0]
    df["angle_deg"] = (df["angle_deg"] - angle0) % 360

    if "rel_kcal" in df.columns:
        df["rel_kcal"] = pd.to_numeric(df["rel_kcal"], errors="coerce")
    elif "energy_hartree" in df.columns:
        e = pd.to_numeric(df["energy_hartree"], errors="coerce")
        e_min = e.min()
        df["rel_kcal"] = (e - e_min) * 627.509474
    elif "energy_kcal" in df.columns:
        e = pd.to_numeric(df["energy_kcal"], errors="coerce")
        e_min = e.min()
        df["rel_kcal"] = e - e_min
    else:
        raise ValueError(
            f"{csv_file} must contain one of: rel_kcal, energy_hartree, energy_kcal"
        )

    df = df[["angle_deg", "rel_kcal"]].dropna().copy()

    df = df.sort_values(["angle_deg", "rel_kcal"]).drop_duplicates(
        subset=["angle_deg"], keep="first"
    )
    df = df.sort_values("angle_deg").reset_index(drop=True)

    if len(df) < 4:
        raise ValueError(f"{csv_file} has too few valid scan points.")

    return df


# =========================
# Sampling methods
# =========================
def metropolis_mc_neighbor(df, temperature, n_steps=50000, burn_in=5000, sample_every=10, seed=42):
    """
    Neighbor MC: only propose left/right neighboring states.
    Good for showing local trapping/barrier effects.
    """
    rng = np.random.default_rng(seed)

    angles = df["angle_deg"].to_numpy()
    energies = df["rel_kcal"].to_numpy()
    n_states = len(angles)

    current_idx = int(np.argmin(energies))
    current_energy = energies[current_idx]

    sampled_indices = []
    accepted = 0

    for step in range(n_steps):
        move = rng.choice([-1, 1])
        proposal_idx = (current_idx + move) % n_states
        proposal_energy = energies[proposal_idx]

        delta_e = proposal_energy - current_energy

        if delta_e <= 0:
            accept = True
        else:
            prob = np.exp(-delta_e / (R_KCAL * temperature))
            accept = rng.random() < prob

        if accept:
            current_idx = proposal_idx
            current_energy = proposal_energy
            accepted += 1

        if step >= burn_in and ((step - burn_in) % sample_every == 0):
            sampled_indices.append(current_idx)

    sampled_indices = np.array(sampled_indices, dtype=int)

    return {
        "angles": angles[sampled_indices],
        "energies": energies[sampled_indices],
        "acceptance_ratio": accepted / n_steps,
        "all_angles": angles,
        "all_energies": energies,
    }


def metropolis_mc_global(df, temperature, n_steps=50000, burn_in=5000, sample_every=10, seed=42):
    """
    Global MC: propose jumps to any angle state.
    Better approximation to equilibrium sampling.
    """
    rng = np.random.default_rng(seed)

    angles = df["angle_deg"].to_numpy()
    energies = df["rel_kcal"].to_numpy()
    n_states = len(angles)

    current_idx = int(np.argmin(energies))
    current_energy = energies[current_idx]

    sampled_indices = []
    accepted = 0

    for step in range(n_steps):
        proposal_idx = rng.integers(0, n_states)
        while proposal_idx == current_idx:
            proposal_idx = rng.integers(0, n_states)

        proposal_energy = energies[proposal_idx]
        delta_e = proposal_energy - current_energy

        if delta_e <= 0:
            accept = True
        else:
            prob = np.exp(-delta_e / (R_KCAL * temperature))
            accept = rng.random() < prob

        if accept:
            current_idx = proposal_idx
            current_energy = proposal_energy
            accepted += 1

        if step >= burn_in and ((step - burn_in) % sample_every == 0):
            sampled_indices.append(current_idx)

    sampled_indices = np.array(sampled_indices, dtype=int)

    return {
        "angles": angles[sampled_indices],
        "energies": energies[sampled_indices],
        "acceptance_ratio": accepted / n_steps,
        "all_angles": angles,
        "all_energies": energies,
    }


def boltzmann_populations(df, temperature):
    """
    Exact equilibrium populations from Boltzmann weighting on discrete scan points.
    """
    energies = df["rel_kcal"].to_numpy()
    angles = df["angle_deg"].to_numpy()

    weights = np.exp(-energies / (R_KCAL * temperature))
    probs = weights / weights.sum()

    trans_prob = 0.0
    cis_prob = 0.0
    intermediate_prob = 0.0

    for angle, p in zip(angles, probs):
        label = classify_angle(angle)
        if label == "trans":
            trans_prob += p
        elif label == "cis":
            cis_prob += p
        else:
            intermediate_prob += p

    return {
        "trans": trans_prob * 100,
        "cis": cis_prob * 100,
        "intermediate": intermediate_prob * 100,
    }


# =========================
# Plot functions
# =========================
def plot_energy_profiles(data_dict, output="energy_profiles_comparison.png"):
    plt.figure(figsize=(8, 5))
    for label, df in data_dict.items():
        plt.plot(df["angle_deg"], df["rel_kcal"], marker="o", label=label)

    plt.xlabel(r"Shifted $\omega$ angle (deg)")
    plt.ylabel("Relative energy (kcal/mol)")
    plt.title("Torsional energy profiles")
    plt.xticks(np.arange(0, 361, 60))
    plt.legend()
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def plot_population_method_comparison(pop_neighbor, pop_global, pop_boltz, system_label, output):
    categories = ["trans", "cis", "intermediate"]
    x = np.arange(len(categories))
    width = 0.22

    plt.figure(figsize=(8, 5))
    vals1 = [pop_neighbor[c] for c in categories]
    vals2 = [pop_global[c] for c in categories]
    vals3 = [pop_boltz[c] for c in categories]

    plt.bar(x - width, vals1, width=width, label="Neighbor MC")
    plt.bar(x,         vals2, width=width, label="Global MC")
    plt.bar(x + width, vals3, width=width, label="Boltzmann")

    plt.xticks(x, categories)
    plt.ylabel("Population (%)")
    plt.title(f"{system_label}: method comparison")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


# =========================
# Main
# =========================
def main():
    data_dict = {}
    neighbor_results = {}
    global_results = {}
    boltz_results = {}

    for label, csv_file in CSV_FILES.items():
        data_dict[label] = load_scan_data(csv_file)

    plot_energy_profiles(data_dict)

    rows = []

    print("\n=== Population summary ===")
    for sys_idx, (label, df) in enumerate(data_dict.items()):
        neighbor_results[label] = {}
        global_results[label] = {}
        boltz_results[label] = {}

        print(f"\nSystem: {label}")
        for temp_idx, temp in enumerate(TEMPERATURES):
            neigh = metropolis_mc_neighbor(
                df,
                temperature=temp,
                n_steps=N_STEPS,
                burn_in=BURN_IN,
                sample_every=SAMPLE_EVERY,
                seed=1000 + sys_idx * 100 + temp_idx,
            )
            glob = metropolis_mc_global(
                df,
                temperature=temp,
                n_steps=N_STEPS,
                burn_in=BURN_IN,
                sample_every=SAMPLE_EVERY,
                seed=2000 + sys_idx * 100 + temp_idx,
            )
            boltz = boltzmann_populations(df, temp)

            neighbor_results[label][temp] = summarize_populations(neigh["angles"])
            global_results[label][temp] = summarize_populations(glob["angles"])
            boltz_results[label][temp] = boltz

            print(f"  Temperature: {temp} K")
            print(f"    Neighbor MC (acc={neigh['acceptance_ratio']:.3f}): "
                  f"trans={neighbor_results[label][temp]['trans']:.2f}%, "
                  f"cis={neighbor_results[label][temp]['cis']:.2f}%, "
                  f"intermediate={neighbor_results[label][temp]['intermediate']:.2f}%")
            print(f"    Global MC   (acc={glob['acceptance_ratio']:.3f}): "
                  f"trans={global_results[label][temp]['trans']:.2f}%, "
                  f"cis={global_results[label][temp]['cis']:.2f}%, "
                  f"intermediate={global_results[label][temp]['intermediate']:.2f}%")
            print(f"    Boltzmann: "
                  f"trans={boltz['trans']:.2f}%, "
                  f"cis={boltz['cis']:.2f}%, "
                  f"intermediate={boltz['intermediate']:.2f}%")

            rows.append({
                "system": label,
                "temperature_K": temp,
                "method": "neighbor_mc",
                "trans_percent": neighbor_results[label][temp]["trans"],
                "cis_percent": neighbor_results[label][temp]["cis"],
                "intermediate_percent": neighbor_results[label][temp]["intermediate"],
            })
            rows.append({
                "system": label,
                "temperature_K": temp,
                "method": "global_mc",
                "trans_percent": global_results[label][temp]["trans"],
                "cis_percent": global_results[label][temp]["cis"],
                "intermediate_percent": global_results[label][temp]["intermediate"],
            })
            rows.append({
                "system": label,
                "temperature_K": temp,
                "method": "boltzmann",
                "trans_percent": boltz["trans"],
                "cis_percent": boltz["cis"],
                "intermediate_percent": boltz["intermediate"],
            })

            plot_population_method_comparison(
                neighbor_results[label][temp],
                global_results[label][temp],
                boltz,
                system_label=f"{label} at {temp} K",
                output=f"{label.lower()}_{temp}K_method_comparison.png",
            )

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv("population_method_comparison.csv", index=False)

    print("\nFiles written:")
    print("- energy_profiles_comparison.png")
    print("- population_method_comparison.csv")
    for label in CSV_FILES.keys():
        for temp in TEMPERATURES:
            print(f"- {label.lower()}_{temp}K_method_comparison.png")


if __name__ == "__main__":
    main()