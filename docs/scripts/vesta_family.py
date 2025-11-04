# ------
# Sample Selection
import rocks

# Load table of asteroid properties (SsODNet/BFT)
bft = rocks.load_bft()

# Select candidate targets: Vesta family members with known diameter
is_vestoid = bft["family.family_name"] == "Vesta"
has_diameter = bft["diameter.value"].notnull()

vestoids = bft[is_vestoid & has_diameter].copy().reset_index(drop=True)

print(f"Vestoids with known diameter: {len(vestoids)}")  # 2041

# ------
# Instrument definition
import jayrock

miri_mrs = jayrock.Instrument("MIRI", mode="MRS")

miri_mrs.detector.nexp = 4  # 4-pt dither
miri_mrs.detector.readout_pattern = "fastr1"  # recommended for MRS

miri_mrs.aperture = "ch2"
miri_mrs.disperser = "long"

miri_lrs = jayrock.Instrument("MIRI", mode="LRSSLIT")

miri_lrs.detector.nexp = 2  # 2-pt dither
miri_lrs.detector.readout_pattern = "fastr1"  # recommended for MRS


# ------
# Get largest Vestoid within diameter range
import pandas as pd

diameter_bins = [1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 96, 128, 256, 512, 1024]

vestoids["diameter_bin"] = pd.cut(
    vestoids["diameter.value"], bins=diameter_bins, right=True, include_lowest=True
)

vestoids_largest_per_diameter_bin = vestoids.loc[
    vestoids.groupby("diameter_bin", observed=True)["diameter.value"].idxmax()
]
print(vestoids_largest_per_diameter_bin[["sso_name", "diameter.value", "diameter_bin"]])


# ------
# Estimate texp to reach SNR 300-400
from pathlib import Path

PATH_CACHE = Path("candidate_vestoides.csv")

if not PATH_CACHE.is_file():
    for idx, asteroid in vestoids_largest_per_diameter_bin.iterrows():
        # Define target and compute ephemeris for JWST Cycle 6
        target = jayrock.Target(asteroid["sso_name"])
        target.compute_ephemeris(cycle=6)
        target.print_ephemeris()

        # Add total length of visibility windows to dataframe
        vestoids.loc[idx, "N_days_observable"] = len(target.ephemeris)

        # Get date of minimum thermal flux during visibility window
        date_obs = target.get_date_obs(at="thermal_min")

        # Compute exposure times for SNR targets
        for inst in [miri_mrs, miri_lrs]:
            # Estimate ngroup and nint
            success = inst.set_snr_target([300, 400], target, date_obs)

            if not success:
                print(f" Could not estimate SNR, skipping...")
                continue

            # Record results
            vestoids.loc[idx, f"nint_{inst.mode}"] = inst.detector.nint
            vestoids.loc[idx, f"ngroup_{inst.mode}"] = inst.detector.ngroup
            vestoids.loc[idx, f"texp_{inst.mode}"] = inst.texp
            vestoids.loc[idx, f"snr_{inst.mode}"] = inst.estimated_snr

    columns = [
        "sso_number",
        "sso_name",
        "diameter.value",
        "N_days_observable",
        "nint_lrsslit",
        "ngroup_lrsslit",
        "texp_lrsslit",
        "snr_lrsslit",
        "nint_mrs",
        "ngroup_mrs",
        "texp_mrs",
        "snr_mrs",
    ]

    possible_targets = vestoids.loc[vestoids.N_days_observable > 0, columns].copy()
    possible_targets.to_csv(PATH_CACHE, index=False)
else:
    vestoids = pd.read_csv(PATH_CACHE)


print(
    vestoids[["sso_name", "diameter.value", "texp_mrs", "texp_lrsslit"]].sort_values(
        "diameter.value"
    )
)

# ------
# Plot results
if False:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    ax.scatter(
        vestoids["diameter.value"],
        vestoids["texp_mrs"] / 60,
        label="MIRI MRS",
        s=50,
    )

    ax.scatter(
        vestoids["diameter.value"],
        vestoids["texp_lrsslit"] / 60,
        label="MIRI LRS",
        s=50,
        marker="s",
    )

    ax.set(
        xscale="log", yscale="log", xlabel="Diameter / km", ylabel="Exposure time / min"
    )
    ax.legend()
    plt.show()

TARGETS_MRS = ["Kollaa", "Koskenniemi", "Mila", "Robelmonte", "Ausonia"]

SNR_targets = {
    "ch2": {"short": 300, "medium": 300, "long": 300},
    "ch3": {"short": 300, "medium": 300, "long": 300},
    "ch4": {"short": 100, "medium": 100, "long": 100},
}

if False:
    MRS_OBSERVING_PARAMETERS = {}

    for target in TARGETS_MRS:
        # Define target and compute ephemeris for JWST Cycle 6
        target = jayrock.Target(target)
        target.compute_ephemeris(cycle=6)

        # Get date of minimum thermal flux during visibility window
        date_obs = target.get_date_obs(at="thermal_min")

        for channel in ["ch1", "ch2", "ch3", "ch4"]:
            for disperser in ["short", "medium", "long"]:
                # Set aperture and disperser
                miri_mrs.aperture = channel
                miri_mrs.disperser = disperser

                if channel == "ch1":
                    # Use minimum ngroup/nint
                    miri_mrs.detector.ngroup = 5
                    miri_mrs.detector.nint = 1
                else:
                    # Get SNR target
                    snr_target = SNR_targets[channel][disperser]
                    # Estimate ngroup and nint based on SNR target
                    miri_mrs.set_snr_target(snr_target, target, date_obs)

                print(
                    f"{target.name:15s} | {channel:3s} | {disperser:6s} | "
                    f"ngroup: {miri_mrs.detector.ngroup:2d} | nint: {miri_mrs.detector.nint:2d} | "
                    f"texp: {miri_mrs.texp / 60:6.2f} min"
                )

                if target.name not in MRS_OBSERVING_PARAMETERS:
                    MRS_OBSERVING_PARAMETERS[target.name] = {}

                MRS_OBSERVING_PARAMETERS[target.name][(channel, disperser)] = {
                    "ngroup": miri_mrs.detector.ngroup,
                    "nint": miri_mrs.detector.nint,
                }
            print(MRS_OBSERVING_PARAMETERS)
else:
    MRS_OBSERVING_PARAMETERS = {
        "Kollaa": {
            ("ch1", "short"): {"ngroup": 5, "nint": 1},
            ("ch1", "medium"): {"ngroup": 5, "nint": 1},
            ("ch1", "long"): {"ngroup": 5, "nint": 1},
            ("ch2", "short"): {"ngroup": 83, "nint": 4},
            ("ch2", "medium"): {"ngroup": 69, "nint": 2},
            ("ch2", "long"): {"ngroup": 74, "nint": 1},
            ("ch3", "short"): {"ngroup": 62, "nint": 1},
            ("ch3", "medium"): {"ngroup": 47, "nint": 1},
            ("ch3", "long"): {"ngroup": 40, "nint": 1},
            ("ch4", "short"): {"ngroup": 28, "nint": 1},
            ("ch4", "medium"): {"ngroup": 48, "nint": 1},
            ("ch4", "long"): {"ngroup": 91, "nint": 3},
        },
        "Koskenniemi": {
            ("ch1", "short"): {"ngroup": 5, "nint": 1},
            ("ch1", "medium"): {"ngroup": 5, "nint": 1},
            ("ch1", "long"): {"ngroup": 5, "nint": 1},
            ("ch2", "short"): {"ngroup": 86, "nint": 2},
            ("ch2", "medium"): {"ngroup": 72, "nint": 1},
            ("ch2", "long"): {"ngroup": 41, "nint": 1},
            ("ch3", "short"): {"ngroup": 35, "nint": 1},
            ("ch3", "medium"): {"ngroup": 27, "nint": 1},
            ("ch3", "long"): {"ngroup": 24, "nint": 1},
            ("ch4", "short"): {"ngroup": 18, "nint": 1},
            ("ch4", "medium"): {"ngroup": 30, "nint": 1},
            ("ch4", "long"): {"ngroup": 69, "nint": 2},
        },
        "Mila": {
            ("ch1", "short"): {"ngroup": 5, "nint": 1},
            ("ch1", "medium"): {"ngroup": 5, "nint": 1},
            ("ch1", "long"): {"ngroup": 5, "nint": 1},
            ("ch2", "short"): {"ngroup": 88, "nint": 2},
            ("ch2", "medium"): {"ngroup": 71, "nint": 1},
            ("ch2", "long"): {"ngroup": 40, "nint": 1},
            ("ch3", "short"): {"ngroup": 33, "nint": 1},
            ("ch3", "medium"): {"ngroup": 26, "nint": 1},
            ("ch3", "long"): {"ngroup": 22, "nint": 1},
            ("ch4", "short"): {"ngroup": 17, "nint": 1},
            ("ch4", "medium"): {"ngroup": 27, "nint": 1},
            ("ch4", "long"): {"ngroup": 62, "nint": 2},
        },
        "Robelmonte": {
            ("ch1", "short"): {"ngroup": 5, "nint": 1},
            ("ch1", "medium"): {"ngroup": 5, "nint": 1},
            ("ch1", "long"): {"ngroup": 5, "nint": 1},
            ("ch2", "short"): {"ngroup": 15, "nint": 1},
            ("ch2", "medium"): {"ngroup": 9, "nint": 1},
            ("ch2", "long"): {"ngroup": 7, "nint": 1},
            ("ch3", "short"): {"ngroup": 7, "nint": 1},
            ("ch3", "medium"): {"ngroup": 6, "nint": 1},
            ("ch3", "long"): {"ngroup": 6, "nint": 1},
            ("ch4", "short"): {"ngroup": 6, "nint": 1},
            ("ch4", "medium"): {"ngroup": 8, "nint": 1},
            ("ch4", "long"): {"ngroup": 20, "nint": 1},
        },
        "Ausonia": {
            ("ch1", "short"): {"ngroup": 5, "nint": 1},
            ("ch1", "medium"): {"ngroup": 5, "nint": 1},
            ("ch1", "long"): {"ngroup": 5, "nint": 1},
            ("ch2", "short"): {"ngroup": 6, "nint": 1},
            ("ch2", "medium"): {"ngroup": 5, "nint": 1},
            ("ch2", "long"): {"ngroup": 5, "nint": 1},
            ("ch3", "short"): {"ngroup": 5, "nint": 1},
            ("ch3", "medium"): {"ngroup": 5, "nint": 1},
            ("ch3", "long"): {"ngroup": 5, "nint": 1},
            ("ch4", "short"): {"ngroup": 5, "nint": 1},
            ("ch4", "medium"): {"ngroup": 5, "nint": 1},
            ("ch4", "long"): {"ngroup": 6, "nint": 1},
        },
    }

for target, params in MRS_OBSERVING_PARAMETERS.items():
    continue
    observations = []

    print(f"Target: {target}")
    target = jayrock.Target(target)
    target.compute_ephemeris(cycle=6)

    date_obs_thermal_min = target.get_date_obs(at="thermal_min")
    date_obs_thermal_max = target.get_date_obs(at="thermal_max")

    for (channel, disperser), settings in params.items():
        ngroup = settings["ngroup"]
        nint = settings["nint"]

        miri_mrs.aperture = channel
        miri_mrs.disperser = disperser
        miri_mrs.detector.ngroup = ngroup
        miri_mrs.detector.nint = nint

        obs_thermal_min = jayrock.observe(target, miri_mrs, date_obs_thermal_min)
        obs_thermal_max = jayrock.observe(target, miri_mrs, date_obs_thermal_max)
        observations.append(obs_thermal_min)
        observations.append(obs_thermal_max)

    # Save plot to file
    jayrock.plot_snr(
        observations, show=False, save_to=f"{target.name}_miri_mrs_snr.png"
    )

TARGETS_LRS = ["2000 WE8", "2000 QA163"]
SNR_TARGET = 300

if False:
    LRS_OBSERVING_PARAMETERS = {}

    for target in TARGETS_LRS:
        # Define target and compute ephemeris for JWST Cycle 6
        target = jayrock.Target(target)
        target.compute_ephemeris(cycle=6)

        # Get date of minimum thermal flux during visibility window
        date_obs = target.get_date_obs(at="thermal_min")

        # Compute exposure times for SNR target at 10 micron
        miri_lrs.set_snr_target(SNR_TARGET, target, date_obs, wave=10)

        print(
            f"{target.name:15s} | "
            f"ngroup: {miri_lrs.detector.ngroup:2d} | nint: {miri_lrs.detector.nint:2d} | "
            f"texp: {miri_lrs.texp / 60:6.2f} min"
        )

        LRS_OBSERVING_PARAMETERS[target.name] = {
            "ngroup": miri_lrs.detector.ngroup,
            "nint": miri_lrs.detector.nint,
        }
    print(LRS_OBSERVING_PARAMETERS)
else:
    LRS_OBSERVING_PARAMETERS = {
        "2000 WE8": {"ngroup": 65, "nint": 2},
        "2000 QA163": {"ngroup": 36, "nint": 1},
    }


for target, params in LRS_OBSERVING_PARAMETERS.items():
    continue
    observations = []

    print(f"Target: {target}")
    target = jayrock.Target(target)
    target.compute_ephemeris(cycle=6)

    date_obs_thermal_min = target.get_date_obs(at="thermal_min")
    date_obs_thermal_max = target.get_date_obs(at="thermal_max")

    miri_lrs.detector.ngroup = params["ngroup"]
    miri_lrs.detector.nint = params["nint"]

    obs_thermal_min = jayrock.observe(target, miri_lrs, date_obs_thermal_min)
    obs_thermal_max = jayrock.observe(target, miri_lrs, date_obs_thermal_max)
    observations.append(obs_thermal_min)
    observations.append(obs_thermal_max)

    # Save plot to file
    jayrock.plot_snr(
        observations, show=False, save_to=f"{target.name}_miri_lrs_snr.png"
    )

TARGETS = TARGETS_MRS + TARGETS_LRS

for target in TARGETS:
    target = jayrock.Target(target)
    target.compute_ephemeris(cycle=6)

    date_obs = target.get_date_obs(at="thermal_min")
    target.export_spectrum(
        f"{target.name}_spectrum_thermal_minimum.txt", date_obs=date_obs
    )
