from pathlib import Path
import pickle
import jayrock

# Target and date_obs
target = jayrock.Target("Ophelia")
target.compute_ephemeris(cycle=6)
date_obs = target.get_date_obs(at="vmag_min")

# Instrument
inst = jayrock.Instrument("MIRI", "MRS")
inst.aperture = "ch2"
inst.disperser = "medium"

inst.detector.nexp = 4
inst.detector.ngroup = 10
inst.detector.nint = 1

# Run simulation
if False:
    obs = jayrock.observe(target, inst, date_obs)

    print(obs.report.keys())
    print(obs.report["scalar"]["sn"])  # snr at reference wavelength
    wave, snr = obs.report["1d"]["sn"]  # snr over wavelength
    print(wave[:5])
    print(snr[:5])

    print(
        obs.partially_saturated
    )  # boolean array over wavelength. True where partially saturated
    print(
        obs.fully_saturated
    )  # boolean array over wavelength. True where fully saturated
    print(
        obs.n_partial_saturated
    )  # number of partially saturated pixels per wavelength
    print(obs.is_saturated)  # any wavelength fully or partially saturated?
    breakpoint()
    # obs.plot_snr()  # plot snr over wavelength
    observations = []

    PATH_CACHE = Path("/tmp/obs_miri_mrs_ophelia.pkl")

    if not PATH_CACHE.is_file():
        for aperture in ["ch1", "ch2", "ch3", "ch4"]:
            for disperser in ["short", "medium", "long"]:
                inst.aperture = aperture
                inst.disperser = disperser

                obs = jayrock.observe(target, inst, date_obs)
                observations.append(obs)
        with open(PATH_CACHE, "wb") as f:
            pickle.dump(observations, f)
    else:
        with open(PATH_CACHE, "rb") as f:
            observations = pickle.load(f)
    jayrock.plot_snr(observations)  # plot SNRs of multiple observations together


observations = []
PATH_CACHE = Path("/tmp/obs_miri_mrs_2_ophelia.pkl")

if not PATH_CACHE.is_file():
    date_thermal_min, date_thermal_max = target.get_date_obs(
        at=["thermal_min", "thermal_max"]
    )
    for date_obs in [date_thermal_min, date_thermal_max]:
        obs = jayrock.observe(target, inst, date_obs)
        observations.append(obs)
    jayrock.plot_snr(observations)  # plot SNRs of multiple observations together
    with open(PATH_CACHE, "wb") as f:
        pickle.dump(observations, f)
else:
    with open(PATH_CACHE, "rb") as f:
        observations = pickle.load(f)
jayrock.plot_snr(observations)  # plot SNRs of multiple observations together
